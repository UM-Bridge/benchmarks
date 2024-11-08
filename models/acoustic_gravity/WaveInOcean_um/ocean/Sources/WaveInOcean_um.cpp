#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include <complex>

#include "OndoMathX.h"
#include "Utility.h" 
#include "Parameters_um.h" 

using namespace OndoMathX;

//////////// For um-bridge 
#include <chrono>
#include <thread>

// Needed for HTTPS, implies the need for openssl, may be omitted if HTTP suffices
// #define CPPHTTPLIB_OPENSSL_SUPPORT

#include "umbridge.h"

class ExampleModel : public umbridge::Model {
public:

  ExampleModel(int test_delay)
   : umbridge::Model("forward"),
     test_delay(test_delay)
  {}

  // Define input and output dimensions of model (here we have a single vector of length 1 for input; same for output)
  std::vector<std::size_t> GetInputSizes(const json& config_json) const override {
    return {1};
  }

  std::vector<std::size_t> GetOutputSizes(const json& config_json) const override {
    return {600};
  }

  std::vector<std::vector<double>> Evaluate(const std::vector<std::vector<double>>& inputs, json config) override {
    // Do the actual model evaluation; here we just multiply the first entry of the first input vector by two, and store the result in the output.
    // In addition, we support an artificial delay here, simulating actual work being done.
    // std::this_thread::sleep_for(std::chrono::milliseconds(test_delay));
    //return {{inputs[0][0] * 2.0}};

    //assert(argc==3);
    // TODO: replace argv by the json config file
    //std::string FileName = argv[2]+NameSimu;
    std::string FileName = "/Code/ocean/Data/"+NameSimu; 
    ofstream infoFile(FileName+"_infos.txt");

    ////////////////////////////////////
    //std::cout << "Running WaveInOcean........................" << std::endl;
    //std::cout << "Simulation of acoustic and gravity waves in an ocean" << std::endl;
    
    Scaling<2> scaling(Lx,Lz);
    Square<Scaling<2>> fespace(FEOrderX,FEOrderZ,Nx,Nz,scaling);
    Center[0] = Lx*Center[0];
    Center[1] = Lz*Center[1];
    //std::cout<<"Center at "<< Center[0] << " " << Center[1]<< endl;
    infoFile<<"Center at "<< Center[0] << " " << Center[1]<< endl;
    Index NDoFs = fespace.getNumDoFs();

    ////////// Get boundary spaces
    auto bottom_fespace = fespace.getBndyFESpace(0);
    auto surface_fespace = fespace.getBndyFESpace(2);
    
    ///// Write outputs for an array of points //////////////////
    ofstream obsPoint_P(FileName+"_P.txt");
    std::vector<Index> iObsPoint_vector; 
    for (Index i = 0; i < nPointX; i++){ 
        for (Index j = 0; j < nPointZ; j++){
            RealVector q{listX[i]/Lx, listZ[j]/Lz};
            Index iObsPoint = findGlobLocIndex(q, fespace);
            // If iObsPoint the given point is really in the domain and we store it
            if(iObsPoint < fespace.getNumDoFs())
            {   
                iObsPoint_vector.push_back(iObsPoint);
                RealVector realQ;
                fespace.getDoFCoordinate(iObsPoint,realQ);
                obsPoint_P << "(" << realQ[0] << ";" << realQ[1] << "):" << iObsPoint << " , " ; 
            }
        }
    }
    // Writes the number of different X and Z coordinates for plots later
    obsPoint_P << nPointX << " , " << nPointZ << " , " ;
    // Writes the number of points for reading the output file later
    obsPoint_P << std::size(iObsPoint_vector) << std::endl ; 
    /////////////////////////////////////////////////////////////////

    ///////////  Construct the parameter fields
    Field::Scalar<Field::Regularity::C0> Rho(_Rho);
    Field::Scalar<Field::Regularity::C0> Rho_cm2(_Rho_cm2);
    Field::ScalarToVector<2, Field::Regularity::C0> scalar_to_vect_ez(_field_ez);
    ///////////////////////////////////////////

    ///////////  Construct the matrices
    // The mass matrices are defined with weighted integration, the weight can be a scalar field or a constant 
    // e.g. CreateMass<1>(fespace, Rho_c2, massMatrixPhi_Sparse, 1.0); second parameter is the scalar field, last parameter is the constant. 
    Field::Scalar unit_field(1.0); 

    //// Mass matrix for phi 
    LAL::DiagonalMatrix massMatrixPhi;
    CreateMass<1>(fespace, Rho_cm2, massMatrixPhi, 1.0); 
    Real rhoSurface = _Rho({0,1,0});
    AssembleMassBndy(surface_fespace, Field::IdentityField, massMatrixPhi, 1./delta2*rhoSurface); 
 
    // Matrices for the bottom source term
    // The mass matrix is created as a sparse matrix because the function CreateMassBndy is defined for a sparse matrix
    LAL::SparseMatrix bottomMassMatrix;

    Real rhoBottom = _Rho({0,0,0});
    CreateMassBndy(bottom_fespace, unit_field, bottomMassMatrix, rhoBottom);

    // Matrices for the velocity
    LAL::DiagonalMatrix massMatrixVelocity; 
    CreateMass<2>(fespace,unit_field,massMatrixVelocity,1.0);

    // Matrices for the pressure
    LAL::DiagonalMatrix massMatrixUnit, massMatrixRho; 
    CreateMass<1>(fespace,unit_field,massMatrixUnit,1.0);
    CreateMass<1>(fespace,Rho,massMatrixRho,1.0);
    ////////////////////////////////////////


    ////////////////////////////////////////
    // Computation of the optimal time step
    Index N_iter;
    Real Residue;
    Real rho;

    rho = LAL::PowerIteration(  [&massMatrixPhi](LAL::Vector &v) {LAL::Solve(massMatrixPhi,v);},
                                [&fespace, &Rho]      (LAL::Vector &u,LAL::Vector &v)
                               { 
                                MltAddStiffness(fespace,Rho,1,u,0,v);
                                },
                               NDoFs,N_iter,Residue);
    Real dt =  0.95 * (2.0 / sqrt(rho));

    //std::cout << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;
    infoFile << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;
   
    //std::cout << "\t" << dt << std::endl;
    infoFile << "\t" << dt << std::endl;
    ////////////////////////////////////////
    
    ////////////////////////////////////////
    // Vectors for the solutions
    LAL::VectorSequence phi(3);
    phi.Allocate(NDoFs);

    LAL::Vector velocity;
    LAL::Allocate(velocity, 2*NDoFs);

    LAL::VectorSequence pressure(2);
    pressure.Allocate(NDoFs);
    ////////////////////////////////////////

    ////////////////////////////////////////
    // Vectors to store the source term
    LAL::Vector sourceVol;
    LAL::Allocate(sourceVol, NDoFs);

    // Source on the seabed
    LAL::Vector sourceCoord, source;
    LAL::Allocate(sourceCoord, bottom_fespace.getNumDoFs());
    LAL::Allocate(source, bottom_fespace.getNumDoFs());
    LAL::Allocate(sourceVol, NDoFs);

    // Compute f(x) at each DoF of the seabed 
    bottom_fespace.setAsFESpace();
    for (Index i=0;i<bottom_fespace.getNumDoFs();++i)
    {
        Index dof = bottom_fespace._glob_bdny2glob_sq(i);
        RealVector P;
        // Get the (x,z) coordinates to evaluate the function sourceSpace
        fespace.getDoFCoordinate(dof,P);
        sourceCoord(i) = functionSpace(P) ;
    }
    // Get FE approximation of f(x) 
    source = bottomMassMatrix * sourceCoord;
    // Store the result in a vector of the size of the whole domain (NDoFs entries) to simplify the sum 
    for (Index i=0;i<bottom_fespace.getNumDoFs();++i)
    {
        Index dof = bottom_fespace._glob_bdny2glob_sq(i);
        sourceVol(dof) = source(i) ;
    }
    ////////////////////////////////////////

    ///////////////// For PML /////////////
    LAL::VectorSequence phi_hat(2);
    phi_hat.Allocate(NDoFs);

    LAL::VectorSequence wx(2);
    wx.Allocate(NDoFs);

    auto subspace_pml = fespace.getSubDomainFESpace([](RealVector xyz){
            if ((xyz[0] > Lx-LxPml_R) || (xyz[0] < LxPml_L))
                return true;
            else return false;
    });
    auto fespace_pml = subspace_pml;
    subspace_pml.setAsSubSpace(); 
    fespace_pml.setAsFESpace(); 

    //Field::ScalarToVector<2, Field::Regularity::C0>  vect_ex(ex);
    Field::Vector cst_vec_ez(ez);
    Field::Vector cst_vec_ex(ex);

    std::map<Index, Real> Dissip_map;
    Dissip_map[0] = 0.0;
    Dissip_map[1] = Dissip_coeff;
    Field::Scalar<Field::Regularity::PiecewiseConstant> Dissip_Disc_Field(Dissip_map);

    LAL::DiagonalMatrix massMatrixPhi_Dissip, massMatrixPhi_Bndy_Dissip;
    CreateMass<1>(subspace_pml, Rho_cm2, massMatrixPhi_Dissip, Dissip_coeff);
    CreateMassBndy(surface_fespace, Dissip_Disc_Field, massMatrixPhi_Bndy_Dissip, 1./delta2*rhoSurface);

    LAL::DiagonalMatrix massMatrix_p;
    LAL::DiagonalMatrix massMatrix_m; 
    massMatrix_p = massMatrixPhi + 0.5*dt*massMatrixPhi_Dissip + 0.5*dt*massMatrixPhi_Bndy_Dissip;
    massMatrix_m = massMatrixPhi - 0.5*dt*massMatrixPhi_Dissip - 0.5*dt*massMatrixPhi_Bndy_Dissip;
        
    /// Mass matrix for wx
    LAL::DiagonalMatrix massMatrix_wx;
    CreateMassDG(fespace_pml, Field::IdentityField, massMatrix_wx,1.0);
        
    //////////////////////////////////////

    ////////////////////////////////////////
    // Time and output 
    Index NOutput = 0;
    Index OutputFreq = floor(OutputDeltat/dt); 
    
    //Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;
    for  (Index iStep = 0; iStep*dt < T_end; iStep++)
    {
        // Operation  K_phi
        Core::MltAdd< Square<Scaling<2>>,Square<Scaling<2>>,true,true,true, Core::QuadratureOpt::Default,
            1,       //DimU
            1,       //DimV
            2,       //DimI
            1,       //DimJ
            Field::Scalar<Field::Regularity::C0>,
            Core::DiffOp::Gradient,
            Core::DiffOp::Gradient
            >
            (fespace,fespace,Rho,-dt*dt,phi.getVector(1),0.0,phi.getVector(0)); 

        // Add bottom source
        Real time = functionTime(dt*iStep);
        phi.getVector(0) += dt*dt*sourceVol*time;
        
        // Without PML contribution: 
        // Solve Phi here 
        //LAL::Solve(massMatrixPhi,phi.getVector(0));
        //phi.getVector(0) +=  2.0*phi.getVector(1);
        //phi.getVector(0) += -1.0*phi.getVector(2);
        // And ignore everything from there to pressure computation

        // PML contributions 
        // Contribution of phi_hat
        Core::MltAdd<decltype(fespace_pml),decltype(subspace_pml),false,true,true, 
            Core::QuadratureOpt::Default,
            1,       //DimU
            1,       //DimV
            1,       //DimI
            1,       //DimJ
            Field::Scalar<Field::Regularity::C0>,
            Core::DiffOp::Gradient,
            Core::DiffOp::Gradient, 
            Field::Vector<2>,
            Field::Vector<2>
            >
            (fespace_pml, subspace_pml, Rho, -dt*dt, phi_hat.getVector(1), 1.0, phi.getVector(0),
                cst_vec_ez, cst_vec_ez); 
        
        // Contribution of wx
        Core::MltAdd<decltype(fespace_pml),decltype(subspace_pml),false,true,true, 
            Core::QuadratureOpt::Default,
            1,       //DimU
            1,       //DimV
            2,       //DimI
            1,       //DimJ
            Field::Scalar<Field::Regularity::C0>,
            Core::DiffOp::Identity,
            Core::DiffOp::Gradient, 
            Field::ScalarToVector<2, Field::Regularity::C0> 
            >
            (fespace_pml,subspace_pml, Rho, dt*dt, wx.getVector(1), 1.0,phi.getVector(0),
                scalar_to_vect_ez); 

        // Surface wx contribution for left and right boundary
        //phi.getVector(0) += -dt * dt * massMatrixBoundaryWx * wx.getVector(1);
        
        // Contribution of previous timesets for phi
        phi.getVector(0) +=   2.0   * massMatrixPhi * phi.getVector(1);
        phi.getVector(0) += - 1.0   * massMatrix_m  * phi.getVector(2);

        // Solve Phi
        LAL::Solve(massMatrix_p,phi.getVector(0));

        /////// Update wx, phi_hat  ////////////
        // stores (phi^{n+1} + phi^n)/2 in phi_hat
        fespace_pml.MltAddRestriction(0.5,phi.getVector(0),0.0,phi_hat.getVector(0)); 
        fespace_pml.MltAddRestriction(0.5,phi.getVector(1),1.0,phi_hat.getVector(0)); 

        Core::MltAdd<decltype(fespace_pml),decltype(fespace_pml),true,true,false,
            Core::QuadratureOpt::Default, 
            1, //DimU 
            1, //DimV
            1, //DimI
            1, //DimJ
            Field::Identity,
            Core::DiffOp::Gradient,
            Core::DiffOp::Identity,
            Field::Vector<2>,
            Field::Identity>
            (fespace_pml, fespace_pml, Field::IdentityField, Dissip_coeff, phi_hat.getVector(0), 0.0, wx.getVector(0),
                cst_vec_ex, Field::IdentityField); 
        LAL::Solve(massMatrix_wx, wx.getVector(0));

        wx.getVector(0) += (wx.getVector(1))*(1.0/dt - 0.5*Dissip_coeff);
        wx.getVector(0) /= (1.0/dt + 0.5*Dissip_coeff);

        phi_hat.getVector(0) *= Dissip_coeff*dt;
        phi_hat.getVector(0) += phi_hat.getVector(1);
        ////////////////////////////////////////////


        // Get pressure TODO: change to displacement formulation 
        Core::MltAdd< Square<Scaling<2>>,Square<Scaling<2>>,true,true,true, Core::QuadratureOpt::Default,
            1,       //DimU
            1,       //DimV
            2,       //DimI
            1,       //DimJ
            Field::Scalar<Field::Regularity::C0>,
            Core::DiffOp::Gradient,
            Core::DiffOp::Identity, 
            Field::Identity, 
            Field::ScalarToVector<2, Field::Regularity::C0>
            >
            (fespace,fespace, Rho, dt * delta2, phi.getVector(0), 0.0, pressure.getVector(0), Field::IdentityField, scalar_to_vect_ez); 
        pressure.getVector(0) += 1./dt * massMatrixRho * (phi.getVector(0) - 2*phi.getVector(1) + phi.getVector(2));
        LAL::Solve(massMatrixUnit, pressure.getVector(0));
        pressure.getVector(0) += pressure.getVector(1); 
        

        // Writing solution.
        if (iStep % OutputFreq == 0)
        {   
            /////// Write file for one point 
            obsPoint_P << std::to_string(iStep*dt) << " , "; 
            for (Index iObsPoint: iObsPoint_vector){ 
                    obsPoint_P << std::to_string(pressure.getVector(0)(iObsPoint))  << " , " ; 
                }
            obsPoint_P << std::endl ; 
            ////////////////////////////////////

            /*// Write VTK files
            if( iStep % (NVtk * OutputFreq) == 0)
            { 

            // Get velocity
            Core::MltAdd< Square<Scaling<2>>,Square<Scaling<2>>,true,true,true, Core::QuadratureOpt::Default,
            1,       //DimU
            2,       //DimV
            2,       //DimI
            1,       //DimJ
            Field::Identity,
            Core::DiffOp::Gradient,
            Core::DiffOp::Identity
            >
            (fespace,fespace,Field::IdentityField,-1,phi.getVector(0),0.0,velocity); 
            LAL::Solve(massMatrixVelocity,velocity);

            std::cout<< "Writing VTK" << NOutput << " " << iStep*dt << " " << T_end/OutputDeltat <<endl;
            infoFile<< "Writing VTK " << NOutput << " " << iStep*dt << " " << T_end/OutputDeltat <<endl;
            ParallelWriteVTK(FileName+"_Phi."+std::to_string(NOutput) + ".vtk", fespace, phi.getVector(0));
            ParallelWriteVTK(FileName+"_P."+std::to_string(NOutput) + ".vtk", fespace, pressure.getVector(0));
            ParallelWriteVTK(FileName+"_U."+std::to_string(NOutput) + ".vtk", fespace, velocity);
            
            NOutput++;
            }
            */
        }

        // Swapping solutions.
        phi.Swap();
        pressure.Swap();
        phi_hat.Swap();
        wx.Swap();
    }
    
    ParallelWriteVTKTerminate();

    //Output of the time spent
    //time(&end);
    //std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    infoFile << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;

    infoFile.close();

    //return pressure.getVector(0); 

    // convert Eigen::matrix to std::vector
    std::vector<std::vector<double>> vec(pressure.getVector(0).data(), pressure.getVector(0).data() + pressure.getVector(0).size());
    return vec; 
  }

  // Specify that our model supports evaluation. Jacobian support etc. may be indicated similarly.
  bool SupportsEvaluate() override {
    return true;
  }

private:
  int test_delay;
};



int main(){

  // Read environment variables for configuration
  char const* port_cstr = std::getenv("PORT");
  int port = 0;
  if ( port_cstr == NULL ) {
    std::cout << "Environment variable PORT not set! Using port 4242 as default." << std::endl;
    port = 4242;
  } else {
    port = atoi(port_cstr);
  }

  char const* delay_cstr = std::getenv("TEST_DELAY");
  int test_delay = 0;
  if ( delay_cstr != NULL ) {
    test_delay = atoi(delay_cstr);
  }
  std::cout << "Evaluation delay set to " << test_delay << " ms." << std::endl;


  // Set up and serve model
  ExampleModel model(test_delay);
  umbridge::serveModels({&model}, "0.0.0.0", port);

  return 0;
} 