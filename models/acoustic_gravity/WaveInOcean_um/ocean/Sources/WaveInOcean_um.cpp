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

  ExampleModel(int test_delay) : umbridge::Model("forward"), 
  test_delay(test_delay),
  scaling(Lx, Lz), fespace(FEOrderX,FEOrderZ,Nx,Nz,scaling), 
  phi(3), phi_hat(2), wx(2), 
  Rho(_Rho), Rho_cm2(_Rho_cm2), scalar_to_vect_ez(_field_ez), 
  cst_vec_ez(ez), cst_vec_ex(ex)
  {
  }

  void Init()
  {
    Center[1] = Lz*Center[1];
    Index NDoFs = fespace.getNumDoFs();

    ////////// Get boundary spaces
    auto bottom_fespace = fespace.getBndyFESpace(0);
    auto surface_fespace = fespace.getBndyFESpace(2);

    for (Index i = 0; i < nPointX; i++){ 
        for (Index j = 0; j < nPointZ; j++){
            RealVector q{listX[i]/Lx, listZ[j]/Lz};
            Index iObsPoint = findGlobLocIndex(q, fespace);
            // If iObsPoint the given point is really in the domain and we store it
            if(iObsPoint < fespace.getNumDoFs())
            {   
                iObsPoint_vector.push_back(iObsPoint);
            }
        }
    }

    ///////////  Construct the matrices
    // The mass matrices are defined with weighted integration, the weight can be a scalar field or a constant 
    // e.g. CreateMass<1>(fespace, Rho_c2, massMatrixPhi_Sparse, 1.0); second parameter is the scalar field, last parameter is the constant. 
    Field::Scalar unit_field(1.0); 

    //// Mass matrix for phi 
    CreateMass<1>(fespace, Rho_cm2, massMatrixPhi, 1.0); 
    AssembleMassBndy(surface_fespace, Field::IdentityField, massMatrixPhi, 1./delta2*_Rho({0,1,0})); 

    // Matrices for the pressure
    CreateMass<1>(fespace,unit_field,massMatrixUnit,1.0);
    ////////////////////////////////////////

    ////////////////////////////////////////
    // Computation of the optimal time step
    // Remark: need to pass explicitely by reference massMatrixPhi since it is a class (private) member
    Index N_iter;
    Real Residue;
    Real rho;
    rho = LAL::PowerIteration([&massMatrixPhi=massMatrixPhi] (LAL::Vector &v) { LAL::Solve(massMatrixPhi,v);},
                              [&fespace=fespace, &Rho=Rho]      (LAL::Vector &u,LAL::Vector &v)
                              {  MltAddStiffness(fespace,Rho,1,u,0,v);},
                              NDoFs,N_iter,Residue);
    dt =  0.95 * (2.0 / sqrt(rho));
    ////////////////////////////////////////
      
    ////////////////////////////////////////
    // Vectors
    LAL::Vector sourceCoord, source;
    phi.Allocate(NDoFs);
    LAL::Allocate(pressure, NDoFs);
    LAL::Allocate(sourceVol, NDoFs); 
    LAL::Allocate(sourceCoord, bottom_fespace.getNumDoFs());
    LAL::Allocate(source, bottom_fespace.getNumDoFs());

    // Compute f(x) at each DoF of the seabed 
    LAL::SparseMatrix bottomMassMatrix;
    CreateMassBndy(bottom_fespace, unit_field, bottomMassMatrix, _Rho({0,0,0}));
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
    phi_hat.Allocate(NDoFs);
    wx.Allocate(NDoFs);

    /*
    std::map<Index, Real> Dissip_map;
    Dissip_map[0] = 0.0;
    Dissip_map[1] = Dissip_coeff;
    Field::Scalar<Field::Regularity::PiecewiseConstant> Dissip_Disc_Field(Dissip_map);
    LAL::DiagonalMatrix massMatrixPhi_Dissip, massMatrixPhi_Bndy_Dissip;
    auto subspace_pml = fespace.getSubDomainFESpace([](RealVector xyz){
            if ((xyz[0] > Lx-LxPml_R) || (xyz[0] < LxPml_L))
                return true;
            else return false;
    });
    auto fespace_pml = subspace_pml;  
    subspace_pml.setAsSubSpace(); 
    fespace_pml.setAsFESpace(); 
    CreateMass<1>(subspace_pml, Rho_cm2, massMatrixPhi_Dissip, Dissip_coeff);
    CreateMassBndy(surface_fespace, Dissip_Disc_Field, massMatrixPhi_Bndy_Dissip, 1./delta2*rhoSurface);
    massMatrix_p = massMatrixPhi + 0.5*dt*massMatrixPhi_Dissip + 0.5*dt*massMatrixPhi_Bndy_Dissip;
    massMatrix_m = massMatrixPhi - 0.5*dt*massMatrixPhi_Dissip - 0.5*dt*massMatrixPhi_Bndy_Dissip;
    CreateMassDG(fespace_pml, Field::IdentityField, massMatrix_wx,1.0);
    */
    //////////////////////////////////////

    ////////////////////////////////////////
    // Time and output 
    OutputFreq = floor(OutputDeltat/dt); 
  }
  

  // Define input and output dimensions of model (here we have a single vector of length 1 for input; same for output)
  std::vector<std::size_t> GetInputSizes(const json& config_json) const override {
    return {1};
  }

  std::vector<std::size_t> GetOutputSizes(const json& config_json) const override {
    return {53};
  }

  std::vector<std::vector<double>> Evaluate(const std::vector<std::vector<double>>& inputs, json config) override {
    // Do the actual model evaluation; here we just multiply the first entry of the first input vector by two, and store the result in the output.
    // In addition, we support an artificial delay here, simulating actual work being done.
    // std::this_thread::sleep_for(std::chrono::milliseconds(test_delay));
    //return {{inputs[0][0] * 2.0}};

    std::cout<< "test access fespace: " << fespace.getNumDoFs() << std::endl; 

    Center[0] = Lx*inputs[0][0];
    std::vector<std::vector<double>> output(1);
    output[0] = {}; 

    // TODO: I don't know how to create subspace_pml in the constructor or in the Init()
    // function, since there is no default constructor and no assignment function 
    // for the class FESpaceTSubdomain
    std::map<Index, Real> Dissip_map;
    Dissip_map[0] = 0.0;
    Dissip_map[1] = Dissip_coeff;
    Field::Scalar<Field::Regularity::PiecewiseConstant> Dissip_Disc_Field(Dissip_map);
    
    auto surface_fespace = fespace.getBndyFESpace(2);
    auto subspace_pml = fespace.getSubDomainFESpace([](RealVector xyz){
            if ((xyz[0] > Lx-LxPml_R) || (xyz[0] < LxPml_L))
                return true;
            else return false;
    });
    auto fespace_pml = subspace_pml;  
    subspace_pml.setAsSubSpace(); 
    fespace_pml.setAsFESpace(); 

    CreateMass<1>(subspace_pml, Rho_cm2, massMatrixPhi_Dissip, Dissip_coeff);
    CreateMassBndy(surface_fespace, Dissip_Disc_Field, massMatrixPhi_Bndy_Dissip, 1./delta2*_Rho({0,1,0}));
    massMatrix_p = massMatrixPhi + 0.5*dt*massMatrixPhi_Dissip + 0.5*dt*massMatrixPhi_Bndy_Dissip;
    massMatrix_m = massMatrixPhi - 0.5*dt*massMatrixPhi_Dissip - 0.5*dt*massMatrixPhi_Bndy_Dissip;
    CreateMassDG(fespace_pml, Field::IdentityField, massMatrix_wx,1.0);
    ////////////////////////////////////  
    std::cout << "Running WaveInOcean........................" << std::endl;

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
        phi.getVector(0) += dt * dt * sourceVol * functionTime(dt*iStep);

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
        
        // Writing solution.
        if (iStep % OutputFreq == 0)
        {   
          // Get pressure, displacement formulation 
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
              (fespace,fespace, Rho, delta2, phi.getVector(0), 0.0, pressure, Field::IdentityField, scalar_to_vect_ez); 
          LAL::Solve(massMatrixUnit, pressure);
          pressure += 1./(dt*dt) * (phi.getVector(0) - 2*phi.getVector(1) + phi.getVector(2));
          output[0].push_back( pressure(iObsPoint_vector[0]) ) ;
        }

        // Swapping solutions.
        phi.Swap();
        phi_hat.Swap();
        wx.Swap();
    } 

    return output; 
  }

  // Specify that our model supports evaluation. Jacobian support etc. may be indicated similarly.
  bool SupportsEvaluate() override {
    return true;
  }

private:
  int test_delay;
  
  Scaling<2> scaling;
  Square<Scaling<2>> fespace;
  //FESpaceTSubdomain< Square<Scaling<2>> ,2>* subspace_pml; //(Square<Scaling<2>> ,2> fespace, std::vector<Index> {0}); 
  //FESpaceTSubdomain< Square<Scaling<2>> ,2> 
  //auto fespace_pml(); 

  // In private attribute are all the element needed in the Evaluate() function
  LAL::DiagonalMatrix massMatrixPhi;
  LAL::DiagonalMatrix massMatrixPhi_Dissip, massMatrixPhi_Bndy_Dissip, massMatrix_wx;
  LAL::DiagonalMatrix massMatrix_m, massMatrix_p;
  LAL::DiagonalMatrix massMatrixUnit; 
  LAL::VectorSequence phi;
  LAL::VectorSequence phi_hat, wx; 
  LAL::Vector sourceVol;
  LAL::Vector pressure; 
  std::vector<Index> iObsPoint_vector;
  
  Field::Scalar<Field::Regularity::C0> Rho_cm2;
  Field::Scalar<Field::Regularity::C0> Rho;
  Field::Vector<2,Field::Regularity::Constant> cst_vec_ez;
  Field::Vector<2,Field::Regularity::Constant> cst_vec_ex;
  Field::ScalarToVector<2, Field::Regularity::C0> scalar_to_vect_ez; 
  
  Real dt;
  Index OutputFreq;

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
  model.Init();
  umbridge::serveModels({&model}, "0.0.0.0", port);

  return 0;
} 