#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <algorithm>
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

  ExampleModel() 
    : umbridge::Model("forward"),
    scaling(Lx, Lz), fespace(FEOrderX,FEOrderZ,Nx,Nz,scaling), 
    phi(3), phi_hat(2), wx(2), // pressure(2), // old pressure computation 
    Rho(_Rho), Rho_cm2(_Rho_cm2), scalar_to_vect_ez(_field_ez), 
    field_ez(ez), pml_dissipation_dir_field(ex), Dissip_Disc_Field({{0,0.0}, {1,Dissip_coeff}}) // for PML 
  {}

  void Init() {
    Index NDoFs = fespace.getNumDoFs();

    ///////////  Construct the matrices
    // The mass matrices are defined with weighted integration, the weight can be a scalar field or a constant 
    // e.g. CreateMass<1>(fespace, Rho_c2, massMatrixPhi_Sparse, 1.0); second parameter is the scalar field, last parameter is the constant. 
    //// Mass matrix for phi 
    auto surface_fespace = fespace.getBndyFESpace(2);
    CreateMass<1>(fespace, Rho_cm2, massMatrixPhi, 1.0); 
    AssembleMassBndy(surface_fespace, Field::IdentityField, massMatrixPhi, 1./delta2*_Rho({0,1,0})); 

    // Matrix for the pressure old version 
    // Field::Scalar unit_field(1.0); 
    // CreateMass<1>(fespace,unit_field,massMatrixUnit,1.0);
    // CreateMass<1>(fespace,Rho,massMatrixRho,1.0);

    // Matrix for the pressure new version 
    CreateMass<1>(fespace, Rho, massMatrixRho, 1.0);
    CreateMass<1>(fespace, Field::IdentityField, massMatrixUnit, 1.0);
    // Compute the product M(1)^{-1} M(rho) 
    invMassMatrixUnit = massMatrixUnit;
    LAL::Invert(invMassMatrixUnit);
    matrixPressure = invMassMatrixUnit * massMatrixRho;


    // Matrices for the PML
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
    CreateMassDG(fespace_pml, Field::IdentityField, massMatrix_wx, 1.0);
    //////////////////////////////////////

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
    phi.Allocate(NDoFs);
    //pressure.Allocate(NDoFs);
    LAL::Allocate(sourceVol, NDoFs); 

    // new pressure computation 
    LAL::Allocate(pressure, NDoFs);
    LAL::Allocate(pressure_int, NDoFs);
    LAL::Allocate(tempPressure, NDoFs);
    ////////////////////////////////////////

    ///////////////// For PML /////////////
    phi_hat.Allocate(NDoFs);
    wx.Allocate(NDoFs);

    ////////////////////////////////////////
    // Time and output 
    OutputFreq = floor(OutputDeltat/dt); 
    ////////////////////////////////////////

    // create dataSpace = the vector of points where the function dataFunction is known
    // For now dataSpace = [i for i in range(Nx)]
    for (int i=0; i<Nx+1; i++)
    {
      dataSpace.push_back(i); 
    }
  }



  std::vector<std::size_t> GetInputSizes(const json& config) const override {
    //if(config.is_null() || config.empty()) {
    if(config["fixedFloor"].is_null() || config["fixedFloor"].empty()) {
      return {Nx}; 
    }
    std::vector<int> fixedFloor = config["fixedFloor"].get<std::vector<int>>();
    return {Nx-fixedFloor.size()};
  }

  std::vector<std::size_t> GetOutputSizes(const json& config) const override {
    std::vector<std::vector<double>> captorCoordinates = config["captors"].get<std::vector<std::vector<double>>>();
    return {159*captorCoordinates.size()};
  }

  std::vector<std::vector<double>> Evaluate(const std::vector<std::vector<double>>& inputs, json config) override {
    // Do the actual model evaluation

    // Get the seabed source from the config file 
    if(config["fixedFloor"].is_null() || config["fixedFloor"].empty())
    {
      dataFunctionSpace = inputs[0];
    }
    else 
    {
      dataFunctionSpace = {};
      int i_input=0;
      std::vector<int> fixedFloor = config["fixedFloor"].get<std::vector<int>>();
      for (Index i=0;i<inputs[0].size() + fixedFloor.size() ;++i)
      {
        if( std::find(fixedFloor.begin(), fixedFloor.end(), i )!=fixedFloor.end())
        { dataFunctionSpace.push_back(0.); }
        else 
        {
          dataFunctionSpace.push_back(inputs[0][i_input]);
          i_input ++; 
        }
      }
    }

    // Get the captors location from the config file
    assert(! config["captors"].is_null() && ! config["captors"].empty());
    std::vector<Index> iObsPoint_vector;
    std::vector<std::vector<double>> captorCoordinates = config["captors"].get<std::vector<std::vector<double>>>();
    for (const std::vector<double>& captor : captorCoordinates)
    { 
      RealVector q{captor[0]/Lx, captor[1]/Lz};
      Index iObsPoint = findGlobLocIndex(q, fespace);
      //std::cout << "(" << captor[0] << "," << captor[1] << ") ; (" << q[0] << "," << q[1] << ") " << iObsPoint << " " << fespace.getNumDoFs() << std::endl; 
      // If iObsPoint the given point is really in the domain and we store it
      if(iObsPoint < fespace.getNumDoFs())
      {   
          iObsPoint_vector.push_back(iObsPoint);
      }
    }
    assert(iObsPoint_vector.size()>0); 
    //std::cout<< iObsPoint_vector.size() << std::endl; 

    
    // Would be nice to have the next lines in Init(), but not possible due to absence of constructor (?)
    auto subspace_pml = fespace.getSubDomainFESpace([](RealVector xyz){
            if ((xyz[0] > Lx-LxPml_R) || (xyz[0] < LxPml_L))
                return true;
            else return false;
    });
    auto fespace_pml = subspace_pml;
    subspace_pml.setAsSubSpace(); 
    fespace_pml.setAsFESpace(); 
    //////////////////////////////////////


    // reset values for phi and pressure vectors
    std::fill(phi.getVector(0).begin(), phi.getVector(0).end(), 0);
    std::fill(phi.getVector(1).begin(), phi.getVector(1).end(), 0);
    std::fill(phi.getVector(2).begin(), phi.getVector(2).end(), 0);

    // Old pressure computation: reset values
    // std::fill(pressure.getVector(0).begin(), pressure.getVector(0).end(), 0);
    // std::fill(pressure.getVector(1).begin(), pressure.getVector(1).end(), 0);

    // New pressure computation: reset values
    std::fill(pressure_int.begin(), pressure_int.end(), 0);
    std::fill(pressure.begin(), pressure.end(), 0);

    std::vector<std::vector<double>> output(1);
    output[0] = {}; 
    //////////////////////////////////// 



    
    ////////// Create source here because it depends on the input
    auto bottom_fespace = fespace.getBndyFESpace(0);
    LAL::Vector sourceCoord, source;
    LAL::Allocate(sourceCoord, bottom_fespace.getNumDoFs());
    LAL::Allocate(source, bottom_fespace.getNumDoFs());
    Field::Scalar unit_field(1.0); 

    //Compute f(x) at each DoF of the seabed 
    LAL::SparseMatrix bottomMassMatrix;
    CreateMassBndy(bottom_fespace, unit_field, bottomMassMatrix, _Rho({0,0,0}));
    bottom_fespace.setAsFESpace();
    for (Index i=0;i<bottom_fespace.getNumDoFs();++i)
    {
        Index dof = bottom_fespace._glob_bdny2glob_sq(i);
        RealVector P;
        // Get the (x,z) coordinates to evaluate the function sourceSpace
        fespace.getDoFCoordinate(dof,P);
        sourceCoord(i) = functionSpace(P);
    }
    // Get FE approximation of f(x)
    source = bottomMassMatrix * sourceCoord;
    // Store the result in a vector of the size of the whole domain (NDoFs entries) to simplify the sum 
    for (Index i=0;i<bottom_fespace.getNumDoFs();++i)
    {
        Index dof = bottom_fespace._glob_bdny2glob_sq(i);
        sourceVol(dof) = source(i) ;
    }
    //////////////////////////////////

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
        phi.getVector(0) += dt*dt * sourceVol * functionTime(dt*iStep);
      

        /*// Without PML
        // LAL::Solve(massMatrixPhi, phi.getVector(0));

        // // Contribution of previous timesets for phi
        // phi.getVector(0) +=   2.0 * phi.getVector(1);
        // phi.getVector(0) += - 1.0 * phi.getVector(2);
        */

        ///////////////////////
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
                field_ez, field_ez); 
        
        // Contribution of wx
        Core::MltAdd<decltype(fespace_pml),decltype(subspace_pml),false,false,true, 
            Core::QuadratureOpt::Default,
            1,       //DimU
            1,       //DimV
            1,       //DimI
            1,       //DimJ
            Field::Scalar<Field::Regularity::C0>,
            Core::DiffOp::Identity,
            Core::DiffOp::Gradient, 
            Field::Identity,
            Field::Vector<2> 
            >
            (fespace_pml,subspace_pml, Rho, dt*dt, wx.getVector(1), 1.0, phi.getVector(0),
                Field::IdentityField, pml_dissipation_dir_field); 
        
        // Contribution of previous timesteps for phi
        phi.getVector(0) +=   2.0 * massMatrixPhi * phi.getVector(1);
        phi.getVector(0) += - 1.0 * massMatrix_m  * phi.getVector(2);

        // Solve Phi
        LAL::Solve(massMatrix_p,phi.getVector(0));


        /////// PML - Update wx, phi_hat  ////////////
        // stores (phi^{n+1} + phi^n)/2 in phi_hat
        fespace_pml.MltAddRestriction(0.5,phi.getVector(0), 0.0,phi_hat.getVector(0)); 
        fespace_pml.MltAddRestriction(0.5,phi.getVector(1), 1.0,phi_hat.getVector(0)); 

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
            (fespace_pml, fespace_pml, Field::IdentityField, -Dissip_coeff, phi_hat.getVector(0), 0.0, wx.getVector(0),
                pml_dissipation_dir_field, Field::IdentityField); 
        LAL::Solve(massMatrix_wx, wx.getVector(0));

        wx.getVector(0) += (wx.getVector(1))*(1.0/dt - 0.5*Dissip_coeff);
        wx.getVector(0) /= (1.0/dt + 0.5*Dissip_coeff);

        phi_hat.getVector(0) *= Dissip_coeff*dt;
        phi_hat.getVector(0) += phi_hat.getVector(1);
        //////////////// END PML


        // Get pressure, velocity formulation, old computation 
        // Core::MltAdd< Square<Scaling<2>>,Square<Scaling<2>>,true,true,true, Core::QuadratureOpt::Default,
        //     1,       //DimU
        //     1,       //DimV
        //     2,       //DimI
        //     1,       //DimJ
        //     Field::Scalar<Field::Regularity::C0>,
        //     Core::DiffOp::Gradient,
        //     Core::DiffOp::Identity, 
        //     Field::Identity, 
        //     Field::ScalarToVector<2, Field::Regularity::C0>
        //     >
        //     (fespace,fespace, Rho, dt * delta2, phi.getVector(0), 0.0, pressure.getVector(0), Field::IdentityField, scalar_to_vect_ez); 
        // pressure.getVector(0) += 1./dt * massMatrixRho * (phi.getVector(0) - 2*phi.getVector(1) + phi.getVector(2));
        // LAL::Solve(massMatrixUnit, pressure.getVector(0));
        // pressure.getVector(0) += pressure.getVector(1); 
        
        // Pressure, new computation 
        // Compute auxiliary variable for the pressure 
        Core::MltAdd< Square<Scaling<2>>,Square<Scaling<2>>,true,true,true, Core::QuadratureOpt::Default,
            1,       //DimU
            1,       //DimV
            1,       //DimI
            1,       //DimJ
            Field::Scalar<Field::Regularity::C0>,
            Core::DiffOp::Gradient,
            Core::DiffOp::Identity, 
            Field::Vector<2>,
            Field::Identity
            >
            (fespace, fespace, Rho, dt , phi.getVector(1), 1, pressure_int,
                field_ez, Field::IdentityField); 

        // Writing solution.
        if (iStep % OutputFreq == 0)
        { 
          // Pressure, new computation 
          // For a velocity input: compute d_t phi and d_t psi
          tempPressure  =   phi.getVector(0);
          tempPressure += - phi.getVector(2);
          tempPressure  =   tempPressure/(2*dt);

          pressure = matrixPressure * (tempPressure + delta2 * pressure_int) ;
          for (Index iObs=0;iObs<iObsPoint_vector.size() ;++iObs)
          { 
            // old pressure computation 
            //output[0].push_back( pressure.getVector(0)(iObsPoint_vector[iObs]) ) ;

            // new pressure computation 
            output[0].push_back( pressure(iObsPoint_vector[iObs]) ) ;
          }
        }
        // Swapping solutions.
        phi.Swap();
    } 
    return output; 
  }

  // Specify that our model supports evaluation. Jacobian support etc. may be indicated similarly.
  bool SupportsEvaluate() override {
    return true;
  }

private:
  Scaling<2> scaling;
  Square<Scaling<2>> fespace;
  //FESpaceTSubdomain< Square<Scaling<2>> ,2>* subspace_pml; //(Square<Scaling<2>> ,2> fespace, std::vector<Index> {0}); 
  //FESpaceTSubdomain< Square<Scaling<2>> ,2> 
  //auto fespace_pml(); 

  std::vector<Real> dataSpace; 
  std::vector<Real> dataFunctionSpace;

  // In private attribute are all the element needed in the Evaluate() function
  LAL::DiagonalMatrix massMatrixPhi; //, massMatrixRho;
  LAL::DiagonalMatrix massMatrixPhi_Dissip, massMatrixPhi_Bndy_Dissip, massMatrix_wx;
  LAL::DiagonalMatrix massMatrix_m, massMatrix_p;
  LAL::DiagonalMatrix massMatrixUnit; 
  LAL::VectorSequence phi;

  // New pressure computation 
  LAL::SparseMatrix massMatrixRho, matrixPressure;
  LAL::DiagonalMatrix invMassMatrixUnit; 
  LAL::Vector pressure; 
  LAL::Vector pressure_int;
  LAL::Vector tempPressure; 

  // Old pressure computation 
  //LAL::VectorSequence pressure;


  LAL::VectorSequence phi_hat, wx; 
  LAL::Vector sourceVol;
  
  Field::Scalar<Field::Regularity::C0> Rho_cm2;
  Field::Scalar<Field::Regularity::C0> Rho;
  Field::Scalar<Field::Regularity::PiecewiseConstant> Dissip_Disc_Field; // for PML
  Field::Vector<2,Field::Regularity::Constant> field_ez;    // For PML
  Field::Vector<2,Field::Regularity::Constant> pml_dissipation_dir_field; // For PML
  Field::ScalarToVector<2, Field::Regularity::C0> scalar_to_vect_ez; 

  std::map<Index, Real> Dissip_map; // for PML

  float dt;
  int OutputFreq;

  Real functionSpace(const RealVector &p){
    // Input: p = (x,z), x \in [0,Lx],  z \in [0,Lz]
    // /!\ This function_interpolate is used differently the rest of OndomathX. 
    // Usually the input is x \in [0,Lx], z \in [0,Lz].
    //    Here the input is x \in [0,Nx], z \in [0,Nz].
    Real result = function_interpolate(dataSpace, dataFunctionSpace, p[0]*Nx/Lx);
    //std::cout<<p[0]<<" "<< p[0]*Nx/Lx<<" "<<result<<std::endl;
    return result;
   }  

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

  // Set up and serve model
  ExampleModel model;
  model.Init();
  umbridge::serveModels({&model}, "0.0.0.0", port);

  return 0;
} 
