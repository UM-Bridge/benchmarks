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
    phi(3), phi_hat(2), wx(2), pressure(2),
    Rho(_Rho), Rho_cm2(_Rho_cm2), scalar_to_vect_ez(_field_ez), 
    cst_vec_ez(ez), cst_vec_ex(ex)
  {}

  void Init()
  {
    Index NDoFs = fespace.getNumDoFs();

    ////////// Get boundary spaces
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
    //// Mass matrix for phi 
    CreateMass<1>(fespace, Rho_cm2, massMatrixPhi, 1.0); 
    AssembleMassBndy(surface_fespace, Field::IdentityField, massMatrixPhi, 1./delta2*_Rho({0,1,0})); 

    // Matrix for the pressure
    Field::Scalar unit_field(1.0); 
    // Displ. formulation
    //CreateMass<1>(fespace,unit_field,massMatrixUnit,1.0);
    // Vel. formulation
    CreateMass<1>(fespace,unit_field,massMatrixUnit,1.0);
    CreateMass<1>(fespace,Rho,massMatrixRho,1.0);

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
    phi.Allocate(NDoFs);
    pressure.Allocate(NDoFs);
    //LAL::Allocate(pressure, NDoFs);
    LAL::Allocate(sourceVol, NDoFs); 
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
  

  // Define input and output dimensions of model (here we have a single vector of length 1 for input; same for output)
  std::vector<std::size_t> GetInputSizes(const json& config) const override {
    if(config.is_null() || config.empty()) {
      return {Nx}; 
    }
    std::vector<int> fixedFloor = config["fixedFloor"].get<std::vector<int>>();
    return {Nx-fixedFloor.size()};
  
  }

  std::vector<std::size_t> GetOutputSizes(const json& config) const override {
    return {105};
  }

  std::vector<std::vector<double>> Evaluate(const std::vector<std::vector<double>>& inputs, json config) override {
    // Do the actual model evaluation
    if(config.is_null() || config.empty()) 
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
        //std::cout<< dataFunctionSpace[i] << std::endl;
      }

      
    }
    
    // reset values for phi and pressure vectors
    std::fill(phi.getVector(0).begin(), phi.getVector(0).end(), 0);
    std::fill(phi.getVector(1).begin(), phi.getVector(1).end(), 0);
    std::fill(phi.getVector(2).begin(), phi.getVector(2).end(), 0);
    std::fill(pressure.getVector(0).begin(), pressure.getVector(0).end(), 0);
    std::fill(pressure.getVector(1).begin(), pressure.getVector(1).end(), 0);

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
      
        LAL::Solve(massMatrixPhi, phi.getVector(0));

        // Contribution of previous timesets for phi
        phi.getVector(0) +=   2.0 * phi.getVector(1);
        phi.getVector(0) += - 1.0 * phi.getVector(2);

        // Get pressure, velocity formulation
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
          output[0].push_back( pressure.getVector(0)(iObsPoint_vector[0]) ) ;
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
  LAL::DiagonalMatrix massMatrixPhi, massMatrixRho;
  LAL::DiagonalMatrix massMatrixPhi_Dissip, massMatrixPhi_Bndy_Dissip, massMatrix_wx;
  LAL::DiagonalMatrix massMatrix_m, massMatrix_p;
  LAL::DiagonalMatrix massMatrixUnit; 
  LAL::VectorSequence phi;
  LAL::VectorSequence pressure;
  LAL::VectorSequence phi_hat, wx; 
  LAL::Vector sourceVol;
  //LAL::Vector pressure; 
  std::vector<Index> iObsPoint_vector;
  
  Field::Scalar<Field::Regularity::C0> Rho_cm2;
  Field::Scalar<Field::Regularity::C0> Rho;
  Field::Vector<2,Field::Regularity::Constant> cst_vec_ez;
  Field::Vector<2,Field::Regularity::Constant> cst_vec_ex;
  Field::ScalarToVector<2, Field::Regularity::C0> scalar_to_vect_ez; 
  
  float dt;
  int OutputFreq;

  Real functionSpace(const RealVector &p){
    // Input: p = (x,z),  x \in [0,Lx], z \in [0,Lz]
    Real result = function_interpolate(dataSpace, dataFunctionSpace, p[0]);
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