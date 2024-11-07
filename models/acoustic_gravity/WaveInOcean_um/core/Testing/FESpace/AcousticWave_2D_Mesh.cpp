#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  
const Index FEOrder = 10;

 

int main(int argc, char *argv[])
{
    std::cout << "Running AcousticWave_2D_Mesh........................" << std::endl;
    
    std::cout << "Simulation of a transversaly isotropic acoustic wave in a potatoe." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Mesh mesh(argv[1],FormatMsh);

    Mesh quad_mesh = mesh.extract(mesh.getElementList(2));

    FEMesh<2> fespace(quad_mesh,FEOrder);
    
    Index N = fespace.getNumDoFs();
    
    //Pre-compute geometric info  
    fespace.ComputeGeometricInfo();

    
    std::string FileName = argv[2];
  
    // Mesh mesh = fespace.getMesh();
    // std::ofstream stream(FileName + ".mesh");
    // print_MESH(mesh, stream);
    // stream.close();

    std::cout << "The simulation has " << N << " degree of freedom" << std::endl;
    
    Field::Scalar density(1);
    
    std::array<std::array<Real,2>, 2> A;
    A[0] = {2.0, 0.0};
    A[1] = {0.0, 1.0};
    
    Field::Matrix<2,2> velocity(A);

    
    LAL::DiagonalMatrix massMatrix;
    CreateMass(fespace,density,massMatrix,1.0);
    
    Index N_iter;
    Real Residue;
    Real rho;
   
    rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&fespace, &velocity](LAL::Vector &u,LAL::Vector &v) {MltAddStiffness(fespace,velocity,1.0,u,0.0,v);},
                              N,N_iter,Residue);
    
    Real dt =  0.95 * (2.0 / sqrt(rho));
    
    std::cout << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;
    
    std::cout << "\t" << dt << std::endl;
    
    
    //Compute the source
    LAL::Vector source;
    LAL::Allocate(source, N);
    
    for (Index i=0;i<fespace.getNumDoFs();++i)
    {
        RealVector P;
        fespace.getDoFCoordinate(i,P);
    
        RealVector Center = { -1.0, 1.5, 0.};
        source(i) = CInfty0(P,Center,0.2,4);
    }
      
    
    LAL::VectorSequence solution(3);
    solution.Allocate(N);
    
    Real T_end = 10.0;
    Real OutputDeltat = 0.1;
    
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput =0;
 

    //Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;
    for (Index iStep = 0; iStep*dt < T_end; iStep++)
    {
       //Leap-frog scheme
        MltAddStiffness(fespace,velocity,-dt*dt,solution.getVector(1),0.0,solution.getVector(0));
        MltAddMass(fespace,density,2.0,solution.getVector(1),1.0,solution.getVector(0));
        MltAddMass(fespace,density,-1.0, solution.getVector(2), 1.0, solution.getVector(0));
     
        solution.getVector(0) += dt*dt*Ricker(iStep*dt,0.20,500)*massMatrix*source;
     
        LAL::Solve(massMatrix,solution.getVector(0));
  
        // Writing solution.
        if (iStep % OutputFreq == 0 && iStep>0)
        {
   
            // Writing solution
            ParallelWriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", fespace,solution.getVector(0));
            
            //ParallelWriteSolutionVIZIR(FileName, fespace, mesh, solution.getVector(0),  NOutput);
            
            NOutput++;
        }
  
        // Swapping solutions.
        solution.Swap();
    }
     
     ParallelWriteVTKTerminate();
    //  ParallelWriteVIZIRTerminate();
     
     //Output of the time spent
     time(&end);
     std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    
    
}
 
