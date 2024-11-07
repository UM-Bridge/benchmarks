#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;


 

int main(int argc, char *argv[])
{
    std::cout << "Running AcousticWave_Square........................" << std::endl;
    
    std::cout << "Simulation of a transversaly isotropic acoustic wave in a square." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Square fespace(15,15,100,100);
    Index N = fespace.getNumDoFs();
    
    std::cout << "The simulation has " << N << " degree of freedom" << std::endl;
     
    RealMatrix2x2 velocity_mat;

    velocity_mat[0][0] = 1.0;  velocity_mat[0][1] = 0.0;
    velocity_mat[1][0] = 0.0;  velocity_mat[1][1] = 1.0;
 
    Field::Matrix<2,2> velocity(velocity_mat);
    Field::Scalar density(1.0);
    

    LAL::DiagonalMatrix massMatrix;
    CreateMass(fespace,density,massMatrix,1.0);
    

    Index N_iter;
    Real Residue;
    Real rho;
   
    rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&fespace,&velocity](LAL::Vector &u,LAL::Vector &v) {MltAddStiffness(fespace,velocity,1.0,u,0.0,v);},
                              N,N_iter,Residue);
    
    Real dt =  0.95 * (2.0 / sqrt(rho));
    
    std::cout << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;
    
    std::cout << "\t" << dt << std::endl;
    
    
    //Compute the source
    LAL::Vector initial_velocity;
    LAL::Allocate(initial_velocity, N);

    Real kx = 200;
    
    for (Index i=0;i<fespace.getNumDoFs();++i)
    {
        RealVector P;
        fespace.getDoFCoordinate(i,P);
    
        
        initial_velocity(i) = Gaussian(P[0],0.9,5000)*Gaussian(P[1],0.5,5000)*cos(2*M_PI*kx*(P[0]+P[1]));
    }


    LAL::VectorSequence solution(3);
    solution.Allocate(N);

    //Init first iterate from initial velocity
    MltAddStiffness(fespace,velocity,-(dt*dt*dt)/6.0,initial_velocity,0.0,solution.getVector(1));
    LAL::Solve(massMatrix,solution.getVector(1));
    solution.getVector(1) += dt*initial_velocity;
    
    Real T_end = 1.0;
    Real OutputDeltat = 0.1;
    
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput =0;
 
    std::string FileName = argv[1];
  
    Mesh mesh = fespace.getMesh();
    std::ofstream stream(FileName + ".mesh");
    print_MESH(mesh, stream);
    stream.close();
    
    //Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;
    for (Index iStep = 0; iStep*dt < T_end; iStep++)
    {
       //Leap-frog scheme
        MltAddStiffness(fespace,velocity,-dt*dt,solution.getVector(1),0.0,solution.getVector(0));

        LAL::Solve(massMatrix,solution.getVector(0));

        solution.getVector(0) += 2*solution.getVector(1)-solution.getVector(2); 
  
       
  
        // Writing solution.
        if (iStep % OutputFreq == 0)
        {
   
            // Writing solution
            ParallelWriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", fespace,solution.getVector(0));
            ParallelWriteSolutionVIZIR(FileName, fespace, mesh, solution.getVector(0),  NOutput);
            
            NOutput++;
        }
  
        // Swapping solutions.
        solution.Swap();
    }
     
    ParallelWriteVTKTerminate();
    ParallelWriteVIZIRTerminate();
     
     //Output of the time spent
     time(&end);
     std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    
    
}
 
