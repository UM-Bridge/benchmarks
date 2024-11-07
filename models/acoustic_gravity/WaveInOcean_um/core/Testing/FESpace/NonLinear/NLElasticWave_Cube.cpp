#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  
const Index FEOrder = 6;

int main(int argc, char *argv[])
{
    std::cout << "Running NLElasticWave_Cube........................" << std::endl;
    
    std::cout << "Simulation of a non linear wave in a cube." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Cube fespace(6,6,6,20,20,20);
    
    Index N = fespace.getNumDoFs();
      
    std::string FileName = argv[1];
  
    std::cout << "The simulation has " << N << " degree of freedom" << std::endl;
    
    Field::Scalar density(1);
    
    Field::IsotropicTensor isotropic_material(2.0,1.0);

    Field::NeoHookeanTensor nl_material(2.0,1.0);
    
 

    LAL::DiagonalMatrix massMatrix;
    CreateMass<3>(fespace,density,massMatrix,1.0);
    
   
    Index N_iter;
    Real Residue;
    Real rho;
   

    // Time step for an isotropic case with the same Lamé Parameter
    rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&fespace, &isotropic_material](LAL::Vector &u,LAL::Vector &v) {MltAddStiffness<3>(fespace,isotropic_material,1.0,u,0.0,v);},
                              3*N,N_iter,Residue);
    
    Real dt =  0.5*0.95 * (2.0 / sqrt(rho));
    
    std::cout << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;
    
    std::cout << "\t" << dt << std::endl;
    ///////////////////////////////////

    
    //Compute the source////////////////////
    LAL::Vector source;
    LAL::Allocate(source, 3*N);
    
    for (Index i=0;i<N;++i)
    {
        RealVector P;
        fespace.getDoFCoordinate(i,P);
    
        RealVector Center = { 0.5, 0.5, 0.5};
        source(3*i+2) = 1000*CInfty0(P,Center,0.05,1); // Push in the Z direction
    }
     ///////////////////////////////////


    LAL::VectorSequence displacement(3);
    displacement.Allocate(3*N);

    LAL::Vector velocity;
    LAL::Allocate(velocity, 3*N);
    
    Real T_end = 1.0;
    Real OutputDeltat = 0.1;
    
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput =0;
 

    //Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;
    for (Index iStep = 0; iStep*dt < T_end; iStep++)
    {
        //Non linear leap-frog scheme
        Core::NonLinearStiffnessMltAdd<3,2>(fespace,nl_material,-dt,displacement.getVector(1),velocity,0.0,displacement.getVector(0),dt);        

        //MltAddStiffness<3>(fespace,isotropic_material,-dt*dt,displacement.getVector(1),0.0,displacement.getVector(0));

        LAL::Solve(massMatrix,displacement.getVector(0));

        displacement.getVector(0) += 2*displacement.getVector(1)-displacement.getVector(2)+dt*dt*Ricker(iStep*dt,0.20,2000)*source;

        velocity = 2*(displacement.getVector(0)-displacement.getVector(1))/dt  - velocity;

        // Writing solution.
        if (iStep % OutputFreq == 0 && iStep>0)
        {
            // Writing solution
            ParallelWriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", fespace,displacement.getVector(0));
            
            NOutput++;
        }
  
        // Swapping solutions.
        displacement.Swap();
    }
     
    ParallelWriteVTKTerminate();

    //Output of the time spent
    time(&end);
    std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
}
 
