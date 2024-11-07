#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
 
Index Order = 4;
Index Order_pressure = 2;

Real lambda = 500;
Index NElements_X = 50;
Index NElements_Y = 50;

int main(int argc, char *argv[])
{
    std::cout << "Running IncompressibleElasticWave_Square........................" << std::endl;
    
    std::cout << "Simulation of an incompressible isotropic elastic wave in a square using or not mixed finite elements." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Scaling<2> scaling(1.0,1.0);
 
    Square<Scaling<2>> fespace(Order,Order,NElements_X,NElements_Y,scaling);
    Square<Scaling<2>> fespace_pressure(Order_pressure,Order_pressure,NElements_X,NElements_Y,scaling);
    
    Index N = 2*fespace.getNumDoFs();
    
    std::cout << "The simulation has " << N << " degree of freedom" << std::endl;

    Field::Scalar density(1.0);

    Field::IsotropicTensor isotropic_material(lambda,1.0);
    Field::IsotropicTensor isotropic_material_no_lambda(0,1.0);

    
    LAL::DiagonalMatrix massMatrix;
    LAL::DiagonalMatrix massMatrixPressure;
   
    
    CreateMass<2>(fespace,density,massMatrix,1.0);
    CreateMassDG(fespace_pressure,density,massMatrixPressure,1.0/lambda);
    
    //GaussLobattoInterpolator aGLI(fespace_pressure.getFE(),fespace.getFE());
    
    Index N_iter;
    Real Residue;
    Real rho;
       
    rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&fespace, &isotropic_material](LAL::Vector &u,LAL::Vector &v) {MltAddStiffness<2>(fespace,isotropic_material,1.0,u,0.0,v);},
                              N,N_iter,Residue);
     
    std::cout << "Spectral radius: " << rho << std::endl;
   
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
    
        RealVector Center = { 0.5, 0.5};
        
        source(2*i) = 5000*CInfty0(P,Center,0.05,1);
    }
      
    LAL::VectorSequence solution(3);
    LAL::VectorSequence solution_mixed(3);
    LAL::Vector pressure_tmp;
    LAL::Vector pressure;

    LAL::Allocate(pressure_tmp,fespace.getNumDoFsDisc());
    LAL::Allocate(pressure,fespace_pressure.getNumDoFsDisc());

    solution.Allocate(N);
    solution_mixed.Allocate(N);

    
    LAL::Vector diff;
    
    Real T_end = 1.0;
    Real OutputDeltat = 0.1;
    
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput=0;
  
    std::string FileName = argv[1];
  
  
    //Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;
    for (Index iStep = 0; iStep*dt < T_end; iStep++)
    {
       //Leap-frog scheme
        MltAddStiffness<2>(fespace,isotropic_material,-dt*dt,solution.getVector(1),0.0,solution.getVector(0));

        LAL::Solve(massMatrix,solution.getVector(0));

        solution.getVector(0) += 2*solution.getVector(1) - solution.getVector(2) + dt*dt*Ricker(iStep*dt,0.20,2000)*source;



        //Leap-frog scheme to compute the solution with mixed operator
        MltAddStiffness<2>(fespace,isotropic_material_no_lambda,-dt*dt,solution_mixed.getVector(1),0.0,solution_mixed.getVector(0));

       
        MltAddDivergence(fespace,fespace_pressure,-dt*dt,solution_mixed.getVector(1),0.0,pressure);

        // MltAddDivergence(fespace,-dt*dt,solution_mixed.getVector(1),0.0,pressure_tmp);
        //InterpolationDG<1,true>(aGLI,fespace_pressure,fespace,pressure_tmp,pressure);
 
        LAL::Solve(massMatrixPressure,pressure);

        //InterpolationDG<1,false>(aGLI,fespace_pressure,fespace,pressure,pressure_tmp);
        //TransposeMltAddDivergence(fespace,1.0,pressure_tmp,1.0,solution_mixed.getVector(0));

        TransposeMltAddDivergence(fespace,fespace_pressure,1.0,pressure,1.0,solution_mixed.getVector(0));    

        LAL::Solve(massMatrix,solution_mixed.getVector(0));

        solution_mixed.getVector(0) += 2*solution_mixed.getVector(1) - solution_mixed.getVector(2) + dt*dt*Ricker(iStep*dt,0.20,2000)*source;


     
        // Writing solution.
        if (iStep % OutputFreq == 0 && iStep>0)
        {
   
            // Writing solution
            ParallelWriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", fespace,solution.getVector(0));
            
            // Writing solution
            ParallelWriteVTK(FileName+"mixed."+std::to_string(NOutput) + ".vtk", fespace,solution_mixed.getVector(0));
            
            diff = solution.getVector(0) - solution_mixed.getVector(0);
            
            // Writing solution
            ParallelWriteVTK(FileName+"diff."+std::to_string(NOutput) + ".vtk", fespace,diff);

            // Writing solution
            ParallelWriteVTK(FileName+"pressure."+std::to_string(NOutput) + ".vtk", fespace_pressure,pressure,false);
            
            
            NOutput++;
        }
  
        // Swapping solutions.
        solution.Swap();
        solution_mixed.Swap();
    }
     
      ParallelWriteVTKTerminate();
     
      //Output of the time spent
      time(&end);
      std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    
}
 
