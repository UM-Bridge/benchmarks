#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;

const Real Thickness = 0.05;

// C^1 function for the transformation
RealVector mapping_cube2shell(const RealVector& p)
{
    RealVector def;
    
    def[0]=p[0]-0.5;
    def[1]=p[1]-0.5;
    def[2]=Thickness*p[2];
    
    def[2]+=sqrt(1-def[0]*def[0]);

    return def;
}

// C^0 function for the gradient of the transformation
std::array<std::array<Real,3>,3> mapping_gradient_cube2shell(const RealVector& p)
{
    RealVector def;
    std::array<std::array<Real,3>,3> gradDef;
    
    def[0]=p[0]-0.5;
 
    gradDef[0][0]=1.0; gradDef[0][1]=0.0; gradDef[0][2]=0.0;
    gradDef[1][0]=0.0; gradDef[1][1]=1.0; gradDef[1][2]=0.0;
                       gradDef[2][1]=0  ; gradDef[2][2]=Thickness;
    
    gradDef[2][0]=-Thickness*def[0]/sqrt(1-def[0]*def[0]);
    
    return gradDef;
}




int main(int argc, char *argv[])
{
    std::cout << "Running ElasticWave_Plate........................" << std::endl;
    
    std::cout << "Simulation of an isotropic elastic wave in a deformed plate." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    
    Diffeomorphism<3> cube2shell(mapping_cube2shell,mapping_gradient_cube2shell);
    
    Cube<Diffeomorphism<3>> fespace(6,6,4,40,40,1,cube2shell);
    
    Index N = 3*fespace.getNumDoFs();
    
    std::cout << "The simulation has " << N << " degree of freedom" << std::endl;

    Field::Scalar density(1.0);
    Field::IsotropicTensor isotropic_material(2.0,1.0);
    
    LAL::DiagonalMatrix massMatrix;
    CreateMass<3>(fespace,density,massMatrix,1.0);
    
    Index N_iter;
    Real Residue;
    Real rho;
   
    rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&fespace, &isotropic_material](LAL::Vector &u,LAL::Vector &v) {MltAddStiffness<3>(fespace,isotropic_material,1.0,u,0.0,v);},
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
    
        RealVector Center = { 0.0, 0.0, 1+Thickness*0.5};
        
        source(3*i) = CInfty0(P,Center,0.05,1);
    }
      
   
    
    LAL::VectorSequence solution(3);
    solution.Allocate(N);
    
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
        MltAddStiffness<3>(fespace,isotropic_material,-dt*dt,solution.getVector(1),0.0,solution.getVector(0));

        LAL::Solve(massMatrix,solution.getVector(0));

        solution.getVector(0) += 2*solution.getVector(1)-solution.getVector(2)+dt*dt*Ricker(iStep*dt,0.20,2000)*source;
  
  
        // Writing solution.
        if (iStep % OutputFreq == 0 && iStep>0)
        {
            // Writing solution
            ParallelWriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", fespace,solution.getVector(0));
        
            NOutput++;
        }
  
        // Swapping solutions.
        solution.Swap();
    }
     
    ParallelWriteVTKTerminate();
     
    //Output of the time spent
    time(&end);
    std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    
}
 
