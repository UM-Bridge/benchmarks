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
    
    Square fespace(8,8,20,20);
    Index N = fespace.getNumDoFs();
    Index M = fespace.getNumDoFsDisc();


    //Create the subdomain where the data are given
    auto fespace_right = fespace.getSubDomainFESpace([](RealVector xyz){
            if (xyz[0] > 0.9)
                return true;
            else return false;
    });

    fespace_right.setAsSubSpace(); // Non optimal version -> auxiliary unknown stored everywhere

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
    
    //Operators for the PML
    Real xi = 50.0;  
    RealMatrix2x2 dissip_mat;
    RealMatrix2x2 coupling_mat;

    coupling_mat[0][0] = -xi;  coupling_mat[0][1] = 0.0;
    coupling_mat[1][0] = 0.0;  coupling_mat[1][1] = +xi;

    dissip_mat[0][0] = xi;   dissip_mat[0][1] = 0.0;
    dissip_mat[1][0] = 0.0;  dissip_mat[1][1] = 0.0;

    
    Field::Scalar dissip_coeff(xi);
    Field::Matrix<2,2> dissip(dissip_mat);
    Field::Matrix<2,2> coupling(coupling_mat);

    LAL::DiagonalMatrix DissipMatrix;
    LAL::DiagonalMatrix MassMatrix_Disc_Dissip;
    LAL::DiagonalMatrix MassMatrix_Disc;

    CreateMass(fespace_right,dissip_coeff,DissipMatrix,1.0);
    CreateMassDG<2>(fespace_right,dissip,MassMatrix_Disc_Dissip,1.0); 
    CreateMassDG<2>(fespace,density,MassMatrix_Disc,1.0); //non optimal -> computed everywhere

    LAL::DiagonalMatrix Matrix_p;
    LAL::DiagonalMatrix Matrix_m;

    Matrix_p = massMatrix + 0.5*dt*DissipMatrix;
    Matrix_m = massMatrix - 0.5*dt*DissipMatrix;

    LAL::DiagonalMatrix MassMatrix_Disc_p;
    LAL::DiagonalMatrix MassMatrix_Disc_m;

    MassMatrix_Disc_p = MassMatrix_Disc + 0.5*dt*MassMatrix_Disc_Dissip;
    MassMatrix_Disc_m = MassMatrix_Disc - 0.5*dt*MassMatrix_Disc_Dissip;

    //Compute the source
    LAL::Vector source;
    LAL::Allocate(source, N);

  
    
    for (Index i=0;i<fespace.getNumDoFs();++i)
    {
        RealVector P;
        fespace.getDoFCoordinate(i,P);
    
        RealVector Center = { 0.5, 0.5, 0.};
        source(i) = CInfty0(P,Center,0.1,1);
    }


    LAL::VectorSequence solution(3);
    solution.Allocate(N);

    LAL::VectorSequence auxilliary_unknown(2);
    auxilliary_unknown.Allocate(2*M);

    LAL::Vector tmp; 
    
    Real T_end = 4.0;
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

        TransposeMltAddGradient(fespace_right,-dt*dt,auxilliary_unknown.getVector(1),1.0,solution.getVector(0));

        solution.getVector(0) -= Matrix_m*solution.getVector(2);
        solution.getVector(0) += massMatrix*(2*solution.getVector(1) + dt*dt*Ricker(iStep*dt,0.20,1000)*source);

        LAL::Solve(Matrix_p,solution.getVector(0));


        //update the auxiliary variable
        tmp = dt*0.5*(solution.getVector(0) + solution.getVector(1));

        Core::MltAdd<decltype(fespace_right),decltype(fespace_right),true,true,false,
                   Core::QuadratureOpt::Default,
                   1,                   //DimU
                   2,                   //DimV
                   1,                   //DimI
                   2,                   //DimJ
                   Field::Matrix<2,2>,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Identity>(fespace_right,fespace_right,coupling,1.0,tmp,0.0,auxilliary_unknown.getVector(0));

        
    
        auxilliary_unknown.getVector(0) += MassMatrix_Disc_m*auxilliary_unknown.getVector(1);

        LAL::Solve(MassMatrix_Disc_p,auxilliary_unknown.getVector(0));

  
  
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
        auxilliary_unknown.Swap();
    }
     
    ParallelWriteVTKTerminate();
    ParallelWriteVIZIRTerminate();
     
     //Output of the time spent
     time(&end);
     std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    
    
}
 
