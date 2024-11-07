#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;

std::array<Real, 2> velocity_field(const RealVector & xyz)
                        {
                            std::array<Real, 2> A;

                            A[0]=1;
                            A[1]=0;

                            return A;
                        }
 

template<class FESpace>
void MltFriedrichs(FESpace & aFESpace,
                     Field::Vector<2,Field::Regularity::C0> & AIJ,
                     Real alpha,
                     LAL::Vector& U,
                     LAL::Vector& V)
{
        Core::MltAdd<FESpace,FESpace,true,false,false,
                   Core::QuadratureOpt::Default, 
                   1, //DimU
                   1, //DimV
                   1, //DimI
                   1, //DimJ
                   Field::Identity,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Identity,
                   Field::Vector<2,Field::Regularity::C0>,
                   Field::Identity>(aFESpace,aFESpace,Field::IdentityField,alpha*0.5,U,0.0,V,AIJ,Field::IdentityField); 
  
        Core::MltAdd<FESpace,FESpace,true,false,false,
                   Core::QuadratureOpt::Default,
                   1, //DimU
                   1, //DimV
                   1, //DimI
                   1, //DimJ
                   Field::Identity,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Gradient,
                   Field::Identity,
                   Field::Vector<2,Field::Regularity::C0>>(aFESpace,aFESpace,Field::IdentityField,-alpha*0.5,U,1.0,V,Field::IdentityField,AIJ); 
}

void flux_op(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 1>, 1> & A)
{
    // Define matrix A(n) 
    Real a_1, a_2;
    a_1 = 1.0;
    a_2 = 0.0;
    A[0][0] = a_1*n[0]+a_2*n[1];
}

void flux_pen(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 1>, 1> & A)
{
    A[0][0] = 1.0;
}

 


int main(int argc, char *argv[])
{
    std::cout << "Running Transport equation........................" << std::endl;
    
    std::cout << "Simulation of a transport equation in a square" << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Scaling<2> scaling(1,1);
    Square<Scaling<2>> fespace(4,4,25,25,scaling);

    auto FESpaceSkeleton = fespace.getSkeletonFESpace(PeriodicityType::XY);
  
    Index N = fespace.getNumDoFsDisc();
    
    std::cout << "The simulation has " << N << " (discontinuous) degree of freedom" << std::endl;
     
    Real dt = 0.5*0.25/100;  

    Field::Vector<2, Field::Regularity::C0> Aij(velocity_field);
    Field::Scalar velocity(1.0);
    
    LAL::SparseMatrix fluxMatrix;
    LAL::SparseMatrix jumpMatrix;
    LAL::SparseMatrix MassMatrixP;
    LAL::SparseMatrix MassMatrixM;
    LAL::FactorizationMatrix invMassMatrixP;
    
    LAL::Allocate(fluxMatrix,N,N);
    LAL::Allocate(jumpMatrix,N,N);

    CreateMassDG(fespace,velocity,MassMatrixP,1.0);
    CreateMassDG(fespace,velocity,MassMatrixM,1.0);

    AssembleCenteredFlux<1>(FESpaceSkeleton,flux_op,fluxMatrix,1.0);
    AssembleJumpFlux<1>(FESpaceSkeleton,flux_pen,jumpMatrix,0.5);
     
    LAL::Compress(fluxMatrix);
    LAL::Compress(jumpMatrix);
    LAL::Compress(MassMatrixM);
    LAL::Factorize(invMassMatrixP,MassMatrixP);
 
    Index N_iter;
    Real Residue;
    Real rho;
   
   /*
    rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&fespace,&velocity](LAL::Vector &u,LAL::Vector &v) {MltAddStiffness(fespace,velocity,1.0,u,0.0,v);},
                              N,N_iter,Residue);
    
    Real dt =  0.95 * (2.0 / sqrt(rho));
    
    std::cout << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;
    
    std::cout << "\t" << dt << std::endl;
    */

    //Compute the source
    LAL::Vector initial_solution;
    LAL::Vector tmp;
    LAL::Vector tmp_tmp;

    LAL::Allocate(initial_solution, N);
    LAL::Allocate(tmp, N);
    LAL::Allocate(tmp_tmp, N);

    Real kx = 2;
    
    for (Index i=0;i<fespace.getNumDoFsDisc();++i)
    {
        RealVector P;
        fespace.getDoFCoordinateDisc(i,P);
        initial_solution(i) = Gaussian(P[0],0.5,25)*dGaussian(P[1],0.5,25);//*cos(2*M_PI*kx*(P[0]));
    }


    LAL::VectorSequence solution(3);
    solution.Allocate(N);

    //Init first iterate from initial velocity
    // U0
    solution.getVector(2) = initial_solution;

    // U1 : U1 = U0 - dt*A*U0
    /*
    MltFriedrichs(fespace,Aij,1.0,initial_solution,tmp);
    tmp += fluxMatrix*initial_solution;
    tmp += jumpMatrix*initial_solution;
    LAL::Solve(massMatrix,tmp);
    solution.getVector(1) = initial_solution - dt*tmp;
    */
    solution.getVector(1) = initial_solution;

    Real T_end = 400;
    Real OutputDeltat = 1;
    
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput =0;
 
    std::string FileName = argv[1];
  
    //Set the time variable
    time_t start, end;
    time(&start);

 

    cout << "Leap-frog iterations starts..." << std::endl;
    for (Index iStep = 0; iStep*dt < T_end; iStep++)
    {

        // Writing solution.
        if (iStep % OutputFreq == 0)
        {
            cout << "Writting output " << NOutput << endl;
            // Writing solution
            WriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", fespace,solution.getVector(1),false);
           
            NOutput++;
        }

    // getVector(0) : U(n+1), getVector(1) : U(n), getVector(2) : U(n-1)

        //Leap-frog scheme
    /*
        // U(n+1) = -2*dt*A*U(n) with A without interface 
        MltFriedrichs(fespace,Aij,-2.0*dt,solution.getVector(1),solution.getVector(0));
    
        // U(n+1) = -2*dt*A*U(n) with A with interface 
        solution.getVector(0) += -2.0*dt*fluxMatrix*solution.getVector(1);

        // U(n+1) = +(M-dt*S)U(n-1)
        solution.getVector(0) += MassMatrixM*solution.getVector(2) - 2*dt*jumpMatrix*solution.getVector(2);
        
        // Compute U(n+1) = M^-1*U(n+1)
        LAL::Solve(invMassMatrixP,solution.getVector(0));
    */

  
        //Euler
        // U(n+1) = -dt*A*U(n) 
   /*    MltFriedrichs(fespace,Aij,-dt,solution.getVector(1),solution.getVector(0));
        solution.getVector(0) += -dt*fluxMatrix*solution.getVector(1);

         // U(n+1) = +(M-dt*S)U(n)
        solution.getVector(0) += MassMatrixM*solution.getVector(1) - dt*jumpMatrix*solution.getVector(1);

        // Compute U(n+1) = M^-1*U(n+1)
        LAL::Solve(invMassMatrixP,solution.getVector(0));
 */
 
 
        //RK2
        //tmp = -dt*M^-1(A+S)*U(n) 
        MltFriedrichs(fespace,Aij,-dt,solution.getVector(1),tmp);
        tmp -= dt*fluxMatrix*solution.getVector(1);
        tmp -= dt*jumpMatrix*solution.getVector(1);
        LAL::Solve(invMassMatrixP,tmp);
        
        //U(n+1) = + 0.5*dt*dt*(M^-1(A+S))^2*U(n)
        MltFriedrichs(fespace,Aij,-dt*0.5,tmp,solution.getVector(0));
        solution.getVector(0) -= dt*0.5*fluxMatrix*tmp;
        solution.getVector(0) -= dt*0.5*jumpMatrix*tmp;

        LAL::Solve(invMassMatrixP,solution.getVector(0));

        // Correction 
        solution.getVector(0) += tmp;
        solution.getVector(0) += solution.getVector(1);
 
/*
        //RK3
        //tmp = -dt*M^-1(A+S)*U(n) 
        MltFriedrichs(fespace,Aij,-dt,solution.getVector(1),tmp);
        tmp -= dt*fluxMatrix*solution.getVector(1);
        tmp -= dt*jumpMatrix*solution.getVector(1);
        LAL::Solve(invMassMatrixP,tmp);
        
        //tmp_tmp = + -dt*dt*(M^-1(A+S))^2*U(n)
        MltFriedrichs(fespace,Aij,-dt,tmp,tmp_tmp);
        tmp_tmp -= dt*fluxMatrix*tmp;
        tmp_tmp -= dt*jumpMatrix*tmp;
        LAL::Solve(invMassMatrixP,tmp_tmp);

        //U(n+1)  = dt*dt*dt*(M^-1(A+S))^3*U(n)/3
        MltFriedrichs(fespace,Aij,-dt*(1.0/6.0),tmp_tmp,solution.getVector(0));
        solution.getVector(0) -= dt*(1.0/6.0)*fluxMatrix*tmp_tmp;
        solution.getVector(0) -= dt*(1.0/6.0)*jumpMatrix*tmp_tmp;
        LAL::Solve(invMassMatrixP,solution.getVector(0));

        // Correction
        solution.getVector(0) += tmp;
        solution.getVector(0) += tmp_tmp*0.5;
        solution.getVector(0) += solution.getVector(1);

*/
        // Swapping solutions.
        solution.Swap();
    }
     
     
     
     //Output of the time spent
     time(&end);
     std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
    
    
}
 
