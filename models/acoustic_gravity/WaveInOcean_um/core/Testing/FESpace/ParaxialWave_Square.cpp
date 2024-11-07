#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;

bool CenteredFlux = false;

// Paraxial wave square 15°
void advection_operator(const RealVector & xyz,
                        const std::array<std::array<Real, 2>, 2> & A,
                        std::array<Real, 2> & b
                        )
                        {
                            // 15° mobile v (to the left) :
                            b[0] = -A[0][0]-A[1][1]/sqrt(2.0);
                            b[1] = -A[0][1]/sqrt(2.0);
                            
                            // 15° solution u (to the top) :
                            /*
                            b[0] = A[0][1]-A[1][0]/sqrt(2.0);
                            b[1] = -A[0][0]/sqrt(2.0);
                            */
                        }

void adj_advection_operator(const RealVector & xyz,
                            const std::array<Real, 2> & b,
                            std::array<std::array<Real, 2>, 2> & A
                        )
                        {
                            // 15° mobile v (to the left) :
                            A[0][0] = -b[0];
                            A[0][1] = -b[1]/sqrt(2.0);
                            A[1][0] = 0.0;
                            A[1][1] = -b[0]/sqrt(2.0);

                            // 15° solution u (to the top) :
                            /*
                            A[0][0] = -b[1]/sqrt(2.0);
                            A[0][1] = b[0];
                            A[1][0] = -b[0]/sqrt(2.0);
                            A[1][1] = 0.0;
                            */
                        }

template<class FESpace>
void MltFriedrichs(FESpace & aFESpace,
                     Field::ThirdOrderTensor<2> & AIJ,
                     Real alpha,
                     LAL::Vector& U,
                     LAL::Vector& V)
{
        Core::MltAdd<FESpace,FESpace,true,false,false,
                   Core::QuadratureOpt::Default,
                   2, //DimU
                   2, //DimV
                   2, //DimI
                   1, //DimJ
                   Field::Identity,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Identity,
                   Field::ThirdOrderTensor<2>,
                   Field::Identity>(aFESpace,aFESpace,Field::IdentityField,alpha*0.5,U,0.0,V,AIJ,Field::IdentityField); 
  
        Core::MltAdd<FESpace,FESpace,true,false,false,
                   Core::QuadratureOpt::Default,
                   2, //DimU
                   2, //DimV
                   2, //DimI
                   1, //DimJ
                   Field::Identity,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Gradient,
                   Field::Identity,
                   Field::ThirdOrderTensor<2>>(aFESpace,aFESpace,Field::IdentityField,-alpha*0.5,U,1.0,V,Field::IdentityField,AIJ); 
}

// Paraxial wave square 15°
void flux_op(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 2>, 2> & A)
{
    // Define matrix A(n) 
    // 15° mobile v (to the left) :
    A[0][0] = -n[0];
    A[0][1] = -n[1]/sqrt(2.0);
    A[1][0] = -n[1]/sqrt(2.0);
    A[1][1] = 0.0;
    
    // 15° solution u (to the top) :
    /*
    A[0][0] = n[1];
    A[0][1] = -n[0]/sqrt(2.0);
    A[1][0] = -n[0]/sqrt(2.0);
    A[1][1] = 0.0;
    */

}

void flux_op_penalized(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 2>, 2> & A)
{
    // Define matrix A(n) 
    
    A[0][0] = abs(n[0])+0.0*abs(n[1])/sqrt(2.0);
    A[0][1] = 0.0;
    A[1][0] = 0.0;
    A[1][1] = 0.0*abs(n[1])/sqrt(2.0);
    

    // For any mesh : 
    /*  
    Real N = sqrt(n[0]*n[0]+2*n[1]*n[1]);
    Real X = n[1]*sqrt(2.0)/(abs(n[0])+N);
    Real X_2 = 1+X*X;
    
    A[0][0] = (N+abs(n[0])+X*X*(N-abs(n[0])))/(2*X_2);
    A[0][1] = -X*n[0]/X_2;
    A[1][0] = A[0][1];
    A[1][1] = (X*X*(N+abs(n[0]))+N-abs(n[0]))/(2*X_2);
    */
}

void bndy_decentered_matrix(std::array<std::array<Real, 2>, 2> & A,std::array<std::array<Real, 2>, 2> & B)
{
    Real tmp_0_0, tmp_0_1, tmp_1_0, tmp_1_1;
    tmp_0_0 = B[0][0]; tmp_0_1 = B[0][1];
    tmp_1_0 = B[1][0]; tmp_1_1 = B[1][1];
    B[0][0] = (A[0][0]-tmp_0_0)*(A[0][0]-tmp_0_0)+(A[1][0]-tmp_1_0)*(A[1][0]-tmp_1_0);
    B[0][1] = (A[0][0]-tmp_0_0)*(A[0][1]-tmp_0_1)+(A[1][0]-tmp_1_0)*(A[1][1]-tmp_1_1);
    B[1][1] = (A[0][1]-tmp_0_1)*(A[0][1]-tmp_0_1)+(A[1][1]-tmp_1_1)*(A[1][1]-tmp_1_1);
    B[1][0] = B[0][1];
}

void bndy_right(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 2>, 2> & B)
{
    std::array<std::array<Real, 2>, 2> A;

    // Define matrix B for a boundary on the right's domain 
    B[0][0] = 0.5*0.5;
    B[0][1] = 0.0*0.5;
    B[1][0] = 0.0*0.5;
    B[1][1] = 0.0*0.5;

    if(!CenteredFlux){
        flux_op(x,n,A);
        bndy_decentered_matrix(A,B);
    }
}

void bndy_left(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 2>, 2> & B)
{
    std::array<std::array<Real, 2>, 2> A;

    // Define matrix B for a boundary on the left's domain 
    B[0][0] = 1.0*0.5;
    B[0][1] = 0.0*0.5;
    B[1][0] = 0.0*0.5;
    B[1][1] = 0.0*0.5;

    if(!CenteredFlux){
        flux_op(x,n,A);
        bndy_decentered_matrix(A,B);
    }
}

void bndy_bottom(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 2>, 2> & B)
{
    std::array<std::array<Real, 2>, 2> A;

    // Define matrix B for a boundary on the bottom's domain 
    B[0][0] = 0.0*0.5;
    B[0][1] = 0.5/sqrt(2.0);
    B[1][0] = -0.5/sqrt(2.0);
    B[1][1] = 0.0*0.5;

    if(!CenteredFlux){
        flux_op(x,n,A);
        bndy_decentered_matrix(A,B);
    }
}

void bndy_top(const RealVector x,const std::array<Real, 2> & n,std::array<std::array<Real, 2>, 2> & B)
{   
    std::array<std::array<Real, 2>, 2> A;

    // Define matrix B for a boundary on the top's domain 
    B[0][0] = 0.0*0.5;
    B[0][1] = -0.5/sqrt(2.0);
    B[1][0] = 0.5/sqrt(2.0);
    B[1][1] = 0.0*0.5;

    if(!CenteredFlux){
        flux_op(x,n,A);
        bndy_decentered_matrix(A,B);
    }
}

int main(int argc, char *argv[])
{
    Real delta_DecenteredFlux = 0.2; // 0 = Centered Flux
    Real alpha_DecenteredFlux = delta_DecenteredFlux+0.5;

    std::cout << "Running ParaxialWave_Square........................" << std::endl;
    
    std::cout << "Simulation of a paraxial equation in a square (15° angle)." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Scaling<2> scaling(3,1);
    Square<Scaling<2>> fespace(8,8,30,10,scaling);

    //std::cout<<"Number of elements in total : "<<fespace.getNumElements()<<std::endl;
    
    auto FESpaceSkeleton = fespace.getSkeletonFESpace();

    auto FESpaceTraceBottom = fespace.getBndyFESpace(0);
    auto FESpaceTraceRight = fespace.getBndyFESpace(1);
    auto FESpaceTraceTop = fespace.getBndyFESpace(2);
    auto FESpaceTraceLeft = fespace.getBndyFESpace(3);

    Index N = 2*fespace.getNumDoFsDisc();

    std::cout << "The simulation has " << N << " (discontinuous) degree of freedom" << std::endl;
     
    Real dt = 0.00005;  

    Field::ThirdOrderTensor<2> Aij(advection_operator,adj_advection_operator);
    Field::Scalar velocity(1.0);
    
    LAL::SparseMatrix MassMatrixP;
    LAL::SparseMatrix MassMatrixM;
    LAL::SparseMatrix fluxMatrix;
    //LAL::SparseMatrix MassMatrixSource;
    LAL::FactorizationMatrix invMassMatrixP;
    //LAL::DiagonalMatrix massMatrix;
   
    LAL::Allocate(fluxMatrix,N,N);

    // To modify velocity c of media : velocity = 1/c

    CreateMassDG<2>(fespace,velocity,MassMatrixP,1);
    CreateMassDG<2>(fespace,velocity,MassMatrixM,1);
    //CreateMassDG<2>(fespace,velocity,MassMatrixSource,1.0);

    // Interface with centered flux : fluxMatrix = A 
    /*
    AssembleCenteredFlux<2>(FESpaceSkeleton,flux_op,fluxMatrix,1.0);
    std::cout<<"fluxMatrix for centered flux assembled"<<std::endl;

    // Boundary condition : centered flux
    // Assemble matrix : (M+dt*B) = MassMatrixP and (M-dt*B) = MassMatrixM

    AssembleBndyCenteredFlux<2>(FESpaceTraceBottom,bndy_bottom,MassMatrixP,dt);
    AssembleBndyCenteredFlux<2>(FESpaceTraceBottom,bndy_bottom,MassMatrixM,-dt);

    AssembleBndyCenteredFlux<2>(FESpaceTraceRight,bndy_right,MassMatrixP,dt);
    AssembleBndyCenteredFlux<2>(FESpaceTraceRight,bndy_right,MassMatrixM,-dt);

    AssembleBndyCenteredFlux<2>(FESpaceTraceTop,bndy_top,MassMatrixP,dt);
    AssembleBndyCenteredFlux<2>(FESpaceTraceTop,bndy_top,MassMatrixM,-dt);

    AssembleBndyCenteredFlux<2>(FESpaceTraceLeft,bndy_left,MassMatrixP,dt);
    AssembleBndyCenteredFlux<2>(FESpaceTraceLeft,bndy_left,MassMatrixM,-dt);

    std::cout<<"massMatrix for centered flux assembled"<<std::endl;

    if(!CenteredFlux){
        // Interface with decentered flux : PenalizedfluxMatrix = S_h
        AssembleJumpFlux<2>(FESpaceSkeleton,flux_op_penalized,MassMatrixP,dt*delta_DecenteredFlux);
        AssembleJumpFlux<2>(FESpaceSkeleton,flux_op_penalized,MassMatrixM,-dt*delta_DecenteredFlux);
        std::cout<<"fluxMatrix for decentered flux assembled"<<std::endl;
        
        // Boundary condition : centered flux
        // Assemble matrix : (M+dt*B+dt*S_Gamma) = MassMatrixP and (M-dt*B-dt*S_Gamma) = MassMatrixM

        AssembleBndyDecenteredFlux<2>(FESpaceTraceBottom,bndy_bottom,MassMatrixP,dt*alpha_DecenteredFlux);
        AssembleBndyDecenteredFlux<2>(FESpaceTraceBottom,bndy_bottom,MassMatrixM,-dt*alpha_DecenteredFlux);

        AssembleBndyDecenteredFlux<2>(FESpaceTraceRight,bndy_right,MassMatrixP,dt*alpha_DecenteredFlux);
        AssembleBndyDecenteredFlux<2>(FESpaceTraceRight,bndy_right,MassMatrixM,-dt*alpha_DecenteredFlux);

        AssembleBndyDecenteredFlux<2>(FESpaceTraceTop,bndy_top,MassMatrixP,dt*alpha_DecenteredFlux);
        AssembleBndyDecenteredFlux<2>(FESpaceTraceTop,bndy_top,MassMatrixM,-dt*alpha_DecenteredFlux);
        
        AssembleBndyDecenteredFlux<2>(FESpaceTraceLeft,bndy_left,MassMatrixP,dt*alpha_DecenteredFlux);
        AssembleBndyDecenteredFlux<2>(FESpaceTraceLeft,bndy_left,MassMatrixM,-dt*alpha_DecenteredFlux);

        std::cout<<"massMatrix for decentered flux assembled"<<std::endl;
    }

    LAL::Compress(fluxMatrix);
    LAL::Compress(MassMatrixM);
    LAL::Factorize(invMassMatrixP,MassMatrixP);

    //Compute the source
    LAL::Vector initial_solution;
    LAL::Vector tmp;
    LAL::Vector source;

    LAL::Allocate(initial_solution, N);
    LAL::Allocate(tmp, N);
    LAL::Allocate(source, N);

    Real kx = 10;
    
    for (Index i=0;i<fespace.getNumDoFsDisc();++i)
    {
        RealVector P;
        fespace.getDoFCoordinateDisc(i,P); 
        
        initial_solution(2*i) = Gaussian(P[0],1.5,150)*dGaussian(P[1],0.5,150)*cos(2*M_PI*kx*(P[0]));
        initial_solution(2*i+1) = - sqrt(2)*dGaussian(P[0],1.5,150)*Gaussian(P[1],0.5,150)*cos(2*M_PI*kx*(P[0])) 
                               - sqrt(2)*2*M_PI*kx*Gaussian(P[0],1.5,150)*Gaussian(P[1],0.5,150)*sin(2*M_PI*kx*(P[0]));
    }
    */

    /*
    for (Index i=0;i<fespace.getNumDoFsDisc();++i)
    {
        RealVector P;
        fespace.getDoFCoordinateDisc(i,P);
        
        source(2*i) = 0.5*Gaussian(P[0],1.5,150)*Gaussian(P[1],0.5,150)*cos(2*M_PI*kx*(P[0]));
    }
    */

    /*
    LAL::VectorSequence solution(3);
    solution.Allocate(N);

    //Init first iterate from initial velocity
    // U0 :
    solution.getVector(2) = initial_solution;

    // U1 : U1 = U0 
    solution.getVector(1) = initial_solution;

    Real T_end = 10.0;
    Real OutputDeltat = 0.1;
    
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput =0;
 
    std::string FileName = argv[1];
  
    // Set the time variable
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

        //Leap-frog scheme
        // getVector(0) : U(n+1), getVector(1) : U(n), getVector(2) : U(n-1)

        // Updating stiffness matrix : 
        // Corresponding to : U(n+1) = -2*dt*A*U(n) with A without interface 
        MltFriedrichs(fespace,Aij,-2.0*dt,solution.getVector(1),solution.getVector(0));
    
        // U(n+1) = -2*dt*A*U(n) with A with interface 
        solution.getVector(0) += -2.0*dt*fluxMatrix*solution.getVector(1);

        // U(n+1) = +(M-dt*B-dt*S)U(n-1)
        solution.getVector(0) += MassMatrixM*solution.getVector(2);

        // U(n+1) += dt*M*F
        //solution.getVector(0) += (dt*Ricker(iStep*dt,0.20,500))*MassMatrixSource*source;

        // Compute U(n+1) = +(M+dt*B+dt*S)^-1*U(n-1)
        LAL::Solve(invMassMatrixP,solution.getVector(0));

        // Swapping solutions.
        solution.Swap();
    }
     
    ParallelWriteVTKTerminate();

    //Output of the time spent
    time(&end);
    std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;

    */
}
 
