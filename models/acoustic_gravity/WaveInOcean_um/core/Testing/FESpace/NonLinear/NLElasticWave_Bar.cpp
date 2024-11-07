#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  

void mat_vec(const RealVector & xyz,
                        const std::array<std::array<Real, 3>, 3> & A,
                        std::array<Real, 9> & b)
                        {
                            b[0] = A[0][0];
                            b[1] = A[0][1];
                            b[2] = A[0][2];
                            b[3] = A[1][0];
                            b[4] = A[1][1];
                            b[5] = A[1][2];
                            b[6] = A[2][0];
                            b[7] = A[2][1];
                            b[8] = A[2][2];
                 
                        }

void adj_mat_vec(const RealVector & xyz,
                            const std::array<Real, 9> & b,
                            std::array<std::array<Real, 3>, 3> & A )
                        {
   
                            A[0][0] = b[0];
                            A[0][1] = b[1];
                            A[0][2] = b[2];
                            A[1][0] = b[3];
                            A[1][1] = b[4];
                            A[1][2] = b[5];
                            A[2][0] = b[6];
                            A[2][1] = b[7];
                            A[2][2] = b[8];
  
                        }


class FollowingPressure {

    public:

        static const Field::Regularity Regularity = Field::Regularity::Constant;

        void Eval(const std::array<Real, 9>& U,
                    std::array<Real, 3>& V, 
                    const RealVector & n, 
                    const RealVector & xyz, 
                    const Index & Label) const
            {
                    std::array<std::array<Real,3>, 3> F; 
                    std::array<std::array<Real,3>, 3> invFT; 
                    Real J;

                    //Compute F = gradU + I
                    F[0][0] = U[0]+1.0;
                    F[0][1] = U[1];
                    F[0][2] = U[2];
                    F[1][0] = U[3];
                    F[1][1] = U[4]+1.0;
                    F[1][2] = U[5];
                    F[2][0] = U[6];
                    F[2][1] = U[7];
                    F[2][2] = U[8]+1.0;

                    //Compute the transpose inverse of F (not sure is the optimal way)
                    ArrayAlgebra::CoMat(F,invFT);
                    J = ArrayAlgebra::Det(F);
                    ArrayAlgebra::Scale(1.0/J,invFT);

            
                    V = n;

                    
                    ArrayAlgebra::MatMlt<3>(invFT,V);               
            }

};


    class Determinant {
        
    public:
        
        static const Field::Regularity Regularity = Field::Regularity::Constant;

        

        // ---------------------------------------------------------------------------------//
        Determinant();


        // ---------------------------------------------------------------------------------//
        void Eval(const std::array<Real, 9>& U,
                                Real& V, 
                                const RealVector & xyz, 
                            const Index & Label) const
        {
            std::array<std::array<Real,3>, 3> F; 
            
            Real J;

            F[0][0] = U[0]+1.0;
            F[0][1] = U[1];
            F[0][2] = U[2];
            F[1][0] = U[3];
            F[1][1] = U[4]+1.0;
            F[1][2] = U[5];
            F[2][0] = U[6];
            F[2][1] = U[7];
            F[2][2] = U[8]+1.0;

            V = ArrayAlgebra::Det(F);
            
        }                  
    };




int main(int argc, char *argv[])
{
    std::cout << "Running NLElasticWave_Bar........................" << std::endl;
    
    std::cout << "Simulation of a non linear deformation of a bar." << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Scaling<3> scaling(1.0,0.1,0.1);
    Cube<Scaling<3>> fespace(4,4,4,10,2,2,scaling);
    auto fespaceLeft = fespace.getBndyFESpace(4);
    auto fespaceRight = fespace.getBndyFESpace(2);
    auto fespaceBottom = fespace.getBndyFESpace(0);
    
    Index N = fespace.getNumDoFs();
    Index N_disc = fespace.getNumDoFsDisc();
      
    std::string FileName = argv[1];
  
    std::cout << "The simulation has " << N << " scalar degree of freedom" << std::endl;
    
    Field::Scalar density(1);
    
    Field::IsotropicTensor isotropic_material(2.0,1.0);

    Field::NeoHookeanTensor nl_material(1.0,1.0);

    Field::ThirdOrderTensor<3,3,9> Voigt_Notation(mat_vec,adj_mat_vec);
    
    FollowingPressure aFollowingPressure;
 
    //Dirichlet boundary condition
    std::vector<Index> DirichletDoFs;
    std::vector<Index> ScalarDirichletDoFs;

    fespaceLeft.setAsFESpace();
    fespaceRight.setAsFESpace();

    for (Index i=0;i<fespaceLeft.getNumDoFs();++i)
    {
        Index dof = fespaceLeft._glob_bdny2glob_cube(i);
    
        DirichletDoFs.push_back(3*dof);
        DirichletDoFs.push_back(3*dof+1);
        DirichletDoFs.push_back(3*dof+2);

        ScalarDirichletDoFs.push_back(dof);
    }

    for (Index i=0;i<fespaceRight.getNumDoFs();++i)
    {
        Index dof = fespaceRight._glob_bdny2glob_cube(i);
    
        DirichletDoFs.push_back(3*dof);
        DirichletDoFs.push_back(3*dof+1);
        DirichletDoFs.push_back(3*dof+2);

        ScalarDirichletDoFs.push_back(dof);
    }




    LAL::DiagonalMatrix massMatrix;
    CreateMass<3>(fespace,density,massMatrix,1.0);
     
    LAL::DiagonalMatrix massMatrixF;
    CreateMassDG<9>(fespace,density,massMatrixF,1.0);
    
    //Compute the source////////////////////
    LAL::Vector source;
    LAL::Allocate(source, 3*N);
    
    for (Index i=0;i<N;++i)
    {
        source(3*i+2) = -1e-5; // Gravity
    }
     ///////////////////////////////////


    LAL::Vector U_tmp;
    LAL::Vector V_tmp;
    LAL::Vector W_tmp;
    LAL::Vector F_tmp;
    LAL::Vector U;



    LAL::Allocate(U_tmp, 3*N);
    LAL::Allocate(F_tmp, 9*N_disc);
    LAL::Allocate(V_tmp, N);
    LAL::Allocate(W_tmp, N_disc);
    LAL::Allocate(U, 3*N);
    
 

    //Set the time variable

    LAL::SparseMatrix stiff_matrix;
    LAL::FactorizationMatrix stiff_matrix_fac;
 
    std::cout << "Assembling preconditionner " << std::endl;
    CreateStiffness(fespace, Field::IdentityField , stiff_matrix, 1.0);
  
    //Modify the matrix to take into account the DBC
    std::cout << "Erasing lines and columns " << std::endl;
    LAL::EraseLines(stiff_matrix,ScalarDirichletDoFs);
    LAL::EraseColumns(stiff_matrix,ScalarDirichletDoFs);
    LAL::AddInteractions(stiff_matrix,ScalarDirichletDoFs,ScalarDirichletDoFs,1.0);

    std::cout << "Factorization begins " << std::endl;
    LAL::Factorize(stiff_matrix_fac,stiff_matrix);
    std::cout << "Factorization ends " << std::endl;

    //Modification of the source term (take into acount Dirchlet BS and multiplication by the mass matrix)
    source = massMatrix*source; 
    for (Index i=0;i<DirichletDoFs.size();++i) 
        source(DirichletDoFs[i]) = 0;


    Index NIter;
    Real J_value;

    Real aPressure = 0.0000001;

    for (Index i=0;i<3;++i)
    {

    LAL::SolveNLPCG(U,
                    [&fespace,&fespaceBottom,&Voigt_Notation,&nl_material,&W_tmp,&U_tmp,&F_tmp,&source,&massMatrixF,&aPressure,&aFollowingPressure](const LAL::Vector &v) -> Real { 
                        Real energy = 0.0;

                      /*  Core::MltAdd<decltype(fespace),decltype(fespace),true,true,false,
                            Core::QuadratureOpt::Default,
                            3,               //DimU
                            1,               //DimV
                            1,               //DimI
                            1,               //DimJ
                            Field::Identity,
                            Core::DiffOp::Gradient,
                            Core::DiffOp::Identity,
                            decltype(nl_material)>(fespace,fespace,Field::IdentityField,1.0,v,0.0,W_tmp,nl_material);*/

 
                        Core::MltAdd<decltype(fespace),decltype(fespace),true,true,false,
                            Core::QuadratureOpt::Default,
                            3,               //DimU
                            9,               //DimV
                            9,               //DimI
                            1,               //DimJ
                            Field::Identity,
                            Core::DiffOp::Gradient,
                            Core::DiffOp::Identity,
                            decltype(Voigt_Notation)>(fespace,fespace,Field::IdentityField,1.0,v,0.0,F_tmp,Voigt_Notation);

                        LAL::Solve(massMatrixF,F_tmp);

                        Core::MltAdd<decltype(fespace),decltype(fespace),true,false,false,
                            Core::QuadratureOpt::Default,
                            9,               //DimU
                            1,               //DimV
                            1,               //DimI
                            1,               //DimJ
                            Field::Identity,
                            Core::DiffOp::Identity,
                            Core::DiffOp::Identity,
                            decltype(nl_material)>(fespace,fespace,Field::IdentityField,1.0,F_tmp,0.0,W_tmp,nl_material);
 
                        energy += LAL::Sum(W_tmp);
/*
                        Core::MltAdd<decltype(fespace),decltype(fespace),true,false,false,
                            Core::QuadratureOpt::Default,
                            9,               //DimU
                            1,               //DimV
                            1,               //DimI
                            1,               //DimJ
                            Field::Identity,
                            Core::DiffOp::Identity,
                            Core::DiffOp::Identity,
                            decltype(nl_det)>(fespace,fespace,Field::IdentityField,aPressure,F_tmp,0.0,W_tmp,nl_det); */

                        energy +=  LAL::Sum(W_tmp);
                       
   

                        energy -= LAL::ScalarProduct(source,v);
                        return energy;},

                    [&fespace,&fespaceBottom,&aPressure,&nl_material,&DirichletDoFs,&source,&F_tmp,&massMatrixF,&Voigt_Notation,&aFollowingPressure](const LAL::Vector &u,LAL::Vector &v){ 
                       // MltAddStiffness<3>(fespace,nl_material,1.0,u,0.0,v); 

                        Core::MltAdd<decltype(fespace),decltype(fespace),true,true,false,
                            Core::QuadratureOpt::Default,
                            3,               //DimU
                            9,               //DimV
                            9,               //DimI
                            1,               //DimJ
                            Field::Identity,
                            Core::DiffOp::Gradient,
                            Core::DiffOp::Identity,
                            decltype(Voigt_Notation)>(fespace,fespace,Field::IdentityField,1.0,u,0.0,F_tmp,Voigt_Notation);

                        LAL::Solve(massMatrixF,F_tmp);

                        Core::MltAdd<decltype(fespace),decltype(fespace),true,false,true,
                            Core::QuadratureOpt::Default,
                            9,               //DimU
                            3,               //DimV
                            3,               //DimI
                            3,               //DimJ
                            Field::Identity,
                            Core::DiffOp::Identity,
                            Core::DiffOp::Gradient,
                            decltype(nl_material)>(fespace,fespace,Field::IdentityField,1.0,F_tmp,0.0,v,nl_material);


                        Core::MltAdd<decltype(fespaceBottom),decltype(fespaceBottom),true,false,true,
                            Core::QuadratureOpt::Default,
                            9,               //DimU
                            3,               //DimV
                            3,               //DimI
                            1,               //DimJ
                            Field::Identity,
                            Core::DiffOp::TraceNormal,
                            Core::DiffOp::Trace,
                            FollowingPressure>(fespaceBottom,fespaceBottom,Field::IdentityField,aPressure,F_tmp,1.0,v,aFollowingPressure);




                        for (Index i=0;i<DirichletDoFs.size();++i) v(DirichletDoFs[i]) = 0;

                        v = v - source;},
                    
                    [&DirichletDoFs,&stiff_matrix_fac,&N,&V_tmp](const LAL::Vector &u,LAL::Vector &v){
                        v = u; 
                        
                        for (Index i=0;i<DirichletDoFs.size();++i) 
                            v(DirichletDoFs[i]) = 0;

                        for (Index d=0;d<3;++d)
                        {
                            Real * data_v = LAL::getData(v);
                            Real * data_v_tmp = LAL::getData(V_tmp);

                            for (Index i=0;i<N;++i)
                                data_v_tmp[i] = data_v[3*i+d];

                            LAL::Solve(stiff_matrix_fac,V_tmp);

                            for (Index i=0;i<N;++i)
                                data_v[3*i+d] = data_v_tmp[i];
                        }

                        },
 
                    3*N,NIter,J_value,400,10e-6,2000,true);
    }

        
        //Non linear stiffness
        //Core::NonLinearStiffnessMltAdd<3,1>(fespace,nl_material,-1.0,displacement.getVector(1),velocity,0.0,displacement.getVector(0),dt);     
 
   
        // Writing solution.
        WriteVTK(FileName+"."+ ".vtk", fespace,U);  
        std::cout << "Done" << std::endl;

        
}
 
