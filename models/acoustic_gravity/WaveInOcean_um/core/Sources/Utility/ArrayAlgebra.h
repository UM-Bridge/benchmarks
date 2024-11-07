#pragma once

// -----------------------------------------------------------------------------------------//
#include <vector>
#include <array>
#include <math.h>
// -----------------------------------------------------------------------------------------//



// -----------------------------------------------------------------------------------------//
namespace OndoMathX {

    namespace ArrayAlgebra {
        
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Computing cofactor matrix.
         \param Mat is a matrix used to compute the corresponding cofactor matrix.
         \param CoMat is the corresponding cofactor matrix.
         */
        template<
        Index Dim
        > inline void CoMat(const std::array<std::array<Real, Dim>, Dim> &Mat, std::array<std::array<Real, Dim>, Dim>& CoMat);
        
        
        
        /*! \brief Specialization of the templated CoMat() function int the 2D case.
         */
        template <
        // Template specialization.
        > inline void CoMat<2>(const std::array<std::array<Real, 2>, 2> &Mat, std::array<std::array<Real, 2>, 2>& CoMat)
        {
            CoMat[0][0] = Mat[1][1];
            CoMat[1][0] = -Mat[0][1];
            CoMat[0][1] = -Mat[1][0];
            CoMat[1][1] = Mat[0][0];
        }
        
        
        /*! \brief Specialization of the templated CoMat() function int the 3D case.
         */
        template <
        // Template specialization.
        > inline void CoMat<3>(const std::array<std::array<Real, 3>, 3> &Mat, std::array<std::array<Real, 3>, 3>& CoMat)
        {
            CoMat[0][0] =  Mat[1][1]*Mat[2][2]-Mat[2][1]*Mat[1][2];
            CoMat[0][1] =  -Mat[1][0]*Mat[2][2]+Mat[2][0]*Mat[1][2];
            CoMat[0][2] =  Mat[1][0]*Mat[2][1]-Mat[2][0]*Mat[1][1];
            
            CoMat[1][0] =  -Mat[0][1]*Mat[2][2]+Mat[2][1]*Mat[0][2];
            CoMat[1][1] =  Mat[0][0]*Mat[2][2]-Mat[2][0]*Mat[0][2];
            CoMat[1][2] =  -Mat[0][0]*Mat[2][1]+Mat[2][0]*Mat[0][1];
            
            CoMat[2][0] =  Mat[0][1]*Mat[1][2]-Mat[1][1]*Mat[0][2];
            CoMat[2][1] =  -Mat[0][0]*Mat[1][2]+Mat[1][0]*Mat[0][2];
            CoMat[2][2] =  Mat[0][0]*Mat[1][1]-Mat[1][0]*Mat[0][1];
        }
        // ---------------------------------------------------------------------------------//
        

        // ---------------------------------------------------------------------------------//
        /*! \brief Computing the determinant of a matrix.
         \param Mat is a matrix used to compute the corresponding graddefobian.
         */
        template<
        Index Dim
        > inline Real Det(const std::array<std::array<Real, Dim>, Dim> &Mat);
       
        
        /*! \brief Specialization of the templated Det() function int the 2D case.
         */
        template <
        // Template specialization.
        > inline Real Det<2>(const std::array<std::array<Real, 2>, 2> &Mat)
        {
            return Mat[0][0] * Mat[1][1] - Mat[0][1] * Mat[1][0];
        }
        
        /*! \brief Specialization of the templated Det() function int the 3D case.
         */
        template <
        // Template specialization.
        > inline Real Det<3>(const std::array<std::array<Real, 3>, 3> &Mat)
        {
            return (Mat[0][0] * (Mat[1][1] * Mat[2][2] - Mat[1][2] * Mat[2][1]) -
                    Mat[0][1] * (Mat[1][0] * Mat[2][2] - Mat[1][2] * Mat[2][0]) +
                    Mat[0][2] * (Mat[1][0] * Mat[2][1] - Mat[1][1] * Mat[2][0]));
        }
        
        // ---------------------------------------------------------------------------------//
        
        
        
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Computing cofactor matrix.
         \param Mat is a matrix used to compute the corresponding cofactor matrix.
         \param CoMat is the corresponding cofactor matrix.
         */
        template<
        Index Dim
        > inline void Inv(const std::array<std::array<Real, Dim>, Dim> &Mat, std::array<std::array<Real, Dim>, Dim>& InvMat);
        
        
        
        /*! \brief Specialization of the templated CoMat() function int the 2D case.
         */
        template <
        // Template specialization.
        > inline void Inv<2>(const std::array<std::array<Real, 2>, 2> &Mat, std::array<std::array<Real, 2>, 2>& InvMat)
        {
            Real inv_det =  1.0/Det<2>(Mat);
            
            InvMat[0][0] = inv_det*Mat[1][1];
            InvMat[1][0] = -inv_det*Mat[1][0];
            InvMat[0][1] = -inv_det*Mat[0][1];
            InvMat[1][1] = inv_det*Mat[0][0];
  
        }
        
        
        /*! \brief Specialization of the templated CoMat() function int the 3D case.
         */
        template <
        // Template specialization.
        > inline void Inv<3>(const std::array<std::array<Real, 3>, 3> &Mat, std::array<std::array<Real, 3>, 3>& InvMat)
        {
            Real inv_det =  1.0/Det<3>(Mat);
            
            InvMat[0][0] =  inv_det*(Mat[1][1]*Mat[2][2]-Mat[2][1]*Mat[1][2]);
            InvMat[1][0] =  inv_det*(-Mat[1][0]*Mat[2][2]+Mat[2][0]*Mat[1][2]);
            InvMat[2][0] =  inv_det*(Mat[1][0]*Mat[2][1]-Mat[2][0]*Mat[1][1]);
            
            InvMat[0][1] =  inv_det*(-Mat[0][1]*Mat[2][2]+Mat[2][1]*Mat[0][2]);
            InvMat[1][1] =  inv_det*(Mat[0][0]*Mat[2][2]-Mat[2][0]*Mat[0][2]);
            InvMat[2][1] =  inv_det*(-Mat[0][0]*Mat[2][1]+Mat[2][0]*Mat[0][1]);
            
            InvMat[0][2] =  inv_det*(Mat[0][1]*Mat[1][2]-Mat[1][1]*Mat[0][2]);
            InvMat[1][2] =  inv_det*(-Mat[0][0]*Mat[1][2]+Mat[1][0]*Mat[0][2]);
            InvMat[2][2] =  inv_det*(Mat[0][0]*Mat[1][1]-Mat[1][0]*Mat[0][1]);
        }
        // ---------------------------------------------------------------------------------//
        
        

                
        /*! \brief Add a matrix into another
         */
        template<
        Index DimI,
        Index DimJ
        > inline void MatAdd(const std::array<std::array<Real, DimI>, DimJ> &A, std::array<std::array<Real, DimI>, DimJ>& B)
        {
            for (Index i=0;i<DimI;++i)
                for (Index j=0;j<DimJ;++j) 
                    B[i][j]+=A[i][j];           
        }

        template<
        Index DimI
        > inline void MatAdd(const std::array<Real, DimI>&A, std::array<Real, DimI>& B)
        {
            for (Index i=0;i<DimI;++i)
                    B[i]+=A[i];           
        }
            
        
        inline void MatAdd(const Real &A, Real& B)
        {
            B+=A;           
        }
            
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Performing the matrix times vector operation, results stored in-place.
         */
        template<
        Index Dim
        > inline void MatMlt(const std::array<std::array<Real, Dim>, Dim> &A,
                      std::array<Real, Dim>& U);
        
        /*! \brief Spcialization of the MatMlt() function in the 1D case.
         */
        template <
        // Template specialization.
        > inline void MatMlt<1>(const std::array<std::array<Real, 1>, 1> &A,
                     std::array<Real, 1>& U)
        {
            U[0] = A[0][0] * U[0];
        }
        
        /*! \brief Spcialization of the MatMlt() function in the 2D case.
         */
        template <
        // Template specialization.
        > inline void MatMlt<2>(const std::array<std::array<Real, 2>, 2> &A,
                         std::array<Real, 2>& U)
        {
            Real Tmp1 = U[0];
            Real Tmp2 = U[1];
            
            U[0] = A[0][0] * Tmp1 + A[0][1] * Tmp2;
            U[1] = A[1][0] * Tmp1 + A[1][1] * Tmp2;
           
        }
        
        
        /*! \brief Spcialization of the MatMlt() function in the 3D case.
         */
        template <
        // Template specialization.
        > inline void MatMlt<3>(const std::array<std::array<Real, 3>, 3> &A,
                         std::array<Real, 3>& U)
        {
          Real Tmp1 = U[0];
          Real Tmp2 = U[1];
          Real Tmp3 = U[2];
          
          U[0] = A[0][0] * Tmp1 + A[0][1] * Tmp2 + A[0][2] * Tmp3;
          U[1] = A[1][0] * Tmp1 + A[1][1] * Tmp2 + A[1][2] * Tmp3;
          U[2] = A[2][0] * Tmp1 + A[2][1] * Tmp2 + A[2][2] * Tmp3;
        }
        // ---------------------------------------------------------------------------------//
        
    // ---------------------------------------------------------------------------------//
 /*! \brief Performing the matrix NxN times matrix NxN operation, results stored in another matrix.
  */
 template<
 Index DimI,
 Index DimJ,
 Index DimK
 > inline void MatMlt(const std::array<std::array<Real, DimJ>, DimI> &A,
               const std::array<std::array<Real, DimK>, DimJ> &B,
               std::array<std::array<Real, DimK>, DimI> &AB)
               {   
                    for (Index i=0;i<DimI;++i)
                        for (Index k=0;k<DimK;++k)
                        {
                            AB[i][k] = 0.0; 
                            for (Index j=0;j<DimJ;++j) AB[i][k]+=A[i][j]*B[j][k];
                        }
               }
 
 
 
 /*! \brief Spcialization of the MatMlt() function in the 2D case.
  */
 template <
 // Template specialization.
 > inline void MatMlt<2,2,2>(const std::array<std::array<Real, 2>, 2> &A,
                  const std::array<std::array<Real, 2>, 2> &B,
                  std::array<std::array<Real, 2>, 2> &AB)
 {
     AB[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0];
     AB[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1];
     AB[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0];
     AB[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1];
 }
 
 
 /*! \brief Spcialization of the MatMlt() function in the 3D case.
  */
 template <
 // Template specialization.
 > inline void MatMlt<3,3,3>(const std::array<std::array<Real, 3>, 3> &A,
                      const std::array<std::array<Real, 3>, 3> &B,
                      std::array<std::array<Real, 3>, 3> &AB)
 {
     AB[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0];
     AB[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1];
     AB[0][2]=A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2];

     AB[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0];
     AB[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1];
     AB[1][2]=A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2];

     AB[2][0]=A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0];
     AB[2][1]=A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1];
     AB[2][2]=A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2];
 }
 // ---------------------------------------------------------------------------------//
     
    

     template<
 Index Dim
 > inline void MatMlt(const std::array<std::array<Real, Dim>, Dim> &A,
                     std::array<std::array<Real, Dim>, Dim> &B);
 
 
 
 /*! \brief Spcialization of the MatMlt() function in the 2D case.
  */
 template <
 // Template specialization.
 > inline void MatMlt<2>(const std::array<std::array<Real, 2>, 2> &A,
                        std::array<std::array<Real, 2>, 2> &B)
 {
    Real b00 = B[0][0];
    Real b01 = B[0][1];
    Real b10 = B[1][0];
    Real b11 = B[1][1];

    B[0][0]=A[0][0]*b00+A[0][1]*b10;
    B[0][1]=A[0][0]*b01+A[0][1]*b11;
    B[1][0]=A[1][0]*b00+A[1][1]*b10;
    B[1][1]=A[1][0]*b01+A[1][1]*b11;
 }
 
 
 /*! \brief Spcialization of the MatMlt() function in the 3D case.
  */
 template <
 // Template specialization.
 > inline void MatMlt<3>(const std::array<std::array<Real, 3>, 3> &A,
                        std::array<std::array<Real, 3>, 3> &B)
 {
    std::array<std::array<Real, 3>, 3> C = B;

    B[0][0]=A[0][0]*C[0][0]+A[0][1]*C[1][0]+A[0][2]*C[2][0];
    B[0][1]=A[0][0]*C[0][1]+A[0][1]*C[1][1]+A[0][2]*C[2][1];
    B[0][2]=A[0][0]*C[0][2]+A[0][1]*C[1][2]+A[0][2]*C[2][2];

    B[1][0]=A[1][0]*C[0][0]+A[1][1]*C[1][0]+A[1][2]*C[2][0];
    B[1][1]=A[1][0]*C[0][1]+A[1][1]*C[1][1]+A[1][2]*C[2][1];
    B[1][2]=A[1][0]*C[0][2]+A[1][1]*C[1][2]+A[1][2]*C[2][2];

    B[2][0]=A[2][0]*C[0][0]+A[2][1]*C[1][0]+A[2][2]*C[2][0];
    B[2][1]=A[2][0]*C[0][1]+A[2][1]*C[1][1]+A[2][2]*C[2][1];
    B[2][2]=A[2][0]*C[0][2]+A[2][1]*C[1][2]+A[2][2]*C[2][2];
 }
 // ---------------------------------------------------------------------------------//


    // ---------------------------------------------------------------------------------//
    /*! \brief Performing the matrix times vector operation, results stored in another vector.
     */
    template<
    Index Dim
    > inline void MatMlt(const std::array<std::array<Real, Dim>, Dim> &A,
                  const std::array<Real, Dim>& U, std::array<Real, Dim>& V);
    
    
    
    /*! \brief Spcialization of the MatMlt() function in the 2D case.
     */
    template <
    // Template specialization.
    > inline void MatMlt<2>(const std::array<std::array<Real, 2>, 2> &A,
                     const std::array<Real, 2>& U, std::array<Real, 2>& V)
    {
        V[0] = A[0][0] * U[0] + A[0][1] * U[1];
        V[1] = A[1][0] * U[0] + A[1][1] * U[1];
    }
    
    
    /*! \brief Spcialization of the MatMlt() function in the 3D case.
     */
    template <
    // Template specialization.
    > inline void MatMlt<3>(const std::array<std::array<Real, 3>, 3> &A,
                     const std::array<Real, 3>& U, std::array<Real, 3>& V)
    {
        V[0] = A[0][0] * U[0] + A[0][1] * U[1] + A[0][2] * U[2];
        V[1] = A[1][0] * U[0] + A[1][1] * U[1] + A[1][2] * U[2];
        V[2] = A[2][0] * U[0] + A[2][1] * U[1] + A[2][2] * U[2];
    }
    // ---------------------------------------------------------------------------------//
        
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Performing the transpose matrix times vector operation, results stored in-place.
         */
        template<
        Index Dim
        > inline void TransposeMatMlt(const std::array<std::array<Real, Dim>, Dim> &A,
                               std::array<Real, Dim>& U);
        
        
        
        /*! \brief Spcialization of the TransposeMatMlt() function in the 2D case.
         */
        template <
        // Template specialization.
        > inline void TransposeMatMlt<2>(const std::array<std::array<Real, 2>, 2> &A,
                                  std::array<Real, 2>& U)
        {
            Real Tmp1 = U[0];
            Real Tmp2 = U[1];
            
            U[0] = A[0][0] * Tmp1 + A[1][0] * Tmp2;
            U[1] = A[0][1] * Tmp1 + A[1][1] * Tmp2;
        }
        
        /*! \brief Spcialization of the TransposeMatMlt() function in the 3D case.
         */
        template <
        // Template specialization.
        > inline void TransposeMatMlt<3>(const std::array<std::array<Real, 3>, 3> &A,
                                  std::array<Real, 3>& U)
        {
            Real Tmp1 = U[0];
            Real Tmp2 = U[1];
            Real Tmp3 = U[2];
            
            U[0] = A[0][0] * Tmp1 + A[1][0] * Tmp2 + A[2][0] * Tmp3;
            U[1] = A[0][1] * Tmp1 + A[1][1] * Tmp2 + A[2][1] * Tmp3;
            U[2] = A[0][2] * Tmp1 + A[1][2] * Tmp2 + A[2][2] * Tmp3;
        }


        template<
        Index Dim
        > inline void TransposeMatMlt(const std::array<std::array<Real, Dim>, Dim> &A,
                               std::array<std::array<Real, Dim>, Dim>& B);


        /*! \brief Spcialization of the MatMlt() function in the 2D case.
        */
        template <
        // Template specialization.
        > inline void TransposeMatMlt<2>(const std::array<std::array<Real, 2>, 2> &A,
                                std::array<std::array<Real, 2>, 2> &B)
        {
            Real b00 = B[0][0];
            Real b01 = B[0][1];
            Real b10 = B[1][0];
            Real b11 = B[1][1];

            B[0][0]=A[0][0]*b00+A[1][0]*b10;
            B[0][1]=A[0][0]*b01+A[1][0]*b11;
            B[1][0]=A[0][1]*b00+A[1][1]*b10;
            B[1][1]=A[0][1]*b01+A[1][1]*b11;
        }
        
        
        /*! \brief Spcialization of the MatMlt() function in the 3D case.
        */
        template <
        // Template specialization.
        > inline void TransposeMatMlt<3>(const std::array<std::array<Real, 3>, 3> &A,
                                std::array<std::array<Real, 3>, 3> &B)
        {
            std::array<std::array<Real, 3>, 3> C = B;

            B[0][0]=A[0][0]*C[0][0]+A[1][0]*C[1][0]+A[2][0]*C[2][0];
            B[0][1]=A[0][0]*C[0][1]+A[1][0]*C[1][1]+A[2][0]*C[2][1];
            B[0][2]=A[0][0]*C[0][2]+A[1][0]*C[1][2]+A[2][0]*C[2][2];

            B[1][0]=A[0][1]*C[0][0]+A[1][1]*C[1][0]+A[2][1]*C[2][0];
            B[1][1]=A[0][1]*C[0][1]+A[1][1]*C[1][1]+A[2][1]*C[2][1];
            B[1][2]=A[0][1]*C[0][2]+A[1][1]*C[1][2]+A[2][1]*C[2][2];

            B[2][0]=A[0][2]*C[0][0]+A[1][2]*C[1][0]+A[2][2]*C[2][0];
            B[2][1]=A[0][2]*C[0][1]+A[1][2]*C[1][1]+A[2][2]*C[2][1];
            B[2][2]=A[0][2]*C[0][2]+A[1][2]*C[1][2]+A[2][2]*C[2][2];
        }

        // ---------------------------------------------------------------------------------//
        
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Performing the transpose matrix times vector operation, results stored in-place.
     */
    template<
    Index Dim
    > inline void TransposeMatMlt(const std::array<std::array<Real, Dim>, Dim> &A,
                           const std::array<Real, Dim>& U,
                                 std::array<Real, Dim>& V);
    
    
    
    /*! \brief Spcialization of the TransposeMatMlt() function in the 2D case.
     */
    template <
    // Template specialization.
    > inline void TransposeMatMlt<2>(const std::array<std::array<Real, 2>, 2> &A,
                              const std::array<Real, 2>& U,
                                    std::array<Real, 2>& V)
    {
        V[0] = A[0][0] * U[0] + A[1][0] * U[1];
        V[1] = A[0][1] * U[0] + A[1][1] * U[1];
    }
    
    /*! \brief Specialization of the TransposeMatMlt() function in the 3D case.
     */
    template <
    // Template specialization.
    > inline void TransposeMatMlt<3>(const std::array<std::array<Real, 3>, 3> &A,
                              const std::array<Real, 3>& U,
                                    std::array<Real, 3>& V)
    {
        V[0] = A[0][0] * U[0] + A[1][0] * U[1] + A[2][0] * U[2];
        V[1] = A[0][1] * U[0] + A[1][1] * U[1] + A[2][1] * U[2];
        V[2] = A[0][2] * U[0] + A[1][2] * U[1] + A[2][2] * U[2];
    }


     template<
 Index DimI,
 Index DimJ,
 Index DimK
 > inline void TransposeMatMlt(const std::array<std::array<Real, DimI>, DimJ> &A,
               const std::array<std::array<Real, DimK>, DimJ> &B,
               std::array<std::array<Real, DimK>, DimI> &AtB)
               {   
                    for (Index i=0;i<DimI;++i)
                        for (Index k=0;k<DimK;++k)
                        {
                            AtB[i][k] = 0.0; 
                            for (Index j=0;j<DimJ;++j) AtB[i][k]+=A[j][i]*B[j][k];
                        }
               }
 


    // ---------------------------------------------------------------------------------//
    
    

        // ---------------------------------------------------------------------------------//
        /*! \brief Performing the matrix times vector operation plus a vector.
         */
        template<
        Index Dim
        > inline void MltAdd(const std::array<std::array<Real, Dim>, Dim> &A,
                      const std::array<Real, Dim> &U,
                      const std::array<Real, Dim> &B,
                      std::array<Real, Dim>& V);
        
        
        
        /*! \brief Spcialization of the MltAdd() function in the 2D case.
         */
        template <
        // Template specialization.
        > inline void MltAdd<2>(const std::array<std::array<Real, 2>, 2> &A,
                         const std::array<Real, 2> &U,
                         const std::array<Real, 2> &B,
                         std::array<Real, 2>& V)
        {
            V[0] = A[0][0] * U[0] + A[0][1] * U[1] + B[0];
            V[1] = A[1][0] * U[0] + A[1][1] * U[1] + B[1];
        }
        // ---------------------------------------------------------------------------------//
        

        // ---------------------------------------------------------------------------------//
        /*! \brief Compute the square of the frobenius norm of a matrix
         */
        template<
        Index DimI,
        Index DimJ
        > inline Real NormFSquared(const std::array<std::array<Real, DimI>, DimJ> &A)
        {
            Real N2 = 0;

            for (Index i=0;i<DimI;++i)
                for (Index j=0;j<DimJ;++j)
                    N2 +=  A[i][j]*A[i][j];

            return N2; 
        }

        // ---------------------------------------------------------------------------------//
        
        


        // ---------------------------------------------------------------------------------//
        /*! \brief Performing addition of matrices with a scaling
         */
        template<
        Index DimI,
        Index DimJ
        > inline void MatScaleAdd(Real alpha,
                                const std::array<std::array<Real, DimJ>, DimI> &A,
                                Real beta,
                                std::array<std::array<Real, DimJ>, DimI> &B)
        {
            for(Index i=0;i<DimI;++i)
                for(Index j=0;j<DimJ;++j)
                {
                    B[i][j] = beta*B[i][j] + alpha*A[i][j];
                }
        }
        // ---------------------------------------------------------------------------------//
        
        
       // ---------------------------------------------------------------------------------//
        /*! \brief Performing computation of the transpose
         */
        template<
        Index DimI,
        Index DimJ
        > inline void Transpose(const std::array<std::array<Real, DimJ>, DimI> &A,
                                std::array<std::array<Real, DimI>, DimJ> &B,
                                Real alpha = 1.0)
        {
            for(Index i=0;i<DimI;++i)
                for(Index j=0;j<DimJ;++j)
                {
                    B[j][i] = alpha*A[i][j];
                }
        }
        // ---------------------------------------------------------------------------------//
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Scaling
         */
        template<
        Index Dim
        > inline void Scale(const Real a,  std::array<Real, Dim>& U);
        
        template <
        // Template specialization.
        > inline void Scale<1> (const Real a,  std::array<Real, 1>& U)
        {
            U[0]*= a;
        }
        
        template <
        // Template specialization.
        > inline void Scale<2> (const Real a,  std::array<Real, 2>& U)
        {
            U[0]*= a;
            U[1]*= a;
        }
        
        template <
        // Template specialization.
        > inline void Scale<3> (const Real a,  std::array<Real, 3>& U)
        {
            U[0]*= a;
            U[1]*= a;
            U[2]*= a;
        }

        /*! \brief Scaling
         */
        template<
        Index DimI,
        Index DimJ
        > inline void Scale(const Real a, const std::array<std::array<Real, DimJ>, DimI>& U
                                 , std::array<std::array<Real, DimJ>, DimI>& V)
                {
                    for (Index i=0;i<DimI;++i)
                        for (Index j=0;j<DimJ;++j)
                            V[i][j] = a*U[i][j];
                }

        /*! \brief Scaling
         */
        template<
        Index DimI,
        Index DimJ
        > inline void Scale(const Real a, std::array<std::array<Real, DimJ>, DimI>& U)
                {
                    for (Index i=0;i<DimI;++i)
                        for (Index j=0;j<DimJ;++j)
                            U[i][j] *= a;
                }

        template<
        Index DimI
        > inline void Scale(const Real a, std::array<Real, DimI>& U)
                {
                    for (Index i=0;i<DimI;++i)
                            U[i] *= a;
                }

        inline void Scale(const Real a, Real & U)
                {
                            U *= a;
                }

        template<
        Index Dim
        > inline void Scale(const Real a, const std::array<Real, Dim>& U
                                 , std::array<Real, Dim>& V)
                {
                    for (Index i=0;i<Dim;++i)
                            V[i] = a*U[i];
                }

        inline void Scale(const Real a, const Real &U, Real &V)
        {
            V = a*U;
        }
        
        
        // ---------------------------------------------------------------------------------//
        





        
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Compute the euclidian norm of a vector
         */
        template<
        Index Dim
        > inline Real Norm(const std::array<Real, Dim> &);
        
        /*! \brief Spcialization of the Norm() function in the 1D case.
         */
        template <
        // Template specialization.
        > inline Real Norm<1>(const std::array<Real, 1> & U)
        {
            return fabs(U[0]);
        }
        
        /*!\brief Spcialization of the Norm() function in the 2D case.
         */
        template <
        // Template specialization.
        > inline Real Norm<2>(const std::array<Real, 2> & U)
        {
            return sqrt(U[0]*U[0]+U[1]*U[1]);
        }
        
        /*!\brief Spcialization of the Norm() function in the 3D case.
         */
        template <
        // Template specialization.
        > inline Real Norm<3>(const std::array<Real, 3> & U)
        {
            return sqrt(U[0]*U[0]+U[1]*U[1]+U[2]*U[2]);
        }
        // ---------------------------------------------------------------------------------//
        
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Compute the euclidian distance of a vector
         */
        template<
        Index Dim
        > inline Real Dist(const std::array<Real, Dim> &, const std::array<Real, Dim> &);
        
        /*! \brief Spcialization of the Dist() function in the 1D case.
         */
        template <
        // Template specialization.
        > inline Real Dist<1>(const std::array<Real, 1> & U, const std::array<Real, 1> & V)
        {
            return fabs(U[0]-V[0]);
        }
        
        /*!\brief Spcialization of the Dist() function in the 2D case.
         */
        template <
        // Template specialization.
        > inline Real Dist<2>(const std::array<Real, 2> & U, const std::array<Real, 2> & V)
        {
            Real Distance = 0;
            for (Index dim=0;dim<2;dim++)
                Distance += (U[dim] - V[dim])*(U[dim] - V[dim]);
            return sqrt(Distance);

        }
        
        /*!\brief Spcialization of the Dist() function in the 3D case.
         */
        template <
        // Template specialization.
        > inline Real Dist<3>(const std::array<Real, 3> & U, const std::array<Real, 3> & V)
        {
            Real Distance = 0;
            for (Index dim=0;dim<3;dim++)
                Distance += (U[dim] - V[dim])*(U[dim] - V[dim]);
            return sqrt(Distance);
        }
        // ---------------------------------------------------------------------------------//
        
        
        
        // ---------------------------------------------------------------------------------//
        /*! \brief Compute the cross product between two vectors
         */
        
        template<
        Index Dim
        > inline void CrossProduct(const std::array<Real, Dim> & U, const std::array<Real, Dim> & V, std::array<Real, Dim> & W)
        {
            
        }
        
        template <
        // Template specialization.
        >
        inline void CrossProduct<3>(const std::array<Real, 3> & U, const std::array<Real, 3> & V, std::array<Real, 3> & W)
        {
            W[0]=U[1]*V[2]-U[2]*V[1];
            W[1]=U[2]*V[0]-U[0]*V[2];
            W[2]=U[0]*V[1]-U[1]*V[0];
        }
    
        template <
        Index Dim
        > inline Real ScalarProduct(const std::array<Real, Dim> & U, const std::array<Real, Dim> & V)
        {
            Real sp = 0;
            for(Index I=0; I < Dim; ++I) sp += U[I]*V[I];
            return sp;
        }

        template <
        Index Dim
        > inline Real ContractionProduct(const std::array<std::array<Real, Dim>,Dim> & U, const std::array<std::array<Real, Dim>,Dim>  & V)
        {
            Real cp = 0;
            for(Index I=0; I < Dim; ++I)
                for(Index J=0; J < Dim; ++J)
                    cp += U[I][J]*V[I][J];
            return cp;
        }

        
        // ---------------------------------------------------------------------------------//


        template<
        Index Dim
        > inline void VecMlt(const std::array<Real, Dim> &A,
                    const std::array<Real, Dim> &B,
                    Real &AB)
                    {   
                        AB = ScalarProduct(A,B);     
                    }
 

         template<
        Index Dim
        > inline void TransposeVecMlt(const std::array<Real, Dim> &A,
                                const double &B,
                                std::array<Real, Dim> &AB)
                    {   
                         Scale(B,A,AB);   
                    }
 


        
    } //ArrayAlgebra

} // OndoMathX

// -----------------------------------------------------------------------------------------//










//
//
//        typedef std::array< std::array<Real, 3>, 3 > inline RealMatrix3x3;
//        typedef std::array< std::array<Real, 2>, 2 > inline RealMatrix2x2;
//        typedef std::array< std::array<Real, 3>, 2 > Matrix2x3;
//        typedef std::array< std::array<Real, 2>, 3 > Matrix3x2;
//
//        typedef std::array< Real, 3 > Vector3;
//        typedef std::array< Real, 2 > Vector2;
//
//
//        /*interpole en 1D la valeur en `p` de la fonction numéro `phi`. Les fonctions
//         * formant une base lagrangienne qui s'appuie sur les `n` points `pos`.
//         * `phi` doit être compris entre 0 et n-1. On ne doit pas avoir deux positions identiques.
//         */
//        inline Real lgrInterp(const std::vector<Real> &pos,Real p,Index phi)
//        {
//            Real s=1.0;
//            for (Index r=0;r<pos.size();++r) if (r!=phi) s*=(p-pos[r])/(pos[phi]-pos[r]);
//            return s;
//        }
//
//
//        /*interpole en 1D la valeur en `p` de la dérivée de la fonction numéro `phi`. Les fonctions
//         * formant une base lagrangienne qui s'appuie sur les `n` points `pos`.
//         * `phi` doit être compris entre 0 et n-1. On ne doit pas avoir deux positions identiques.
//         */
//        inline Real dLgrInterp(const std::vector<Real> &pos,Real p,Index phi)
//        {
//            Real c;
//            Real s = 0.0;
//            Real d = 1.0;
//
//            for (Index r=0;r<pos.size();++r) if (r!=phi) d*=(pos[phi]-pos[r]);
//
//            for (Index r1=0;r1<pos.size();++r1)
//            {
//                if (r1!=phi)
//                {
//                    c = 1.0;
//                    for (Index r2=0;r2<pos.size();++r2) if (r2!=r1 && r2!=phi) c*=(p-pos[r2]);
//                    s += c;
//                }
//            }
//
//            return s/d;
//        }
//
//        //Numerical tools for 2x2 matrices////////////////////////////////////////////
//        inline void prod2x2(const RealMatrix2x2 & A, const RealMatrix2x2 & B, RealMatrix2x2 & AB)
//        {
//            AB[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0];
//            AB[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1];
//            AB[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0];
//            AB[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1];
//        }
//
//        inline void prod2x2(const RealMatrix2x2 & A, const RealMatrix2x2 & B, const RealMatrix2x2 & C, RealMatrix2x2 & ABC)
//        {
//            RealMatrix2x2 tmp;
//            prod2x2(B,C,tmp);
//            prod2x2(A,tmp,ABC);
//        }
//
//
//        inline void prod2x2(const RealMatrix2x2 & A,const Vector2 &x, Vector2 &Ax)
//        {
//            Ax[0]=A[0][0]*x[0]+A[0][1]*x[1];
//            Ax[1]=A[1][0]*x[0]+A[1][1]*x[1];
//        }
//
//        inline void transp2x2(const RealMatrix2x2 & mat, RealMatrix2x2 & tr)
//        {
//            tr[0][0] =  mat[0][0];
//            tr[1][0] =  mat[0][1];
//            tr[0][1] =  mat[1][0];
//            tr[1][1] =  mat[1][1];
//        }
//
//        inline void scale2x2(Real alpha, RealMatrix2x2 & mat)
//        {
//            mat[0][0] = alpha*mat[0][0];
//            mat[0][1] = alpha*mat[0][1];
//            mat[1][0] = alpha*mat[1][0];
//            mat[1][1] = alpha*mat[1][1];
//        }
//
//        inline Real det2x2(const RealMatrix2x2 & mat)
//        {
//            return mat[0][0] * mat[1][1] - mat[1][0] * mat[0][1];
//        }
//
//        inline void com2x2(const RealMatrix2x2 & mat, RealMatrix2x2 & com)
//        {
//            com[0][0] =  mat[1][1];
//            com[1][0] = -mat[0][1];
//            com[0][1] = -mat[1][0];
//            com[1][1] =  mat[0][0];
//        }
//
//        inline void tcom2x2(const RealMatrix2x2 & mat, RealMatrix2x2 & tcom)
//        {
//            tcom[0][0] =  mat[1][1];
//            tcom[1][0] = -mat[1][0];
//            tcom[0][1] = -mat[0][1];
//            tcom[1][1] =  mat[0][0];
//        }
//
//        inline void inv2x2(const RealMatrix2x2 & mat, RealMatrix2x2 & inv)
//        {
//            Real det = det2x2(mat);
//            tcom2x2(mat,inv);
//            scale2x2(1.0/det,inv);
//        }
//
//        inline Real norm2(const Vector2 &x)
//        {
//            return sqrt(x[0]*x[0]+x[1]*x[1]);
//        }
//
//        inline Real prodScal2(const Vector2 &x, const Vector2 &y)
//        {
//            return x[0]*y[0] + x[1]*y[1];
//        }
//
//        inline Real prodScal2(const RealMatrix2x2 & A, const Vector2 &x, const Vector2 &y)
//        {
//            return  A[0][0]*x[0]*y[0] + A[0][1]*x[1]*y[0] + A[1][0]*x[0]*y[1] + A[1][1]*x[1]*y[1];
//        }
//
//        inline void prodve(const Vector3 &a, const Vector3 &b, Vector3 &c)
//        {
//            c[2] = a[0] * b[1] - a[1] * b[0];
//            c[1] = -a[0] * b[2] + a[2] * b[0];
//            c[0] = a[1] * b[2] - a[2] * b[1];
//        }
//
//        //Numerical tools for 3x3 matrices and vectors////////////////////////////////////////////////////
//
//        inline void prod3x3(const RealMatrix3x3 & A,const RealMatrix3x3 & B,RealMatrix3x3 & AB)
//        {
//            AB[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0];
//            AB[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1];
//            AB[0][2]=A[0][0]*B[0][2]+A[0][1]*B[1][2]+A[0][2]*B[2][2];
//
//            AB[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0];
//            AB[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1];
//            AB[1][2]=A[1][0]*B[0][2]+A[1][1]*B[1][2]+A[1][2]*B[2][2];
//
//            AB[2][0]=A[2][0]*B[0][0]+A[2][1]*B[1][0]+A[2][2]*B[2][0];
//            AB[2][1]=A[2][0]*B[0][1]+A[2][1]*B[1][1]+A[2][2]*B[2][1];
//            AB[2][2]=A[2][0]*B[0][2]+A[2][1]*B[1][2]+A[2][2]*B[2][2];
//        }
//
//        inline void prod3x3(const RealMatrix3x3 & A,const Vector3 &x, Vector3 &Ax)
//        {
//            Ax[0]=A[0][0]*x[0]+A[0][1]*x[1]+A[0][2]*x[2];
//            Ax[1]=A[1][0]*x[0]+A[1][1]*x[1]+A[1][2]*x[2];
//            Ax[2]=A[2][0]*x[0]+A[2][1]*x[1]+A[2][2]*x[2];
//        }
//
//        inline Real norm3(const Vector3 &x)
//        {
//            return sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
//        }
//
//        inline Real dist3squared(const Vector3 &x,const Vector3 &y)
//        {
//            return ((x[0] - y[0])*(x[0] - y[0])+(x[1] - y[1])*(x[1] - y[1])+(x[2]-y[2])*(x[2]-y[2]));
//        }
//
//        inline Real dist3(const Vector3 &x,const Vector3 &y)
//        {
//            return sqrt(dist3squared(x,y));
//        }
//
//        inline Real prodScal3(const Vector3 &x,const Vector3 &y)
//        {
//            return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
//        }
//
//        inline void prod3x2(const Matrix3x2 & A,const Vector2 &x, Vector3 &Ax)
//        {
//            Ax[0]=A[0][0]*x[0]+A[0][1]*x[1];
//            Ax[1]=A[1][0]*x[0]+A[1][1]*x[1];
//            Ax[2]=A[2][0]*x[0]+A[2][1]*x[1];
//        }
//
//        inline void prod2x3(const Matrix2x3 & A,const Vector3 &x, Vector3 &Ax)
//        {
//            Ax[0]=A[0][0]*x[0]+A[0][1]*x[1]+A[0][2]*x[2];
//            Ax[1]=A[1][0]*x[0]+A[1][1]*x[1]+A[1][2]*x[2];
//        }
//
//        inline void prod2x3(const Matrix2x3 & A,const Matrix3x2 & B, RealMatrix2x2 & AB)
//        {
//            AB[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0]+A[0][2]*B[2][0];
//            AB[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1]+A[0][2]*B[2][1];
//            AB[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0]+A[1][2]*B[2][0];
//            AB[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1]+A[1][2]*B[2][1];
//        }
//
//        inline void transp3x2(const Matrix3x2 & mat, Matrix2x3 & tr)
//        {
//            tr[0][0] =  mat[0][0];
//            tr[0][1] =  mat[1][0];
//            tr[0][2] =  mat[2][0];
//
//            tr[1][0] =  mat[0][1];
//            tr[1][1] =  mat[1][1];
//            tr[1][2] =  mat[2][1];
//        }
//
//        inline void prod3x3(const Matrix3x2 & A,const Matrix2x3 & B,RealMatrix3x3 & AB)
//        {
//            AB[0][0]=A[0][0]*B[0][0]+A[0][1]*B[1][0];
//            AB[0][1]=A[0][0]*B[0][1]+A[0][1]*B[1][1];
//            AB[0][2]=A[0][0]*B[0][2]+A[0][1]*B[1][2];
//
//            AB[1][0]=A[1][0]*B[0][0]+A[1][1]*B[1][0];
//            AB[1][1]=A[1][0]*B[0][1]+A[1][1]*B[1][1];
//            AB[1][2]=A[1][0]*B[0][2]+A[1][1]*B[1][2];
//
//            AB[2][0]=A[2][0]*B[0][0]+A[2][1]*B[1][0];
//            AB[2][1]=A[2][0]*B[0][1]+A[2][1]*B[1][1];
//            AB[2][2]=A[2][0]*B[0][2]+A[2][1]*B[1][2];
//        }
//
//        inline void prod3x3(const RealMatrix3x3 & A,const RealMatrix3x3 & B,const RealMatrix3x3 & C,RealMatrix3x3 & ABC)
//        {
//            RealMatrix3x3 BC;
//            prod3x3(B,C,BC);
//            prod3x3(A,BC,ABC);
//        }
//
//        inline void transp3x3(const RealMatrix3x3 & mat,RealMatrix3x3 & tr)
//        {
//            tr[0][0] =  mat[0][0];
//            tr[0][1] =  mat[1][0];
//            tr[0][2] =  mat[2][0];
//
//            tr[1][0] =  mat[0][1];
//            tr[1][1] =  mat[1][1];
//            tr[1][2] =  mat[2][1];
//
//            tr[2][0] =  mat[0][2];
//            tr[2][1] =  mat[1][2];
//            tr[2][2] =  mat[2][2];
//        }
//
//
//        inline void scale3x3(Real alpha, RealMatrix3x3 & mat)
//        {
//            mat[0][0] = alpha*mat[0][0];
//            mat[0][1] = alpha*mat[0][1];
//            mat[0][2] = alpha*mat[0][2];
//
//            mat[1][0] = alpha*mat[1][0];
//            mat[1][1] = alpha*mat[1][1];
//            mat[1][2] = alpha*mat[1][2];
//
//            mat[2][0] = alpha*mat[2][0];
//            mat[2][1] = alpha*mat[2][1];
//            mat[2][2] = alpha*mat[2][2];
//        }
//
//        inline Real det3x3(const RealMatrix3x3 & mat)
//        {
//            return (mat[0][0] * (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) -
//                    mat[0][1] * (mat[1][0] * mat[2][2] - mat[1][2] * mat[2][0]) +
//                    mat[0][2] * (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]));
//        }
//
//        inline void com3x3(const RealMatrix3x3 & mat, RealMatrix3x3 & com)
//        {
//            com[0][0] =  mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2];
//            com[0][1] =  -mat[1][0]*mat[2][2]+mat[2][0]*mat[1][2];
//            com[0][2] =  mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1];
//
//            com[1][0] =  -mat[0][1]*mat[2][2]+mat[2][1]*mat[0][2];
//            com[1][1] =  mat[0][0]*mat[2][2]-mat[2][0]*mat[0][2];
//            com[1][2] =  -mat[0][0]*mat[2][1]+mat[2][0]*mat[0][1];
//
//            com[2][0] =  mat[0][1]*mat[1][2]-mat[1][1]*mat[0][2];
//            com[2][1] =  -mat[0][0]*mat[1][2]+mat[1][0]*mat[0][2];
//            com[2][2] =  mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1];
//        }
//
//        inline void tcom3x3(const RealMatrix3x3 & mat, RealMatrix3x3 & tcom)
//        {
//            tcom[0][0] =  mat[1][1]*mat[2][2]-mat[2][1]*mat[1][2];
//            tcom[1][0] =  -mat[1][0]*mat[2][2]+mat[2][0]*mat[1][2];
//            tcom[2][0] =  mat[1][0]*mat[2][1]-mat[2][0]*mat[1][1];
//
//            tcom[0][1] =  -mat[0][1]*mat[2][2]+mat[2][1]*mat[0][2];
//            tcom[1][1] =  mat[0][0]*mat[2][2]-mat[2][0]*mat[0][2];
//            tcom[2][1] =  -mat[0][0]*mat[2][1]+mat[2][0]*mat[0][1];
//
//            tcom[0][2] =  mat[0][1]*mat[1][2]-mat[1][1]*mat[0][2];
//            tcom[1][2] =  -mat[0][0]*mat[1][2]+mat[1][0]*mat[0][2];
//            tcom[2][2] =  mat[0][0]*mat[1][1]-mat[1][0]*mat[0][1];
//        }
//
//        inline void inv3x3(const RealMatrix3x3 & mat, RealMatrix3x3 & inv)
//        {
//            Real det = det3x3(mat);
//            tcom3x3(mat,inv);
//            scale3x3(1.0/det,inv);
//        }
//
//        inline void prodvect3x3(const Vector3 u,const Vector3 v,Vector3 uXv)
//        {
//            uXv[0]=u[1]*v[2]-u[2]*v[1];
//            uXv[1]=u[2]*v[0]-u[0]*v[2];
//            uXv[2]=u[0]*v[1]-u[1]*v[0];
//        }
//
 
