#pragma once

// -----------------------------------------------------------------------------------------//
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <memory>
#include <map>
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
namespace OndoMathX {

namespace LAL {

// ----------------------------------------------------------------------------------//
/*! \brief Definition of the object for matrix and vector subscripting and data types.
 */

/*! \brief Definition of the object storing the vector of degree of freedom.
 */
typedef ::Eigen::Matrix<Real, ::Eigen::Dynamic, 1> Vector;

/*! \brief Definition of the object storing the vector of degree of freedom.
 */
typedef ::Eigen::Matrix<Complex, ::Eigen::Dynamic, 1> ComplexVector;


/*! \brief Definition of sparse matrix.
 */
typedef ::Eigen::SparseMatrix<Real, ::Eigen::RowMajor, int> SparseMatrix;

/*! \brief Definition of diagonal matrix.
 */
typedef ::Eigen::DiagonalMatrix<Real, ::Eigen::Dynamic> DiagonalMatrix;

/*! \brief Definition of real dense matrix.
 */
typedef ::Eigen::Matrix<Real, ::Eigen::Dynamic, ::Eigen::Dynamic> DenseMatrix;

/*! \brief Definition of complex dense matrix.
 */
typedef ::Eigen::Matrix<Complex, ::Eigen::Dynamic, ::Eigen::Dynamic> ComplexDenseMatrix;

/*! \brief Definition of Factorization matrix.
 */
typedef ::Eigen::SparseLU<::Eigen::SparseMatrix<Real, ::Eigen::ColMajor, int>, ::Eigen::COLAMDOrdering<int>> FactorizationMatrix;
// ----------------------------------------------------------------------------------//



// ---------------------------------------------------------------------------------//
/*! \brief Allocate a Vector.
 \param aVector is a Vector to allocate.
 \param Dim is the dimension of the Vector after successful allocation.
 */
inline void Allocate(Vector & aVector, Index Dim)
{
    // Resize vector and fill it with 0.
    aVector.setZero(Dim);
}


/*! \brief Extracting the dimension of a Vector.
 \param aVector is the vector of interest.
 \return The dimension of the Vector as an Index object.
 */
inline Index getDimension(const Vector& aVector)
{
    return  (Index)aVector.size();
}


/*! \brief Extracting pointer on dataVector.
 \param aVector is the vector of interest.
 \return The pointer on data
 */
inline Real * getData(Vector& aVector)
{
    return aVector.data();
}

inline const Real * getData(const Vector& aVector)
{
    return aVector.data();
}


//  ---------------------------------------------------------------------------------//



// ---------------------------------------------------------------------------------//
/*! \brief Allocate a ComplexVector.
 \param aVector is a ComplexVector to allocate.
 \param Dim is the dimension of the Vector after successful allocation.
 */
inline void Allocate(ComplexVector & aVector, Index Dim)
{
    // Resize vector and fill it with 0.
    aVector.setZero(Dim);
}


/*! \brief Extracting the dimension of a Vector.
 \param aVector is the vector of interest.
 \return The dimension of the Vector as an Index object.
 */
inline Index getDimension(const ComplexVector& aVector)
{
    return  (Index)aVector.size();
}


/*! \brief Extracting pointer on dataVector.
 \param aVector is the vector of interest.
 \return The pointer on data
 */
inline Complex * getData(ComplexVector& aVector)
{
    return aVector.data();
}

inline const Complex * getData(const ComplexVector& aVector)
{
    return aVector.data();
}


//  ---------------------------------------------------------------------------------//




//  ---------------------------------------------------------------------------------//
/*! \brief Allocate a sparse matrix.
 \param aSparseMatrix is a sparse matrix of interest.
 \param NRow is the number of rows to allocate.
 \param NCol is the number of columns to allocate.
 */
inline void Allocate(SparseMatrix & aSparseMatrix,
              Index NRow,Index NCol)
{
    // Resizing.
    aSparseMatrix.resize(NRow, NCol);
}



// Create the block matrix
//  [ A B ]
//  [ C D ]
inline void CreateBlockMatrix(const SparseMatrix &A,
                       const SparseMatrix &B,
                       const SparseMatrix &C,
                       const SparseMatrix &D,
                       SparseMatrix &ABCD)
{
    Index MA = A.cols();
    Index MB = B.cols();
    Index NA = A.rows();
    Index NC = C.rows();
    
    LAL::Allocate(ABCD,MA+MB, NA+NC);
            
    std::vector<::Eigen::Triplet<double> > tripletList;
    
    for (int k = 0; k < A.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(A, k); it; ++it)
            tripletList.push_back(::Eigen::Triplet<double>(it.row(),it.col(),it.value()));
    
    for (int k = 0; k < B.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(B, k); it; ++it)
        {
            tripletList.push_back(::Eigen::Triplet<double>(it.row(),it.col()+MA,it.value()));
        }
    
    for (int k = 0; k < C.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(C, k); it; ++it)
        {
            tripletList.push_back(::Eigen::Triplet<double>(it.row()+NA,it.col(),it.value()));
        }
    
    for (int k = 0; k < D.outerSize(); ++k)
        for (SparseMatrix::InnerIterator it(D, k); it; ++it)
        {
            tripletList.push_back(::Eigen::Triplet<double>(it.row()+NA,it.col()+MA,it.value()));
        }
    
    ABCD.setFromTriplets(tripletList.begin(), tripletList.end());
}

inline void EraseLines(SparseMatrix &A,const std::vector<Index> &Lines)
{
    std::map<Index,bool> Lines_map;

    for (Index k = 0; k < Lines.size(); ++k) Lines_map[Lines[k]] = true;

    for (Index k = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it)
        {
            if (Lines_map.count(it.row()) == 1) it.valueRef() = 0;
        }
    }
}

inline void EraseColumns(SparseMatrix &A,std::vector<Index> Columns)
{
    std::map<Index,bool> Columns_map;

    for (Index k = 0; k < Columns.size(); ++k) Columns_map[Columns[k]] = true;

    for (Index k = 0; k < A.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(A, k); it; ++it)
        {
            if (Columns_map.count(it.row()) == 1) it.valueRef() = 0;
        }
    }
}



 

/*! \brief Compute a factorization of a matrix for inversion.
 \param anInvertibleMatrix will bear the result of the factorization
 \param aSparseMatrix is the matrix to invert.
 */
inline void Factorize(FactorizationMatrix & aFactorization,
               SparseMatrix& aSparseMatrix)
{

    //Make a copy in ColMajor format to optimize factorization
    ::Eigen::SparseMatrix<Real, ::Eigen::ColMajor, int> aSparseMatrixCopy(aSparseMatrix);

    //Compress the matrix
    aSparseMatrixCopy.makeCompressed();

    // Compute the ordering permutation vector from the structural pattern of the matrix.
    aFactorization.analyzePattern(aSparseMatrixCopy);

    // Creating factorization.
    aFactorization.factorize(aSparseMatrixCopy);
}



/*! \brief Adding interaction in a sparse matrix, using a type-relevant reference for efficiency.
 \param aSparseMatrix is the sparse matrix of interest.
 \param IdxRow is the index used as the matrix subscript in row.
 \param IdxCol is the index used as the matrix subscript in column.
 \param Value is the value of the interaction to add.
 */
inline void AddInteraction(SparseMatrix& aSparseMatrix,
                    Index IdxRow, Index IdxCol, Real Value)
{
    aSparseMatrix.coeffRef(IdxRow, IdxCol) += Value;
}

inline void AddInteractions(SparseMatrix &A,const std::vector<Index> & Lines,const std::vector<Index> & Columns, Real value)
{
    assert(Lines.size() == Columns.size());

    for (Index k = 0; k < Lines.size(); ++k) 
    {
        AddInteraction(A,Lines[k],Columns[k],value);
    }
 
}




/*! \brief Accessing non-zero value in a column of a sparse matrix.
 \param aSparseMatrix is the sparse matrix of interest.
 \param RowIdx index of the row.
 \return A pointer on the array containing the data.
 */
//Real* GetRowData(SparseMatrix& aSparseMatrix, Index RowIdx)
//{
//    return aSparseMatrix.valuePtr() + aSparseMatrix.outerIndexPtr()[RowIdx];
//}



//            /*! \brief Accessing row data indices.
//             \param aSparseMatrix is the sparse matrix of interest.
//             \param RowIdx index of the row.
//             \return A pointer on the array containing the row indices
//             */
//            Index* GetColumnInd(SparseMatrix& aSparseMatrix, Index RowIdx)
//            {
//                return aSparseMatrix.innerIndexPtr() + aSparseMatrix.outerIndexPtr()[RowIdx];
//            }



/*! \brief Accessing number of non-zero value in a column.
 \param aSparseMatrix is the sparse matrix of interest.
 \param RowIdx index of the row.
 \return Number of non-zero value in a column.
 */
//Index GetNRowData(SparseMatrix& aSparseMatrix, Index RowIdx)
//{
//    if (aSparseMatrix.isCompressed())
//    {
//        return aSparseMatrix.outerIndexPtr()[RowIdx + 1] - aSparseMatrix.outerIndexPtr()[RowIdx];
//    }
//    else
//    {
//        return aSparseMatrix.innerNonZeroPtr()[RowIdx];
//    }
//}
//

/*! \brief Compress the matrix for faster MltAdd but slower modification
 \param aSparseMatrix is the sparse matrix to compress.
 */
inline void Compress(SparseMatrix& aSparseMatrix)
{
    aSparseMatrix.makeCompressed();
}

 


 
 

/*! \brief Performing
 
 \f$ \mathcal{O}^{-1} U \longrightarrow  U \f$,
 
 where \f$\mathcal{O}\f$ is a factorized matrix.
 \param aFactorizedMatrix is the factorization of the matrix appearing in the linear system.
 \param U is a Vector bearing the result and the rhs.
 */
inline void Solve(FactorizationMatrix & aFactorization,Vector& U)
{
    U = aFactorization.solve(U);
}

inline void Allocate(DiagonalMatrix & aDiagMatrix, Index nRow, Index nColumn)
{
    assert(nRow == nColumn);
    aDiagMatrix.resize(nRow);
    aDiagMatrix.setZero();
}

inline void AddInteraction(DiagonalMatrix & aDiagMatrix,Index iRow, Index iColumn, Real Value)
{
    assert(iRow == iColumn);
    getData(aDiagMatrix.diagonal())[iRow] += Value;
}

inline Real getValue(DiagonalMatrix & aDiagMatrix,Index iRow, Index iColumn)
{
    assert(iRow == iColumn);
    return getData(aDiagMatrix.diagonal())[iRow];
}

inline void setValue(DiagonalMatrix & aDiagMatrix,Index iRow, Index iColumn, Real Value)
{
    assert(iRow == iColumn);
    getData(aDiagMatrix.diagonal())[iRow] = Value;
}

 
 

// ---------------------------------------------------------------------------------//
/*! \brief Perform
 \f$ \mathcal{O} ^ {-1} U \longrightarrow U\f$.
 \param V is a vector of double representing the right hand side of the linear system and
 bearing the result.
 */
inline void Solve( DiagonalMatrix & aDiagMatrix, Vector & refU)
{
    Real * D = LAL::getData(aDiagMatrix.diagonal());
    Real * U = LAL::getData(refU);
    Index Dim = getDimension(aDiagMatrix.diagonal());
    
#pragma omp parallel for
    for (Index iDoF = 0; iDoF < Dim; iDoF++)
    {
        
        U[iDoF] /= D[iDoF];
    }
}
// ---------------------------------------------------------------------------------//

inline void Invert( DiagonalMatrix & aDiagMatrix)
{
    Real * D = LAL::getData(aDiagMatrix.diagonal());
    Index Dim = getDimension(aDiagMatrix.diagonal());
    
    for (Index iDoF = 0; iDoF < Dim; iDoF++)
    {
        D[iDoF] = 1.0/D[iDoF];
    }
    
}


/*! \brief Accessing number of rows.
 \return The number of rows in the matrix.
 */
inline Index getNumRows(DiagonalMatrix & aDiagMatrix)
{
    return getDimension(aDiagMatrix.diagonal());
}

/*! \brief Accessing number of columns.
 \return The number of columns in the matrix.
 */
inline Index getNumColumns(DiagonalMatrix & aDiagMatrix)
{
    return getNumRows(aDiagMatrix);
}
// ---------------------------------------------------------------------------------//


inline void Scale(const Real & aScaling, Vector & aVector)
{
    Real * data = getData(aVector);
    Index NDoF = getDimension(aVector);
    
    for (Index iDoF = 0.0; iDoF < NDoF; iDoF++)
        data[iDoF] *= aScaling;
}


/*! \brief  Get the infinite norm
 \param aVector is a Vector to assign to bVector
 \return the infinite norm of the vector
 */
inline Real GetInfiniteNorm(
                                   Vector& aVector
                                   )
{
    //Recovers the data
    Real * data = getData(aVector);
    Index NDoF = getDimension(aVector);

    Real Norm = fabs(data[0]);
    

    for (Index iDoF = 1; iDoF < NDoF; iDoF++)
        if (fabs(data[iDoF]) > Norm) Norm = fabs(data[iDoF]);
    
    return Norm;
}


// ---------------------------------------------------------------------------------//
/*! \brief Scaling a complex vector.
 \param aScale is the coefficient used for scaling.
 \param aVector is a Vector to scale.
 */


inline void Copy(const Vector & aVector1, Vector & aVector2)
{
    const Real * data1 = getData(aVector1);
    Real * data2 = getData(aVector2);
    
    Index NDoF1 = getDimension(aVector1);
    Index NDoF2 = getDimension(aVector2);
    
    assert(NDoF1==NDoF2);
    
    for (Index iDoF = 0.0; iDoF < NDoF1; iDoF++)
        data2[iDoF] = data1[iDoF];
}





/*! \brief Get the euclidian norm
 \param aVector is the Vector which norm is computed
 \return the euclidian norm of the vector
 */
inline static Real EuclidianNorm(const Vector& aVector)
{
    //Recovers the data
    const Real * data = getData(aVector);
    Index NDoF = getDimension(aVector);
    
    Real Norm = 0.0;
    

    for (Index iDoF = 0; iDoF < NDoF; iDoF++)
        Norm += data[iDoF] * data[iDoF];
    
    return sqrt(Norm);
}


/*! \brief Get the scalar product
 \param aVector1 is the first Vector to use for the scalar product
 \param aVector2 is the second Vector to use for the scalar product
 \return the  scalar product of the two vectors
 */
inline Real ScalarProduct(const Vector& aVector1, const Vector& aVector2)
{
    //Recovers the data
    const Real * data1 = getData(aVector1);
    const Real * data2 = getData(aVector2);
    Index NDoF = getDimension(aVector1);
    
    Real sc = 0.0;
    
    for (Index iDoF = 0; iDoF < NDoF; iDoF++)
        sc += data1[iDoF] * data2[iDoF];
    
    return sc;
}



inline Real Sum(const Vector& aVector1)
{
    //Recovers the data
    const Real * data = getData(aVector1);
    Index NDoF = getDimension(aVector1);
    
    Real s = 0.0;
    
    for (Index iDoF = 0; iDoF < NDoF; iDoF++)
        s += data[iDoF];
    
    return s;
}




// Methods for real  dense matrices ----------------------------------------------------------------------------------//


inline void Allocate(DenseMatrix & aDenseMatrix, Index nRow, Index nColumn)
{
    aDenseMatrix.resize(nRow,nColumn);
    aDenseMatrix.setZero();
}


inline Real getValue(DenseMatrix & aDenseMatrix,Index iRow, Index iColumn)
{
    return aDenseMatrix(iRow,iColumn);
}

inline void setValue(DenseMatrix & aDenseMatrix,Index iRow, Index iColumn, Real Value)
{
    aDenseMatrix(iRow,iColumn) = Value;
}


// Methods for complex dense matrices ----------------------------------------------------------------------------------//


inline void Allocate(ComplexDenseMatrix & aComplexDenseMatrix, Index nRow, Index nColumn)
{
    aComplexDenseMatrix.resize(nRow,nColumn);
    aComplexDenseMatrix.setZero();
}

inline Complex getValue(ComplexDenseMatrix & aComplexDenseMatrix,Index iRow, Index iColumn)
{
    return aComplexDenseMatrix(iRow,iColumn);
}

inline void setValue(ComplexDenseMatrix & aComplexDenseMatrix,Index iRow, Index iColumn, Real Value)
{
    aComplexDenseMatrix(iRow,iColumn) = Value;
}

//void GeneralizedSelfAdjointEigenvalueDecomposition(ComplexDenseMatrix &A, aComplexDenseMatrix &B,std:vector<ComplexVector> &eigenvectors, //std:vector<Real> &eigenvalues)
//{
//    Index N = A.rows();
//
//    GeneralizedSelfAdjointEigenSolver<ComplexDenseMatrix> es(A,B);
//
//    eigenvalues.resize(N);
//    eigenvectors.resize(N);
//
//    for (Index n=0;n<N;n++)
//    {
//        eigenvalues[n]=es.eigenvalues()[n];
//        eigenvectors[n]=es.eigenvectors().col(n);
//    }
//}
// ----------------------------------------------------------------------------------//

} // LAL

} // OndoMathX


