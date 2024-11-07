#pragma once



// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{


template<Index SysDim = 1, class FESpace, class ParameterField, class SparseMatrix> void AssembleMassBndy(FESpace & aFESpace,  ParameterField & aParameterField, SparseMatrix & aMatrix, Real Alpha)
{
    Core::Assemble<FESpace,FESpace,true,true,true,
                   Core::QuadratureOpt::Default,
                   SysDim,  //DimU
                   SysDim,  //DimV
                   1,       //DimI
                   SysDim,  //DimJ
                   ParameterField,
                   SparseMatrix,
                   Core::DiffOp::Trace,
                   Core::DiffOp::Trace>(aFESpace,aFESpace,aParameterField,aMatrix,Alpha);

}

template<Index SysDim = 1, class FESpace, class ParameterField, class SparseMatrix> void CreateMassBndy(FESpace & aFESpace,  ParameterField & aParameterField, SparseMatrix & aMatrix, Real Alpha)
{
    LAL::Allocate(aMatrix, SysDim*aFESpace.getNumDoFs(), SysDim*aFESpace.getNumDoFs());

    AssembleMassBndy<SysDim>(aFESpace,aParameterField,aMatrix,Alpha);
}




 
 
 

 
} // OndoMathX
// -----------------------------------------------------------------------------------------//




