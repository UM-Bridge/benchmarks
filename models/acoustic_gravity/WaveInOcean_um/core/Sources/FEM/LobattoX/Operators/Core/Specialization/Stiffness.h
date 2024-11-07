#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 
template<Index SysDim = 1, class FESpace, class ParameterField, class SparseMatrix> void CreateStiffness(FESpace & aFESpace,  ParameterField & aParameterField, SparseMatrix & aMatrix, Real Alpha)
{
    LAL::Allocate(aMatrix, SysDim*aFESpace.getNumDoFs(), SysDim*aFESpace.getNumDoFs());

/*
    Core::Assemble<FESpace,true,true,
                   SysDim,       //DimU
                   SysDim,       //DimV
                   SysDim,       //DimI
                   FESpace::Dim, //DimJ
                   ParameterField,
                   SparseMatrix,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Gradient>(aFESpace,aParameterField,aMatrix,Alpha); */

    Core::Assemble<FESpace,FESpace,true,true,true,
                Core::QuadratureOpt::Default,
                SysDim,       //DimU
                SysDim,       //DimV
                SysDim,       //DimI
                FESpace::Dim, //DimJ
                ParameterField,
                SparseMatrix,
                Core::DiffOp::Gradient,
                Core::DiffOp::Gradient>(aFESpace,aFESpace,aParameterField,aMatrix,Alpha);

}


template<Index SysDim = 1, class FESpace, class ParameterField>
void MltAddStiffness(FESpace & aFESpace,
                     ParameterField & aParameterField,
                     Real Alpha,
                     const LAL::Vector& U,
                     Real Beta,
                     LAL::Vector& V)
{

/*
       Core::MltAdd<FESpace,true,true,
                   SysDim,       //DimU
                   SysDim,       //DimV
                   SysDim,       //DimI
                   FESpace::Dim, //DimJ
                   ParameterField,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Gradient>(aFESpace,aParameterField,Alpha,U,Beta,V); */

       Core::MltAdd<FESpace,FESpace,true,true,true,
                   Core::QuadratureOpt::Default,
                   SysDim,       //DimU
                   SysDim,       //DimV
                   SysDim,       //DimI
                   FESpace::Dim, //DimJ
                   ParameterField,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Gradient>(aFESpace,aFESpace,aParameterField,Alpha,U,Beta,V);
}

 
} //OndoMathX
// -----------------------------------------------------------------------------------------//


