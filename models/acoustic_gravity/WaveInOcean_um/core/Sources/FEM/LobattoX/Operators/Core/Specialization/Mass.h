#pragma once



// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{
 
template<Index SysDim = 1, class FESpace, class ParameterField, class SparseMatrix> void AssembleMass(FESpace & aFESpace,  ParameterField & aParameterField, SparseMatrix & aMatrix, Real Alpha)
{
    Core::Assemble<FESpace,FESpace,true,true,true,
                   Core::QuadratureOpt::Default,
                   SysDim,  //DimU
                   SysDim,  //DimV
                   1,       //DimI
                   SysDim,  //DimJ
                   ParameterField,
                   SparseMatrix,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Identity>(aFESpace,aFESpace,aParameterField,aMatrix,Alpha);

}

template<Index SysDim = 1, class FESpace, class ParameterField, class SparseMatrix> void CreateMass(FESpace & aFESpace,  ParameterField & aParameterField, SparseMatrix & aMatrix, Real Alpha)
{
    LAL::Allocate(aMatrix, SysDim*aFESpace.getNumDoFs(), SysDim*aFESpace.getNumDoFs());

    AssembleMass<SysDim>(aFESpace,aParameterField,aMatrix,Alpha);
}

template<Index SysDim = 1, class FESpace, class ParameterField, class SparseMatrix> void CreateMassDG(FESpace & aFESpace,  ParameterField & aParameterField, SparseMatrix & aMatrix, Real Alpha)
{
    LAL::Allocate(aMatrix, SysDim*aFESpace.getNumDoFsDisc(), SysDim*aFESpace.getNumDoFsDisc());

    Core::Assemble<FESpace,FESpace,true,false,false,
                   Core::QuadratureOpt::Default,
                   SysDim,  //DimU
                   SysDim,  //DimV
                   1,       //DimI
                   SysDim,  //DimJ
                   ParameterField,
                   SparseMatrix,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Identity>(aFESpace,aFESpace,aParameterField,aMatrix,Alpha);

}


template<Index SysDim = 1, class FESpace, class ParameterField>
void MltAddMass(FESpace & aFESpace,
                    ParameterField & aParameterField,
                    Real Alpha,
                    LAL::Vector& U,
                    Real Beta,
                    LAL::Vector& V)
{
        Core::MltAdd<FESpace,FESpace,true,true,true,
                Core::QuadratureOpt::Default,
                SysDim,  //DimU
                SysDim,  //DimV
                1,       //DimI
                SysDim,  //DimJ
                ParameterField,
                Core::DiffOp::Identity,
                Core::DiffOp::Identity>(aFESpace,aFESpace,aParameterField,Alpha,U,Beta,V); 
}


 
} // OndoMathX
// -----------------------------------------------------------------------------------------//




