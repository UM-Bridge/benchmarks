#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

template<class FESpace>
void MltAddGradient(FESpace & aFESpace,
                Real Alpha,
                LAL::Vector& U,
                Real Beta,
                LAL::Vector& V)
{
       Core::MltAdd<FESpace,FESpace,true,true,false,
                   Core::QuadratureOpt::Default,
                   1,               //DimU
                   FESpace::Dim,    //DimV
                   1,               //DimI
                   FESpace::Dim,    //DimJ
                   Field::Identity,
                   Core::DiffOp::Gradient,
                   Core::DiffOp::Identity>(aFESpace,aFESpace,Field::IdentityField,Alpha,U,Beta,V);
}

template<class FESpace>
void TransposeMltAddGradient(FESpace & aFESpace,
                Real Alpha,
                LAL::Vector& U,
                Real Beta,
                LAL::Vector& V)
{

       Core::MltAdd<FESpace,FESpace,true,false,true,
                   Core::QuadratureOpt::Default,
                   FESpace::Dim,        //DimU
                   1,                   //DimV
                   1,                   //DimI
                   FESpace::Dim,        //DimJ
                   Field::Identity,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Gradient>(aFESpace,aFESpace,Field::IdentityField,Alpha,U,Beta,V);
}





 
} //OndoMathX
// -----------------------------------------------------------------------------------------//


