#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 

template<class FESpace>
void MltAddDivergence(FESpace & aFESpace,
                     Real Alpha,
                     LAL::Vector& U,
                     Real Beta,
                     LAL::Vector& V)
{
       Core::MltAdd<FESpace,FESpace,true,true,false,
                   Core::QuadratureOpt::Default,
                   FESpace::Dim,    //DimU
                   1,               //DimV
                   1,               //DimI
                   1,               //DimJ
                   Field::Identity,
                   Core::DiffOp::Divergence,
                   Core::DiffOp::Identity>(aFESpace,Field::IdentityField,Alpha,U,Beta,V);
}

template<class FESpace_Displ, class FESpace_Press>
void MltAddDivergence(FESpace_Displ & aFESpaceDispl,
                      FESpace_Press & aFESpacePress,
                      Real Alpha,
                      LAL::Vector& Displ,
                      Real Beta,
                      LAL::Vector& Press)
{
       Core::MltAdd<FESpace_Displ,FESpace_Press,false,true,false,
                   Core::QuadratureOpt::Default,
                   FESpace_Displ::Dim,    //DimU
                   1,                     //DimV
                   1,                     //DimI
                   1,                     //DimJ
                   Field::Identity,
                   Core::DiffOp::Divergence,
                   Core::DiffOp::Identity>(aFESpaceDispl,aFESpacePress,Field::IdentityField,Alpha,Displ,Beta,Press);
}


template<class FESpace>
void TransposeMltAddDivergence(FESpace & aFESpace,
                            Real Alpha,
                            LAL::Vector& U,
                            Real Beta,
                            LAL::Vector& V)
{
       Core::MltAdd<FESpace,FESpace,true,false,true,
                   Core::QuadratureOpt::Default,
                   1,               //DimU
                   FESpace::Dim,    //DimV
                   1,               //DimI
                   1,               //DimJ
                   Field::Identity,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Divergence>(aFESpace,Field::IdentityField,Alpha,U,Beta,V);
}



template<class FESpace_Displ, class FESpace_Press>
void TransposeMltAddDivergence(FESpace_Displ & aFESpaceDispl,
                               FESpace_Press & aFESpacePress,
                               Real Alpha,
                               LAL::Vector& Press,
                               Real Beta,
                               LAL::Vector& Displ)
{
       Core::MltAdd<FESpace_Press,FESpace_Displ,false,false,true,
                   Core::QuadratureOpt::Trial,
                   1,                     //DimU
                   FESpace_Displ::Dim,    //DimV
                   1,                     //DimI
                   1,                     //DimJ
                   Field::Identity,
                   Core::DiffOp::Identity,
                   Core::DiffOp::Divergence>(aFESpacePress,aFESpaceDispl,Field::IdentityField,Alpha,Press,Beta,Displ);
}




 
} //OndoMathX
// -----------------------------------------------------------------------------------------//


