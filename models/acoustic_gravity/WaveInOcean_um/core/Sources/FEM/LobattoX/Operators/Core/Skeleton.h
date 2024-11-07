#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 
// ---------------------------------------------------------------------------------//
template<class FESkeletonSpace>
void AssembleNormalJump(FESkeletonSpace & aFESkeletonSpace,
                     LAL::SparseMatrix & aMatrix,
                     Real Alpha
                     )
{
    constexpr Index SysDim = FESkeletonSpace::Dim;

    std::array<Real,SysDim> n;
 
    const auto & FE = aFESkeletonSpace.getFE();

    //Loop on skeleton
    for (Index iElt= 0; iElt < aFESkeletonSpace.getNumElements(); iElt++)
    {
        
        for (Index iLoc = 0; iLoc < FE.getNumPoints(); iLoc++)
        {
            Index iEltK;
            Index iEltL;
            Index iLocK;
            Index iLocL;

            Real Weight = FE.getQuadratureWeight(iLoc);

            aFESkeletonSpace.getNeighbours(iElt,iEltK,iEltL);
            iLocK = aFESkeletonSpace.Loc2Loc(iElt,iLoc,true);
            iLocL = aFESkeletonSpace.Loc2Loc(iElt,iLoc,false);

            Index IdxK = (aFESkeletonSpace.getFESpace().template Loc2Glob<false>(iEltK,iLocK))*SysDim;
            Index IdxL = (aFESkeletonSpace.getFESpace().template Loc2Glob<false>(iEltL,iLocL))*SysDim;
 
            aFESkeletonSpace.getNormal(iElt,iLoc,n);
            Real J = ArrayAlgebra::Norm(n);

            if constexpr(SysDim==2)
            {
                Real tmp;

                tmp = Alpha*Weight*n[0]*n[0]/J;
                LAL::AddInteraction(aMatrix,IdxK,IdxK,tmp);
                LAL::AddInteraction(aMatrix,IdxL,IdxL,tmp);
                LAL::AddInteraction(aMatrix,IdxK,IdxL,-tmp);
                LAL::AddInteraction(aMatrix,IdxL,IdxK,-tmp);

                tmp = Alpha*Weight*n[0]*n[1]/J;
                LAL::AddInteraction(aMatrix,IdxK,IdxK+1,tmp);
                LAL::AddInteraction(aMatrix,IdxL,IdxL+1,tmp);
                LAL::AddInteraction(aMatrix,IdxK,IdxL+1,-tmp);
                LAL::AddInteraction(aMatrix,IdxL,IdxK+1,-tmp);

                LAL::AddInteraction(aMatrix,IdxK+1,IdxK,tmp);
                LAL::AddInteraction(aMatrix,IdxL+1,IdxL,tmp);
                LAL::AddInteraction(aMatrix,IdxK+1,IdxL,-tmp);
                LAL::AddInteraction(aMatrix,IdxL+1,IdxK,-tmp);

                tmp = Alpha*Weight*n[1]*n[1]/J;
                LAL::AddInteraction(aMatrix,IdxK+1,IdxK+1,tmp);
                LAL::AddInteraction(aMatrix,IdxL+1,IdxL+1,tmp);
                LAL::AddInteraction(aMatrix,IdxK+1,IdxL+1,-tmp);
                LAL::AddInteraction(aMatrix,IdxL+1,IdxK+1,-tmp);
            }  
        }
    }
}



// Assemble matrix : aMarix = A (from centered flux) + delta*S_h (from decentered flux)

template<Index SysDim,class FESkeletonSpace>
void AssembleCenteredFlux(FESkeletonSpace & aFESkeletonSpace,
                          std::function<void(const RealVector & x,
                                             const std::array<Real, FESkeletonSpace::Dim> & n,
                                             std::array<std::array<Real, SysDim>, SysDim> & A)> flux_op,
                          LAL::SparseMatrix & aMatrix,
                          Real Alpha
                        )
{
    constexpr Index Dim = FESkeletonSpace::Dim;

    std::array<Real,Dim> n;
    RealVector x;
    std::array<std::array<Real, SysDim>, SysDim> A;
 
    const auto & FE = aFESkeletonSpace.getFE();

    //Loop on skeleton
    for (Index iElt= 0; iElt < aFESkeletonSpace.getNumElements(); iElt++)
    {
        Index iEltK;
        Index iEltL;
        aFESkeletonSpace.getNeighbours(iElt,iEltK,iEltL);
        
        for (Index iLoc = 0; iLoc < FE.getNumPoints(); iLoc++)
        {
            Index iLocK;
            Index iLocL;

            Real Weight = FE.getQuadratureWeight(iLoc);

            iLocK = aFESkeletonSpace.Loc2Loc(iElt,iLoc,true);
            iLocL = aFESkeletonSpace.Loc2Loc(iElt,iLoc,false);

            Index IdxK = (aFESkeletonSpace.getFESpace().Loc2GlobDisc(iEltK,iLocK))*SysDim;
            Index IdxL = (aFESkeletonSpace.getFESpace().Loc2GlobDisc(iEltL,iLocL))*SysDim;
 
            aFESkeletonSpace.getNormal(iElt,iLoc,n);
            Real J = ArrayAlgebra::Norm(n);
            for (Index d=0; d<Dim;++d) n[d]/=J;

            aFESkeletonSpace.getDoFCoordinate(iElt,iLoc,x);

            flux_op(x,n,A);

            for (Index u=0; u<SysDim;++u)
            {
                for (Index v=0; v<SysDim;++v)
                {
                    Real tmp = 0.5*Alpha*Weight*A[u][v]*J;

                    LAL::AddInteraction(aMatrix,IdxK+u,IdxL+v,tmp);
                    LAL::AddInteraction(aMatrix,IdxL+u,IdxK+v,-tmp);

                } 
            } 
        }
    }
}


template< Index SysDim,class FESkeletonSpace>
void AssembleJumpFlux(FESkeletonSpace & aFESkeletonSpace,
                          std::function<void(const RealVector & x,
                                             const std::array<Real, FESkeletonSpace::Dim> & n,
                                             std::array<std::array<Real, SysDim>, SysDim> & A)> flux_op,
                          LAL::SparseMatrix & aMatrix,
                          Real Alpha
                        )
{
    constexpr Index Dim = FESkeletonSpace::Dim;

    std::array<Real,Dim> n;
    RealVector x;
    std::array<std::array<Real, SysDim>, SysDim> A;

    const auto & FE = aFESkeletonSpace.getFE();

    //Loop on skeleton
    for (Index iElt= 0; iElt < aFESkeletonSpace.getNumElements(); iElt++)
    {
        Index iEltK;
        Index iEltL;
        aFESkeletonSpace.getNeighbours(iElt,iEltK,iEltL);
        
        for (Index iLoc = 0; iLoc < FE.getNumPoints(); iLoc++)
        {
            Index iLocK;
            Index iLocL;

            Real Weight = FE.getQuadratureWeight(iLoc);

            iLocK = aFESkeletonSpace.Loc2Loc(iElt,iLoc,true);
            iLocL = aFESkeletonSpace.Loc2Loc(iElt,iLoc,false);

            Index IdxK = (aFESkeletonSpace.getFESpace().Loc2GlobDisc(iEltK,iLocK))*SysDim;
            Index IdxL = (aFESkeletonSpace.getFESpace().Loc2GlobDisc(iEltL,iLocL))*SysDim;
 
            aFESkeletonSpace.getNormal(iElt,iLoc,n);
            Real J = ArrayAlgebra::Norm(n);
            for (Index d=0; d<Dim;++d) n[d]/=J;

            aFESkeletonSpace.getDoFCoordinate(iElt,iLoc,x);

            flux_op(x,n,A);

            for (Index u=0; u<SysDim;++u)
            {
                for (Index v=0; v<SysDim;++v)
                {
                    Real tmp = Alpha*Weight*A[u][v]*J;

                    LAL::AddInteraction(aMatrix,IdxK+u,IdxK+v,tmp);
                    LAL::AddInteraction(aMatrix,IdxK+u,IdxL+v,-tmp);
                    LAL::AddInteraction(aMatrix,IdxL+u,IdxK+v,-tmp);
                    LAL::AddInteraction(aMatrix,IdxL+u,IdxL+v,tmp);
                } 
            } 
        }
    }
}




 
} //OndoMathX
// -----------------------------------------------------------------------------------------//


