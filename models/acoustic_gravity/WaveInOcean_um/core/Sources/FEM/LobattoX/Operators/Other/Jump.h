#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

namespace GGL {

// ---------------------------------------------------------------------------------//
template<class FESkeletonSpace>
void AssembleNormalJump(FESkeletonSpace & aFESkeletonSpace,
                     LAL::SparseMatrix & aMatrix,
                     Real Alpha
                     )
{
    //allocating matrix
    constexpr Index SysDim = FESkeletonSpace::Dim;

    std::array<Real,SysDim> n;
 
    const GaussLobattoElement & GLE = aFESkeletonSpace.getFE();

    //Loop on skeleton
    for (Index iElt= 0; iElt < aFESkeletonSpace.getNumElements(); iElt++)
    {
        
        for (Index iLoc = 0; iLoc < GLE.getNumPoints(); iLoc++)
        {
            Index iEltK;
            Index iEltL;
            Index iLocK;
            Index iLocL;

            Real Weight = GLE.getQuadratureWeight(iLoc);

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




} //GGL
} //OndoMathX
// -----------------------------------------------------------------------------------------//


