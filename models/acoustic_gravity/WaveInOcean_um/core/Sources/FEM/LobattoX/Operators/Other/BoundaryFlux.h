#pragma once



// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{



template<Index SysDim,class TraceSpace>
void AssembleBndyCenteredFlux(TraceSpace & aTraceSpace,
                              std::function<void(const RealVector & x,
                                             const std::array<Real, TraceSpace::Dim+1> & n,
                                             std::array<std::array<Real, SysDim>, SysDim> & A)> bndy_func,
                              LAL::SparseMatrix & aMatrix,
                              Real Alpha 
                        )
{
    std::array<std::array<Real, SysDim>, SysDim> B;
    std::array<Real,TraceSpace::Dim+1> n;
    RealVector x;

    aTraceSpace.setAsTraceSpace();

    const GaussLobattoElement & GLE = aTraceSpace.getFE();


    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        Index iEltGlob = aTraceSpace.getVolumeElement(iElt);

        for (Index iLoc = 0; iLoc < GLE.getNumPoints(); iLoc++)
        {
            aTraceSpace.getNormal(iElt,iLoc,n);
            aTraceSpace.getDoFCoordinate(iElt,iLoc,x);

            Real J = ArrayAlgebra::Norm(n);
            for (Index d=0; d<SysDim;++d) n[d]/=J;

            Real Weight = GLE.getQuadratureWeight(iLoc);

            Index iLocVol = aTraceSpace.Loc2Loc(iLoc);
            Index Idx = aTraceSpace.getFESpace().Loc2GlobDisc(iEltGlob,iLocVol)*SysDim;

            //cout << iElt << " " << iLoc << "///" << iEltGlob << " " << iLocVol << " " << Idx << endl;

            // Get the matrix B at position x (n is the normal) 
            bndy_func(x,n,B);

            // Assembling procedure
            for (Index u=0; u<SysDim;++u)
            {
                for (Index v=0; v<SysDim;++v)
                {
                    Real tmp = Alpha*Weight*B[u][v]*J;
                    //cout << "tmp : "<<tmp<<",u : "<<u<<",v : "<<v<< endl;
                    LAL::AddInteraction(aMatrix,Idx+u,Idx+v,tmp);               
                } 
            } 
        }
    } 
}

template<Index SysDim,class TraceSpace>
void AssembleBndyDecenteredFlux(TraceSpace & aTraceSpace,
                              std::function<void(const RealVector & x,
                                             const std::array<Real, TraceSpace::Dim+1> & n,
                                             std::array<std::array<Real, SysDim>, SysDim> & A)> bndy_func,
                              LAL::SparseMatrix & aMatrix,
                              Real Alpha 
                        )
{
    std::array<std::array<Real, SysDim>, SysDim> B;
    std::array<Real,TraceSpace::Dim+1> n;
    RealVector x;

    aTraceSpace.setAsTraceSpace();

    const GaussLobattoElement & GLE = aTraceSpace.getFE();


    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        Index iEltGlob = aTraceSpace.getVolumeElement(iElt);

        for (Index iLoc = 0; iLoc < GLE.getNumPoints(); iLoc++)
        {
            Index iLocVol = aTraceSpace.Loc2Loc(iLoc);

            aTraceSpace.getNormal(iElt,iLoc,n);
            aTraceSpace.getDoFCoordinate(iElt,iLoc,x);

            Real J = ArrayAlgebra::Norm(n);
            for (Index d=0; d<SysDim;++d) n[d]/=J;

            Real Weight = GLE.getQuadratureWeight(iLoc);

            Index Idx = aTraceSpace.getFESpace().Loc2GlobDisc(iEltGlob,iLocVol)*SysDim;

            // Get the matrix B at position x (n is the normal) 
            bndy_func(x,n,B);

            // Assembling procedure
            for (Index u=0; u<SysDim;++u)
            {
                for (Index v=0; v<SysDim;++v)
                {
                    Real tmp = Alpha*Weight*B[u][v]*J;
                    LAL::AddInteraction(aMatrix,Idx+u,Idx+v,tmp);               
                } 
            } 
        }
    } 
}








 
} // OndoMathX
// -----------------------------------------------------------------------------------------//




