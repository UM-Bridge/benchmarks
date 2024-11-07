#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

namespace GGL {

template<class TraceSpace>
void CreateBndySourceNormal(TraceSpace & aTraceSpace,
                          std::function<Real(const RealVector&)> aFunc,
                          LAL::Vector& U)
{
    constexpr Index SysDim = TraceSpace::Dim+1;
    
    aTraceSpace.setAsTraceSpace();
    LAL::Allocate(U,aTraceSpace.getNumDoFs()*SysDim);
      
    const auto & FE = aTraceSpace.getFE();
    
    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < FE.getNumPoints() ; iLoc++)
        {
            Index iGlob = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            std::array<Real,SysDim> n;
            aTraceSpace.getNormal(iElt,iLoc,n);
            
            
            Real weight = FE.getQuadratureWeight(iLoc);
            
            RealVector xyz;
            aTraceSpace.getDoFCoordinate(iElt,iLoc,xyz);
            
            Real function_value = aFunc(xyz);
            
            for (Index jLoc=0;jLoc<SysDim;jLoc++)
                U(SysDim*iGlob+jLoc) += n[jLoc]*weight*function_value;
        }
    }
    
}

/*
template<class TraceSpace>
void CreateBndyScalarSource(TraceSpace & aTraceSpace,
                          std::function<Real(const RealVector&)> aFunc,
                          LAL::Vector& U)
{   
    LAL::Allocate(U,aTraceSpace.getNumDoFs());
      
    const GaussLobattoElement & GLE = aTraceSpace.getFE();
    
    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < GLE.getNumPoints() ; iLoc++)
        {
            Index iGlob = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            Real jacobian = aTraceSpace.getJacobian(iElt,iLoc);
            Real weight = GLE.getQuadratureWeight(iLoc);
            
            RealVector xyz;
            aTraceSpace.getDoFCoordinate(iElt,iLoc,xyz);
            
            Real function_value = aFunc(xyz);
            
             U(iGlob) += weight*jacobian*function_value;
        }
    }
    
}

template<Index SysDim,class TraceSpace>
void CreateBndySource(TraceSpace & aTraceSpace,
                        std::function<std::array<Real,SysDim>(const RealVector&)> aFunc,
                        LAL::Vector& U)
{
    LAL::Allocate(U,SysDim*aTraceSpace.getNumDoFs());
      
    const GaussLobattoElement & GLE = aTraceSpace.getFE();
    
    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < GLE.getNumPoints() ; iLoc++)
        {
            Index iGlob = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            Real jacobian = aTraceSpace.getJacobian(iElt,iLoc);
            Real weight = GLE.getQuadratureWeight(iLoc);
            
            RealVector xyz;
            aTraceSpace.getDoFCoordinate(iElt,iLoc,xyz);
            
            std::array<Real,SysDim> function_value = aFunc(xyz);
            
            for (Index i=0;i<SysDim;++i)
                U(iGlob*SysDim+i) += weight*jacobian*function_value[i];
        }
    }
    
}
*/

} //GGL
}//OndoMathX
