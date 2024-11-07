#pragma once



// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{

namespace GGL
{

template<Index SysDim = 1, class FESpaceMaster, class FESpaceSlave>
void CreateMortarMass(FESpaceMaster & aFESpaceM,
                      FESpaceSlave & aFESpaceS, 
                      LAL::SparseMatrix & CM, LAL::SparseMatrix & CS)
{
    assert(FESpaceMaster::Dim==1);
    assert(FESpaceSlave::Dim==1);
    
    const GaussLobattoElement & GLE_M = aFESpaceM.getFE();
    const GaussLobattoElement & GLE_S = aFESpaceS.getFE();
    
    GaussLobattoElement GE;
    
    if (GLE_M.getNI()>2)
        GE = GaussLobattoElement(GLE_M.getNI()-2,0,0,false);
    else
        GE = GaussLobattoElement(1,0,0,false);
  
    const Mesh & meshM = aFESpaceM.getMesh();
    const Mesh & meshS = aFESpaceS.getMesh();
    
    std::vector<Index> intersection_matrix = intersectionMatrix(meshS,meshM);
    
    //Allocating CM and CS
    LAL::Allocate(CM, SysDim*aFESpaceM.getNumDoFs(), SysDim*meshM.getNumElements()*GE.getNumPoints());
    LAL::Allocate(CS, SysDim*aFESpaceS.getNumDoFs(), SysDim*meshM.getNumElements()*GE.getNumPoints());
    
    
    // Assembling of CM
    for (Index el_M = 0; el_M < meshM.getNumElements(); ++el_M)
    {
        for (Index I=0;I<GLE_M.getNumPoints();++I)
        {
            Index iGlob = aFESpaceM.Loc2Glob(el_M, I);
            RealVector PGL = GLE_M.getPointCoordinate(I);
            Real jacobian = aFESpaceM.getJacobian(el_M, I);
            Real weight = GLE_M.getQuadratureWeight(I);
            
            for (Index J=0;J<GE.getNumPoints();++J)
            {
                Index jGlob = GE.getNumPoints()*el_M+J;
                
                Real value = jacobian*weight*GE.getValueBasisFunc(J,PGL);
                
                if (value != 0.0) for (Index k = 0; k < SysDim; ++k)
                    LAL::AddInteraction(CM, SysDim*iGlob+k, SysDim*jGlob+k,value);
        
            }
        }
    }
    
    
    
    
    // Assembling of CS
    for (Index el_S = 0; el_S < meshS.getNumElements(); el_S++)
    {
        RealVector PGL_uvw;
        RealVector G_uvw;
        RealVector PGL_xyz;
        
        Index el_M = intersection_matrix[el_S];
        
        for (Index I=0;I<GLE_S.getNumPoints();++I)
        {
            Index iGlob = aFESpaceS.Loc2Glob(el_S, I);
            
            PGL_uvw = GLE_S.getPointCoordinate(I);
            meshS.getElement(el_S).uvw2xyz(PGL_uvw,PGL_xyz);
            meshM.getElement(el_M).xyz2uvw(PGL_xyz,G_uvw);
            
            Real jacobian = aFESpaceS.getJacobian(el_M, I);
            Real weight = GLE_S.getQuadratureWeight(I);
            
            for (Index J=0;J<GE.getNumPoints();++J)
            {
                Index jGlob = GE.getNumPoints()*el_M+J;
                
                Real value = jacobian*weight
                            *GE.getValueBasisFunc(J,G_uvw);
                
                if (value != 0.0) for (Index k = 0; k < SysDim; ++k)
                    LAL::AddInteraction(CS, SysDim*iGlob+k, SysDim*jGlob+k,value);
            }
        }
    }
}








} //GGL

} // Ondolab
// -----------------------------------------------------------------------------------------//




