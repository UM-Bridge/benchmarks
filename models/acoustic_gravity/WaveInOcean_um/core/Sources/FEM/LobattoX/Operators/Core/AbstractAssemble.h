#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
namespace Core {

template<
    class FESpaceU,  
    class FESpaceV,  
    bool Same_FESpace,
    bool ContinuousU, 
    bool ContinuousV,
    QuadratureOpt QuadOpt,
    Index DimU,
    Index DimV,
    Index DimI,
    Index DimJ,
    class ParameterField_C,
    class SparseMatrix,
    DiffOp E1,
    DiffOp F1,
    class ParameterField_L1 = OndoMathX::Field::Identity,
    class ParameterField_M1 = OndoMathX::Field::Identity,
    DiffOp E2 = Zero,
    DiffOp F2 = Zero,
    class ParameterField_L2 = OndoMathX::Field::Identity,
    class ParameterField_M2 = OndoMathX::Field::Identity,
    DiffOp E3 = Zero,
    DiffOp F3 = Zero,
    class ParameterField_L3 = OndoMathX::Field::Identity,
    class ParameterField_M3 = OndoMathX::Field::Identity,
    class FEQuad = DefaultFE
    >
    void Assemble(FESpaceU & aFESpaceU,
                  FESpaceV & aFESpaceV,
            ParameterField_C & aParameterField_C,
            SparseMatrix & aMatrix,
            Real Alpha,
            ParameterField_L1 & aParameterField_L1 = Field::IdentityField,
            ParameterField_M1 & aParameterField_M1 = Field::IdentityField,
            ParameterField_L2 & aParameterField_L2 = Field::IdentityField,
            ParameterField_M2 & aParameterField_M2 = Field::IdentityField,
            ParameterField_L3 & aParameterField_L3 = Field::IdentityField,
            ParameterField_M3 & aParameterField_M3 = Field::IdentityField,
            FEQuad & FE_Q = _DefaultFE)
            {
                constexpr bool lumped_matrix = 
                             (   E1 == Identity
                            &&   F1 == Identity
                            &&   (E2 == Zero || E2 == Identity) 
                            &&   (F2 == Zero || F2 == Identity) 
                            &&   (E3 == Zero || E3 == Identity) 
                            &&   (F3 == Zero || F3 == Identity));

                //Todo : use the knowledge that the matrix is lumped to simplify some loops

                const auto & GLE_U = aFESpaceU.getFE();
                const auto & GLE_V = aFESpaceV.getFE();

                Index NLocU = GLE_U.getNumPoints();
                Index NLocV = GLE_V.getNumPoints();
                Index NElem = aFESpaceV.getNumElements();

                LAL::Vector U;
                LAL::Vector V;

                Index nU = aFESpaceU.getNumDoFsDisc()*DimU;
                Index nV = aFESpaceV.getNumDoFsDisc()*DimV; 

                LAL::Allocate(U,nU);
                LAL::Allocate(V,nV);

                Real * pU = LAL::getData(U);
                Real * pV = LAL::getData(V);

                // Loop on the dimension of the unknowns
                for (Index u=0; u<DimU; ++u)
                {
                    // Loop on local number of dof
                    for (Index iLoc=0; iLoc<NLocU; ++iLoc)   
                    {
                        {
                            //Loop on the elements to assign 1 to the corresponding dof
                            for (Index iElt=0;iElt<NElem;++iElt)
                            {
                                Index IdxU = aFESpaceU.Loc2GlobDisc(iElt,iLoc)*DimU+u;
                                pU[IdxU] = Alpha;
                            }
                            
                            //MltAdd in order to recover one row of the elementary matrix of each elements
                            Core::MltAdd<FESpaceU,FESpaceV,Same_FESpace,false,false,
                            QuadOpt,
                            DimU,
                            DimV,
                            DimI,
                            DimJ,
                            ParameterField_C,
                            E1,F1,ParameterField_L1,ParameterField_M1,
                            E2,F2,ParameterField_L2,ParameterField_M2,
                            E3,F3,ParameterField_L3,ParameterField_M3,
                            false,
                            FEQuad> (aFESpaceU, aFESpaceV, aParameterField_C, 1.0, U, 0.0, V,
                                                                        aParameterField_L1, aParameterField_M1,
                                                                        aParameterField_L2, aParameterField_M2,
                                                                        aParameterField_L3, aParameterField_M3,
                                                                        FE_Q);
                            
                            //Loop on the elements and assembling procedure
                            for (Index iElt=0;iElt<NElem;++iElt)
                            {
                                Index IdxU;

                                //Pre-compute the global indices.
                                if constexpr(ContinuousU)
                                {
                                    IdxU = aFESpaceU.Loc2Glob(iElt,iLoc)*DimU+u;
                                }
                                else
                                {
                                    IdxU = aFESpaceU.Loc2GlobDisc(iElt,iLoc)*DimU+u;
                                }

                                for (Index jLoc=0; jLoc<NLocV; ++jLoc)
                                {
                                    for (Index v=0; v<DimV; ++v)
                                    {
                                        Index IdxV;
                                        Index IdxV_D;
                                        
                                        IdxV_D = aFESpaceV.Loc2GlobDisc(iElt,jLoc)*DimV+v;
                                        
                                        if constexpr(ContinuousV == true)
                                            IdxV = aFESpaceV.Loc2Glob(iElt,jLoc)*DimV+v;
                                        else
                                            IdxV = IdxV_D;
                                        
                                        Real value = pV[IdxV_D];
                                        
                                        if (value != 0.0)
                                        {
                                            LAL::AddInteraction(aMatrix,IdxV,IdxU,value);
                                        }
                                    }
                                }
                            }
                            
                            //Loop on the elements to assign 1 to the corresponding dof
                            for (Index iElt=0;iElt<NElem;++iElt)
                            {
                                Index IdxU = aFESpaceU.Loc2GlobDisc(iElt,iLoc)*DimU+u;
                                pU[IdxU] = 0.0;
                            }
                        }
                    }
                }
            }




} //Core
 
}//OndoMathX
