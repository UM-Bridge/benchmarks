#pragma once

 
namespace OndoMathX {

 
template<Index SysDim = 1,    
         bool Transpose = false,
         class FESpace>
    void InterpolationDG(GaussLobattoInterpolator & aGLI,
                FESpace & aFESpaceIn,
                FESpace & aFESpaceOut,
                const LAL::Vector& U,
                LAL::Vector& V)
    {
        const Real * dataU = LAL::getData(U);
        Real * dataV = LAL::getData(V);
        
        const GaussLobattoElement & GLEIn = aFESpaceIn.getFE();
        const GaussLobattoElement & GLEOut = aFESpaceOut.getFE();
        
        Index N_Loc_In = GLEIn.getNumPoints();
        Index N_Loc_Out = GLEOut.getNumPoints();

        {
#pragma omp parallel
            {
                // Temporary variables for storing local solution.
                std::array<std::vector<Real>,SysDim> TmpSol_In;
                for (Index d = 0; d < SysDim; d++)
                    TmpSol_In[d].resize(N_Loc_In);
                
                // Temporary variables for storing local projected solution.
                std::array<std::vector<Real>,SysDim> TmpSol_Out;
                for (Index d = 0; d < SysDim; d++)
                    TmpSol_Out[d].resize(N_Loc_Out);

                // Temporary variables for storing global index.
                std::vector<Index> IdxGlobIn(N_Loc_In);
                std::vector<Index> IdxGlobOut(N_Loc_Out);
        
                Index NumElem = aFESpaceIn.getNumElements();

                // Loop on elements.
#pragma omp for
                for (Index iElt = 0; iElt < NumElem; iElt++)
                {
                    //Pre-compute the global indices.
                    aFESpaceIn.Loc2GlobDisc(iElt,IdxGlobIn);
                    aFESpaceOut.Loc2GlobDisc(iElt,IdxGlobOut);

                    if constexpr(!Transpose)
                    {
                        for (Index iLoc = 0; iLoc < N_Loc_In ; iLoc++)
                        {
                            Index Idx =  IdxGlobIn[iLoc];
                            for (Index d = 0; d < SysDim; d++)
                            {
                                TmpSol_In[d][iLoc] = dataU[Idx];  
                                Idx++;
                            }
                        }

                        for (Index d = 0; d < SysDim; d++)
                        {
                            aGLI.Interpolation(TmpSol_In[d],TmpSol_Out[d]);
                        }

                        for (Index iLoc = 0; iLoc < N_Loc_Out ; iLoc++)
                        {
                            Index Idx =  IdxGlobOut[iLoc];
                            for (Index d = 0; d < SysDim; d++)
                            {
                                dataV[Idx] = TmpSol_Out[d][iLoc];
                                Idx++;
                            }
                        }
                    }

                    if constexpr(Transpose)
                    {
                        for (Index iLoc = 0; iLoc < N_Loc_Out ; iLoc++)
                        {
                            Index Idx =  IdxGlobOut[iLoc];
                            for (Index d = 0; d < SysDim; d++)
                            {
                                TmpSol_Out[d][iLoc] = dataU[Idx];  
                                Idx++;
                            }
                        }

                        for (Index d = 0; d < SysDim; d++)
                        {
                            aGLI.TransposeInterpolation(TmpSol_Out[d],TmpSol_In[d]);
                        }

                        for (Index iLoc = 0; iLoc < N_Loc_In ; iLoc++)
                        {

                            Index Idx =  IdxGlobIn[iLoc];
                            for (Index d = 0; d < SysDim; d++)
                            {
                                dataV[Idx] = TmpSol_In[d][iLoc];
                                Idx++;
                            }
                        }
                    }

                }
            }
        }

         
    }

 



} // OndoMathX
// -----------------------------------------------------------------------------------------//



     
