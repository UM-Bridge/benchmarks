#pragma once



// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{



/*


template<class TraceSpace, class SparseMatrix>
void AssembleBndyMassNormal(TraceSpace aTraceSpace, SparseMatrix& aMatrix, Real Coef)
{
    constexpr Index SysDim = TraceSpace::Dim+1;
    
    aTraceSpace.setAsTraceSpace();
    Index n_tracespace = aTraceSpace.getNumDoFs();
    
    const GaussLobattoElement & GLE = aTraceSpace.getFE();

    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < GLE.getNumPoints() ; iLoc++)
        {
            Index iGlob_tracespace = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            std::array<Real,SysDim> n;
            aTraceSpace.getNormal(iElt,iLoc,n);

            Real J = ArrayAlgebra::Norm(n);
            
            Real weight_J = GLE.getQuadratureWeight(iLoc)/J;
    
            for (Index id=0;id<SysDim;id++)
            {
                for (Index jd=0;jd<SysDim;jd++)
                {
                    Real value = Coef*n[id]*n[jd]*weight_J;

                    if (abs(value) > Zero_Machine)
                        LAL::AddInteraction(aMatrix,
                                        SysDim*iGlob_tracespace+id,
                                        SysDim*iGlob_tracespace+jd,
                                        value);
                }
            }
        }
    }
}

 

//Careful here the matrix is rectangular
template<class TraceSpace, class SparseMatrix>
void CreateBndyMassNormal(TraceSpace aTraceSpace, SparseMatrix& aMatrix, Real Coef)
{
    constexpr Index SysDim = TraceSpace::Dim+1;
    
    aTraceSpace.setAsFESpace();
    Index n_fespace = aTraceSpace.getNumDoFs();
    
    aTraceSpace.setAsTraceSpace();
    Index n_tracespace = aTraceSpace.getNumDoFs();
    
    const GaussLobattoElement & GLE = aTraceSpace.getFE();

    LAL::Allocate(aMatrix, n_fespace, SysDim*n_tracespace);
    
    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < GLE.getNumPoints() ; iLoc++)
        {
            aTraceSpace.setAsFESpace();
            Index iGlob_fespace = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            aTraceSpace.setAsTraceSpace();
            Index iGlob_tracespace = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            std::array<Real,SysDim> n;
            aTraceSpace.getNormal(iElt,iLoc,n);
            
           
            Real weight = GLE.getQuadratureWeight(iLoc);
 
            for (Index jLoc=0;jLoc<SysDim;jLoc++)
                LAL::AddInteraction(aMatrix,
                                    iGlob_fespace,
                                    SysDim*iGlob_tracespace+jLoc,
                                    Coef*n[jLoc]*weight);
        
        }
    }
}



template<class TraceSpace>
void AssembleBndyMass(TraceSpace & aTraceSpace, LAL::DiagonalMatrix & aMatrix, Real Coef)
{
    aTraceSpace.setAsTraceSpace();
    Index n_tracespace = aTraceSpace.getNumDoFs();

    const GaussLobattoElement & GLE = aTraceSpace.getFE();

    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < GLE.getNumPoints() ; iLoc++)
        {
            Index iGlob_tracespace = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            std::array<Real,TraceSpace::Dim+1> n;
            aTraceSpace.getNormal(iElt,iLoc,n);
            Real J = ArrayAlgebra::Norm(n);
           
            Real weight = GLE.getQuadratureWeight(iLoc);
 
            LAL::AddInteraction(aMatrix,
                                    iGlob_tracespace,
                                    iGlob_tracespace,
                                    Coef*J*weight);
        
        }
    }
}

template<class TraceSpace, class SparseMatrix>
void CreateBndyMass(TraceSpace aTraceSpace, SparseMatrix& aMatrix, Real Coef)
{    
    aTraceSpace.setAsFESpace();
    Index n_fespace = aTraceSpace.getNumDoFs();

    LAL::Allocate(aMatrix, n_fespace, n_fespace);
    
    AssembleBndyMass(aTraceSpace,aMatrix,Coef);
}



*/




template<class TraceSpace, class SparseMatrix>
void AssembleBndyMassTangent(TraceSpace aTraceSpace, SparseMatrix& aMatrix, Real Coef)
{
    constexpr Index SysDim = TraceSpace::Dim+1;
    assert(SysDim == 2);
    
    aTraceSpace.setAsTraceSpace();
    Index n_tracespace = aTraceSpace.getNumDoFs();
    
    const auto & FE = aTraceSpace.getFE();

    for (Index iElt = 0; iElt < aTraceSpace.getNumElements(); iElt++)
    {
        for (Index iLoc = 0; iLoc < FE.getNumPoints() ; iLoc++)
        {
            Index iGlob_tracespace = aTraceSpace.Loc2Glob(iElt,iLoc);
            
            std::array<Real,SysDim> n;
            std::array<Real,SysDim> t;
            aTraceSpace.getNormal(iElt,iLoc,n);

            Real J = ArrayAlgebra::Norm(n);
            
            t[0] = n[1];
            t[1] = -n[0];


            Real weight_J = FE.getQuadratureWeight(iLoc)/J;
    
            for (Index id=0;id<SysDim;id++)
            {
                for (Index jd=0;jd<SysDim;jd++)
                {
                    Real value = Coef*t[id]*t[jd]*weight_J;

                    if (abs(value) > Zero_Machine)
                        LAL::AddInteraction(aMatrix,
                                        SysDim*iGlob_tracespace+id,
                                        SysDim*iGlob_tracespace+jd,
                                        value);
                }
            }
        }
    }
}



 
} // OndoMathX
// -----------------------------------------------------------------------------------------//




