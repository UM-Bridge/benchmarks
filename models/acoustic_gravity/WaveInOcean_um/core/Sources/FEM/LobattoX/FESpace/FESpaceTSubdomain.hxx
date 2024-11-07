#pragma once
 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

template<class FESpace,Index DIM>
class FESpaceTSubdomain 
{
public: 

    static const Index Dim = DIM;

    typedef std::array<std::array<Real,Dim>,Dim> RealMatrix;

    typedef GaussLobattoInterpolator Interpolator;

private:

    const FESpace & _FESpace;
    std::vector<Index> _Elements;

    std::vector<Index> _NElemByColors;
    std::vector<std::vector<Index>> _ElemByColors;
     
    bool _isSubSpace;

    Index _NDoF;
    Index _NLocDoF;
    Index _NColor;
    Index _NElem;
   

    std::vector<Index> _glob2glob;
    std::vector<std::vector<Index>> _loc2glob;
    
public:

    FESpaceTSubdomain() = delete; 

    FESpaceTSubdomain(const FESpace &aFESpace, std::vector<Index> Elements) : _FESpace(aFESpace), _Elements(Elements)
    {
        _NElem = _Elements.size();
        _NColor = _FESpace.getNumColors();
        
        assert(_NElem>0);

        _NElemByColors.resize(_NColor,0);
        _ElemByColors.resize(_NColor);

        _isSubSpace = true;

        Index FESpaceNDoF = aFESpace.getNumDoFs();
        _NLocDoF = aFESpace.getNumLocDoFs();

     

        //Construct the glob2glob array
        _NDoF = 0;
        std::vector<Index> Space2SubSpace(FESpaceNDoF,0);

        std::vector<bool> isInSubdomain(FESpaceNDoF,false);

        for (Index e = 0; e<_NElem;++e)
        {
            Index Elem = _Elements[e];
            Index color;
        
            color = aFESpace.EltGlobalNum2Color(e);

            _NElemByColors[color]++;
            _ElemByColors[color].push_back(e);

            for (Index i=0;i<_NLocDoF;++i)
            {
                Index dof =  aFESpace._Loc2Glob(Elem,i);

                isInSubdomain[dof] = true;
            }
        }

        for (Index dof=0;dof<FESpaceNDoF;++dof)
        {
            if (isInSubdomain[dof])
            {
                _NDoF++;
                Space2SubSpace[dof] = _NDoF;
            }
        }        
 
        _glob2glob.resize(_NDoF);

        for (Index n=0;n<FESpaceNDoF;++n)
        {
            if (Space2SubSpace[n] > 0)
            {
                _glob2glob[Space2SubSpace[n]-1] = n;
            }
        }

        //Construct the loc2glob array
        _loc2glob.resize(_NElem);

        for (Index e = 0; e<_NElem;++e)
        {
            _loc2glob[e].resize(_NLocDoF);

            Index Elem = _Elements[e];

            for (Index i=0;i<_NLocDoF;++i)
            {
                _loc2glob[e][i] = Space2SubSpace[aFESpace._Loc2Glob(Elem,i)]-1;
            }
        }

    }

    
    Index _glob_subdomain2glob(Index iDof) const
    {
        return _glob2glob[iDof];
    }

    void setAsSubSpace()
    {
        _isSubSpace = true;
    }
    
    void setAsFESpace()
    {
        _isSubSpace = false;
    }


    // ---------------------------------------------------------------------------------//

    Index Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        if (_isSubSpace) return _FESpace.Loc2Glob(_Elements[iEltGlob],iLoc);
        else return _loc2glob[iEltGlob][iLoc];
    } 

    void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        if (_isSubSpace) _FESpace.Loc2Glob(_Elements[iEltGlob],loc2glob);
        else loc2glob = _loc2glob[iEltGlob]; 
    }

    Index Loc2GlobDisc(Index iEltGlob, Index iLoc) const
    {
        if (_isSubSpace) return _FESpace.Loc2GlobDisc(_Elements[iEltGlob],iLoc);
        else return iEltGlob*_NLocDoF+iLoc;
    }
    
    void Loc2GlobDisc(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
         if (_isSubSpace) _FESpace.Loc2GlobDisc(_Elements[iEltGlob],loc2glob);
         else 
         {
            loc2glob.resize(_NLocDoF);
            
            for (Index iLoc=0; iLoc < _NLocDoF; iLoc++)
                loc2glob[iLoc] = iEltGlob*_NLocDoF+iLoc;
         }
    }

    Index getNumDoFs() const
    {
        if (_isSubSpace) return _FESpace.getNumDoFs();
        else return _NDoF;
    }

    Index getNumLocDoFs() const
    {
        return _FESpace.getNumLocDoFs();
    }

    Index getNumElements() const
    {
        return _NElem;
    }

    Index getNumDoFsDisc() const
    {
        if (_isSubSpace) return _FESpace.getNumDoFsDisc();
        return _NLocDoF*_NElem;
    }

    Real getJacobian(Index iEltGlob, Index iLoc) const
    {
        RealMatrix GradDef;
        _FESpace.getGradDef(_Elements[iEltGlob],iLoc,GradDef);
        return ArrayAlgebra::Det(GradDef);
    }

    // Get the label of the element 
    Index getLabel(const Index iEltGlob) const
    {
        return _FESpace.getLabel(_Elements[iEltGlob]);
    }

    bool isAffine(const Index iEltGlob) const
    {
        return _FESpace.isAffine(_Elements[iEltGlob]);
    }


    template <Index SysDim = 1> void MltAddExtension(Real Alpha, const LAL::Vector & in, Real Beta, LAL::Vector & out) const
    {
        for (Index n=0;n<_NDoF;++n) 
        {
            Index m = _glob2glob[n];

            for (Index d=0;d<SysDim;++d) 
                out[SysDim*m+d]=Beta*out[SysDim*m+d] + Alpha*in[SysDim*n+d];
        }
    }


    template <Index SysDim = 1> void MltAddRestriction(Real Alpha, const LAL::Vector & in, Real Beta, LAL::Vector & out) const
    {
        for (Index n=0;n<_NDoF;++n) 
        {
            Index m = _glob2glob[n];

            for (Index d=0;d<SysDim;++d) 
                out[SysDim*n+d] = Beta*out[SysDim*n+d] + Alpha*in[SysDim*m+d];
        }
    }

    template <Index SysDim = 1> void MltAddExtensionDisc(Real Alpha, const LAL::Vector & in, Real Beta, LAL::Vector & out) const
    {
        for (Index e=0;e<_NElem;++e) 
        {
            Index glob_elem = _Elements[e];

            for (Index iLoc=0; iLoc<_NLocDoF; ++iLoc)
            {
                Index iglob_fespace = SysDim*_FESpace.Loc2GlobDisc(glob_elem,iLoc);
                Index iglob_proper_fespace = SysDim*e*_NLocDoF+iLoc;


                for (Index d=0;d<SysDim;++d) 
                {
                    out[iglob_fespace+d]*=Beta;
                    out[iglob_fespace+d]+=Alpha*in[iglob_proper_fespace+d];
                }
            }
        }
    }


    template <Index SysDim = 1> void MltAddRestrictionDisc(Real Alpha, const LAL::Vector & in, Real Beta, LAL::Vector & out) const
    {
        for (Index e=0;e<_NElem;++e) 
        {
            Index glob_elem = _Elements[e];

            for (Index iLoc=0; iLoc<_NLocDoF; ++iLoc)
            {
                Index iglob_fespace = SysDim*_FESpace.Loc2GlobDisc(glob_elem,iLoc);
                Index iglob_proper_fespace = SysDim*e*_NLocDoF+iLoc;

                for (Index d=0;d<SysDim;++d) 
                {
                    out[iglob_proper_fespace+d]*=Beta;
                    out[iglob_proper_fespace+d]+=Alpha*in[iglob_fespace+d];
                }
            }
        }
    }

    void getGradDef(Index iEltGlob, Index iLoc, RealMatrix& GradDef) const
    {          
        _FESpace.getGradDef(_Elements[iEltGlob],iLoc,GradDef);
    }

    Index getNumColors() const
    {          
        return _NColor;
    }

    void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector& xyz) const
	{
        Index Elem = _Elements[iEltGlob];
		_FESpace.getDoFCoordinate(Elem,iLoc,xyz);
	}


    void getDoFCoordinateDisc(const Index iGlob, RealVector& xyz) const
    {
        assert(false);
    }

    void getDoFCoordinate(const Index iGlob, RealVector& xyz) const
    {
        Index n = _glob2glob[iGlob];
        _FESpace.getDoFCoordinate(n,xyz);
    }

    Index getNumElements(Index c) const
	{
         assert(c<_NColor);

         return _NElemByColors[c];
    }

    Index EltColorNum2EltGlobalNum(Index c, Index e) const
	{
        assert(c<_NColor);

        return _ElemByColors[c][e];
    }

    void MltCoMatGradDef(Index iEltGlob, Index iLoc, std::array<Real,Dim> &U) const
    {
        RealMatrix GradDef;
        RealMatrix CoMatGradDef;
        _FESpace.getGradDef(_Elements[iEltGlob],iLoc,GradDef);
        ArrayAlgebra::CoMat(GradDef,CoMatGradDef);
        ArrayAlgebra::MatMlt(CoMatGradDef,U);
    }

    void TransposeMltCoMatGradDef(Index iEltGlob, Index iLoc, std::array<Real,Dim> &U) const
    {
        RealMatrix GradDef;
        RealMatrix CoMatGradDef;
        _FESpace.getGradDef(_Elements[iEltGlob],iLoc,GradDef);
        ArrayAlgebra::CoMat(GradDef,CoMatGradDef);
        ArrayAlgebra::TransposeMatMlt(CoMatGradDef,U);
    }

    // ---------------------------------------------------------------------------------//

    const GaussLobattoElement & getFE() const { return _FESpace.getFE();}

    const FESpace & getFESpace() const {return _FESpace;}

    // ---------------------------------------------------------------------------------//

};


 
}//OndoMathX
 
