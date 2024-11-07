#pragma once
 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {


template <typename GeoTransform> class SquareBndy
{

private:

    Line<Identity<1>> _Line;
    const Square<GeoTransform> & _Square;
    
    //Number of the boundary
    Index _iBoundary;
    
    //Boolean indicating if the space should be seen as a trace space
    bool _isTraceSpace;
    
    std::vector<Index> _Labels;

public:

    typedef GaussLobattoInterpolator Interpolator;
    
    void _loc_bdny2loc_sq(Index iElt_bndy,Index iloc_bndy, Index & iElt_sq,Index & iloc_sq) const
    {
        switch(_iBoundary)
        {
            case 0:
                iElt_sq = iElt_bndy;
                iloc_sq = iloc_bndy;
            break;
            case 1:
                iElt_sq = _Square._NX      *(iElt_bndy+1)-1;
                iloc_sq = _Square._NLocDoFX*(iloc_bndy+1)-1;
            break;
            case 2:
                iElt_sq   = _Square._NElem-iElt_bndy-1;
                iloc_sq = _Square._NLocDoF-iloc_bndy-1;
            break;
            case 3:
                iElt_sq   = _Square._NElem -  _Square._NX*(iElt_bndy+1);
                iloc_sq   = _Square._NLocDoF- _Square._NLocDoFX*(iloc_bndy+1);
            break;
                
            default: {assert(false);}
        }
    }
    
    Index _glob_bdny2glob_sq(Index iDof) const
    {
        switch(_iBoundary)
        {
            case 0: return iDof; break;
            case 1: return _Square._NDoFX*(iDof+1)-1; break;
            case 2: return _Square._NDoF-iDof-1; break;
            case 3: return _Square._NDoF-_Square._NDoFX*(iDof+1); break;
            default: {assert(false);}
        }
        
        return 0;
    }

  

    /*! \brief Definition of the dimension as a static variable.*/
    static const Index Dim = 1;
   
    SquareBndy(Index iBoundary, const Square<GeoTransform> & aSquare) : _Line( Line<Identity<1>>(1,1)), _Square(aSquare) 
    {
        _iBoundary = iBoundary;
 
        if (_iBoundary == 0 || _iBoundary == 2)
            _Line = Line<Identity<1>>(_Square._PX,_Square._NX);
        else if (_iBoundary == 1 || _iBoundary == 3)
            _Line = Line<Identity<1>>(_Square._PY,_Square._NY);
        else {assert(false);}

        //at initialisation the space is seen as a trace space
        _isTraceSpace =true;

    }
    
    void setAsTraceSpace()
    {
        _isTraceSpace = true;
    }
    
    void setAsFESpace()
    {
        _isTraceSpace = false;
    }
    

    // ---------------------------------------------------------------------------------//
    Index getNumDoFs() const
    {
        if(_isTraceSpace)
            return _Square.getNumDoFs();
        else
            return _Line.getNumDoFs();
    }

    Index getNumDoFsDisc() const
    {
        if(_isTraceSpace)
            return _Square.getNumDoFsDisc();
        else
            return getNumElements()*getNumLocDoFs();
    }




    
    Index getNumLocDoFs() const { return _Line.getNumLocDoFs();}
    Index getNumColors() const {return _Line.getNumColors();}
    Index getNumElements(Index iColor) const {return _Line.getNumElements(iColor);}
    Index getNumElements() const {return _Line.getNumElements();}
    const GaussLobattoElement & getFE() const { return _Line.getFE();}
    Index getNumVertices() const {return _Line.getNumVertices();}
    void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
    {_Line.getElementVertices(iEltGlob,numbering);}
    Index EltColorNum2EltGlobalNum(Index c, Index e) const {return _Line.EltColorNum2EltGlobalNum(c,e);}

    const Square<GeoTransform> & getFESpace() const {return _Square;}
    // ---------------------------------------------------------------------------------//

    void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
    {
        RealVector uvw;
        
        uvw[2]=0.0;
        
        switch(_iBoundary)
        {
            case 0: uvw[0]=iVertex*_Square._hX; uvw[1]=0.0; break;
            case 1: uvw[0]=1.0; uvw[1]=iVertex*_Square._hY; break;
            case 2: uvw[0]=1.0-iVertex*_Square._hX; uvw[1]=1.0; break;
            case 3: uvw[0]=0.0; uvw[1]=1.0-iVertex*_Square._hY; break;
        }
        
        // Applying transformation.
        _Square._GeoTransform.Eval(uvw, xyz);
    }

    Index Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        Index iTmp = _Line.Loc2Glob(iEltGlob,iLoc);
        
        if (!_isTraceSpace) return iTmp;
        else return _glob_bdny2glob_sq(iTmp);
    }

    void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        Index NLocDoF = _Line.getNumLocDoFs();
        loc2glob.resize(NLocDoF);
        
        for (Index iLoc=0; iLoc < NLocDoF; iLoc++)
            loc2glob[iLoc] = Loc2Glob(iEltGlob,iLoc);
    }

        // Tested : ok 
    Index Loc2Loc(Index iLoc) const
    {
        switch(_iBoundary)
        {
            case 0: 
                return iLoc; 
                break;
            case 1: 
                return (_Square._PX+1)*iLoc+_Square._PX; 
                break;
            case 2: 
                return (_Square._PX+1)*(_Square._PY+1)-1-iLoc; 
                break;
            case 3:  
                return (_Square._PX+1)*(_Square._PY+1)-1-_Square._PX-(_Square._PX+1)*iLoc; 
                break;
            default: {assert(false);}
        }
        return 0;
    }

        // Tested : ok
    Index getVolumeElement(Index iElt) const
    {
        switch(_iBoundary)
        { 
            case 0: return iElt; break;
            case 1: return _Square._NX*(iElt+1)-1; break;
            case 2: return _Square._NElem-iElt-1; break;
            case 3: return _Square._NElem-_Square._NX*(iElt+1); break;
            default: {assert(false);}
        }
        
        return 0;
    }


    Index Loc2GlobDisc(Index iEltGlob, Index iLoc) const
    {
        Index NLocDoF = _Line.getNumLocDoFs();
        
        if (!_isTraceSpace)
        {
            return iEltGlob*NLocDoF+iLoc;
        } 
        else 
        {
            Index ve = getVolumeElement(iEltGlob);
            Index viloc = Loc2Loc(iLoc);
            Index NLocDoF = _Square.getNumLocDoFs();

            return  NLocDoF*ve + viloc;
        }
    }

    void Loc2GlobDisc(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        loc2glob.resize(getNumLocDoFs());

        for (Index iLoc=0; iLoc < getNumLocDoFs(); iLoc++)
            loc2glob[iLoc] = Loc2GlobDisc(iEltGlob,iLoc);
            
    }

    
    
    // Jacobian of the transform from the reference line to the element iEltGlob at the quadrature position iLoc.
    void getNormal(Index iEltGlob, Index iLoc,  std::array<Real, 2> &n) const
    {
        Index iEltGlob_sq;
        Index iLoc_sq;
        
        _loc_bdny2loc_sq(iEltGlob,iLoc,iEltGlob_sq,iLoc_sq);
        
        std::array<Real, 2> n0;
        
        switch(_iBoundary)
        {
            case 0: n0[0]= 0.0;n0[1]=-1.0; break;
            case 1: n0[0]= 1.0;n0[1]= 0.0; break;
            case 2: n0[0]= 0.0;n0[1]= 1.0; break;
            case 3: n0[0]=-1.0;n0[1]= 0.0;break;
            default: {assert(false);}
        }

        RealMatrix2x2 DF;
        RealMatrix2x2 ComatDF;
        
        _Square.getGradDef(iEltGlob_sq,iLoc_sq,DF);
        ArrayAlgebra::CoMat(DF,ComatDF);
        ArrayAlgebra::MatMlt(ComatDF,n0);

         n[0] = n0[0];
         n[1] = n0[1];
    }
    
    void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
    {
        RealVector uvw;
     
        uvw[0] = 0.0;
        uvw[1] = 0.0;
        uvw[2] = 0.0;
        
        switch(_iBoundary)
        {
            case 0:
                uvw[0]=(iEltGlob+0.5)*_Square._hX;
            break;
            case 1:
                uvw[0]=1.0;
                uvw[1]=(iEltGlob+0.5)*_Square._hY;
            break;
            case 2:
                uvw[0]=1-(iEltGlob+0.5)*_Square._hX;
                uvw[1]=1.0;
            break;
            case 3:
                uvw[1]=1.0-(iEltGlob+0.5)*_Square._hY;
            break;
            default: {assert(false);}
        }

        _Square._GeoTransform.Eval(uvw, xyzCenter);
    }
    
    
    void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector& xyz) const
    {
        const RealVector& p = _Line.getFE().getPointCoordinate(iLoc);
    
        RealVector uvw;
     
        uvw[0] = 0.0;
        uvw[1] = 0.0;
        uvw[2] = 0.0;
        
        switch(_iBoundary)
        {
            case 0:
                uvw[0]=(iEltGlob+p[0])*_Square._hX;
            break;
            case 1:
                uvw[0]=1.0;
                uvw[1]=(iEltGlob+p[0])*_Square._hY;
            break;
            case 2:
                uvw[0]=1.0-(iEltGlob+p[0])*_Square._hX;
                uvw[1]=1.0;
            break;
            case 3:
                uvw[1]=1.0-(iEltGlob+p[0])*_Square._hY;
            break;
            default: {assert(false);}
        }

        _Square._GeoTransform.Eval(uvw, xyz);
    }
    
    
    
    Mesh getMesh()
    {
        std::vector<std::shared_ptr<Point> > coords;
        std::vector<std::shared_ptr<const Point> > local_coords(2);
        std::vector<std::shared_ptr<Element> > elements;
        
        RealVector xyz;
        for (Index v = 0; v<getNumVertices(); v++)
        {
            getVertexCoordinate(v,xyz);
            coords.push_back(std::make_shared<Point>(xyz[0],xyz[1],xyz[2],v));
        }
        
        for (Index e = 0; e<getNumElements(); e++)
        {
            std::vector<Index> numbering;
            getElementVertices(e,numbering);
                
            local_coords[0]=coords[numbering[0]];
            local_coords[1]=coords[numbering[1]];
            
            elements.push_back(std::make_shared<Segment2>(local_coords,0));
        }
        
        return Mesh(coords,elements);
    }



    void setLabels(std::function<Index(RealVector)> label_from_barycenter_pos)
    {
        _Labels.resize(getNumElements());

        for (Index e = 0; e<getNumElements(); e++)
        {
            RealVector xyz; 
            getBarycenter(e,xyz);

            _Labels[e] = label_from_barycenter_pos(xyz);
        }     
    }

    Index getLabel(Index iEltGlob) const 
    {
        if (_Labels.size()==0) return 0;
        else return _Labels[iEltGlob];
    }

    bool isAffine(const Index iEltGlob) const
    {
        return false;
    }
    



};


 
}//OndoMathX
 
