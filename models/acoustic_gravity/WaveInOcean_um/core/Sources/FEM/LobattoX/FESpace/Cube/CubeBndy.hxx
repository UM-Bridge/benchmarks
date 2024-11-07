#pragma once
 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
  
template <typename GeoTransform> class CubeBndy
{

private:

    Square<Identity<2>> _Square;
    Cube<GeoTransform> _Cube;
    
    //Number of the boundary
    Index _iBoundary;
    
    //Boolean indicating if the space should be seen as a trace space
    bool _isTraceSpace;
    
    
public:
    
    typedef GaussLobattoInterpolator Interpolator;
    
    void _loc_bdny2loc_cube(Index iElt_bndy,Index iloc_bndy, Index & iElt_cube,Index & iloc_cube) const
    {
        Index iElt_bndy_X = iElt_bndy%_Square.getNX();
        Index iElt_bndy_Y = iElt_bndy/_Square.getNX();
        
        Index iloc_bndy_X = iloc_bndy%(_Square.getPX()+1);
        Index iloc_bndy_Y = iloc_bndy/(_Square.getPX()+1);
        
        switch(_iBoundary)
        {
            case 0:
                iElt_cube = iElt_bndy;
                iloc_cube = iloc_bndy;
            break;
            case 1:
                iElt_cube = _Cube._NXY      *iElt_bndy_Y+iElt_bndy_X;
                iloc_cube = _Cube._NLocDoFXY*iloc_bndy_Y+iloc_bndy_X;
            break;
            case 2:
                iElt_cube = _Cube._NXY      *iElt_bndy_Y+(iElt_bndy_X+1)*_Cube._NX      -1;
                iloc_cube = _Cube._NLocDoFXY*iloc_bndy_Y+(iloc_bndy_X+1)*_Cube._NLocDoFX-1;
            break;
            case 3:
                iElt_cube = _Cube._NXY      *(iElt_bndy_Y+1)+iElt_bndy_X-_Cube._NX;
                iloc_cube = _Cube._NLocDoFXY*(iloc_bndy_Y+1)+iloc_bndy_X-_Cube._NLocDoFX;
            break;
            case 4:
                iElt_cube = _Cube._NXY      *(iElt_bndy_Y)+(iElt_bndy_X)*_Cube._NX;
                iloc_cube = _Cube._NLocDoFXY*(iloc_bndy_Y)+(iloc_bndy_X)*_Cube._NLocDoFX;
            break;
            case 5:
                iElt_cube = _Cube._NElem  -_Cube._NXY      +iElt_bndy;
                iloc_cube = _Cube._NLocDoF-_Cube._NLocDoFXY+iloc_bndy;
            break;
                
            default: {assert(false);}
        }
    }
    
    Index _glob_bdny2glob_cube(Index iDof) const
    {
        Index nDoFX = _Square.getNX()*_Square.getPX()+1;
        Index iDoFX = iDof%nDoFX;
        Index iDoFY = iDof/nDoFX;
        
        switch(_iBoundary)
        {
            case 0: return iDof; break;
            case 1: return _Cube._NDoFXY*iDoFY+iDoFX; break;
            case 2: return _Cube._NDoFXY*iDoFY+(iDoFX+1)*_Cube._NDoFX-1; break;
            case 3: return _Cube._NDoFXY*(iDoFY+1)+iDoFX-_Cube._NDoFX; break;
            case 4: return _Cube._NDoFXY*iDoFY+iDoFX*_Cube._NDoFX; break;
            case 5: return _Cube._NDoF-_Cube._NDoFXY + iDof; break;
            default: {assert(false);}
        }
       
        return 0;
    }
    
 

    static const Index Dim = 2;
   
    CubeBndy(Index iBoundary, Cube<GeoTransform> aCube) : _Square( Square<Identity<2>>(1,1,1,1)), _Cube(aCube)
    {
        _iBoundary = iBoundary;
 
        if (_iBoundary == 0 || _iBoundary == 5)
            _Square = Square<Identity<2>>(_Cube._PX,_Cube._PY,_Cube._NX,_Cube._NY);
        else if (_iBoundary == 1 || _iBoundary == 3)
            _Square = Square<Identity<2>>(_Cube._PX,_Cube._PZ,_Cube._NX,_Cube._NZ);
        else if (_iBoundary == 2 || _iBoundary == 4)
            _Square = Square<Identity<2>>(_Cube._PY,_Cube._PZ,_Cube._NY,_Cube._NZ);
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
            return _Cube.getNumDoFs();
        else
            return _Square.getNumDoFs();
    }
    
    Index getNumLocDoFs() const { return _Square.getNumLocDoFs();}
    Index getNumColors() const {return _Square.getNumColors();}
    Index getNumElements(Index iColor) const {return _Square.getNumElements(iColor);}
    Index getNumElements() const {return _Square.getNumElements();}
    const GaussLobattoElement & getFE() const { return _Square.getFE();}
    Index getNumVertices() const {return _Square.getNumVertices();}
    void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
    {_Square.getElementVertices(iEltGlob,numbering);}
    Index EltColorNum2EltGlobalNum(Index c, Index e) const {return _Square.EltColorNum2EltGlobalNum(c,e);}
    // ---------------------------------------------------------------------------------//

    void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
    {
        RealVector uv;
        RealVector uvw;
    
        _Square.getVertexCoordinate(iVertex,uv);
        
        switch(_iBoundary)
        {
            case 0: uvw[0]=uv[0]; uvw[1]=uv[1]; uvw[2]=  0.0; break;
            case 1: uvw[0]=uv[0]; uvw[1]=  0.0; uvw[2]=uv[1]; break;
            case 2: uvw[0]=  1.0; uvw[1]=uv[0]; uvw[2]=uv[1]; break;
            case 3: uvw[0]=uv[0]; uvw[1]=  1.0; uvw[2]=uv[1]; break;
            case 4: uvw[0]=  0.0; uvw[1]=uv[0]; uvw[2]=uv[1]; break;
            case 5: uvw[0]=uv[0]; uvw[1]=uv[1]; uvw[2]=  1.0; break;
            default: {assert(false);}
        }
        
        // Applying transformation.
        _Cube._GeoTransform.Eval(uvw, xyz);
    }

    Index Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        Index iTmp = _Square.Loc2Glob(iEltGlob,iLoc);
        
        if (!_isTraceSpace) return iTmp;
        else return _glob_bdny2glob_cube(iTmp);
    }
    
    Index BndyElement2VolElement(Index iEltGlob) const
    {
        if (!_isTraceSpace) return iEltGlob;
        else
        {
            Index iLoc =0;
            Index iEltGlob_cube;
            Index iLoc_cube;
             
            _loc_bdny2loc_cube(iEltGlob,iLoc,iEltGlob_cube,iLoc_cube);
            
            return iEltGlob_cube;
        }   
    }
    
    void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        Index NLocDoF = _Square.getNumLocDoFs();
        loc2glob.resize(NLocDoF);
        
        for (Index iLoc=0; iLoc < NLocDoF; iLoc++)
            loc2glob[iLoc] = Loc2Glob(iEltGlob,iLoc);
    }

/*
    Index getVolumeElement(Index iElt) const
    {
        switch(_iBoundary)
        { 
            case 0: return iElt; break;
            case 5: return _Cube._NElem -_Square._NElem  + iElt; break;

            default: {assert(false);}
        }
        
        return 0;
    }


    Index Loc2Loc(Index iLoc) const
    {
        switch(_iBoundary)
        {
            case 0: 
                return iLoc; 
                break;

            default: {assert(false);}
        }
        return 0;
    }*/


    Index Loc2GlobDisc(Index iEltGlob, Index iLoc) const
    {
        Index NLocDoF = _Square.getNumLocDoFs();
        
        if (!_isTraceSpace)
        {
            return iEltGlob*NLocDoF+iLoc;
        } 
        else 
        {
            Index ve;
            Index viloc;
            _loc_bdny2loc_cube(iEltGlob,iLoc,ve,viloc);


            //ve = getVolumeElement(iEltGlob);
            //viloc = Loc2Loc(iLoc);
            Index NLocDoF = _Cube.getNumLocDoFs();

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
    void getNormal(Index iEltGlob, Index iLoc,  std::array<Real, 3> &n) const
    {
        Index iEltGlob_cube;
        Index iLoc_cube;
        
        _loc_bdny2loc_cube(iEltGlob,iLoc,iEltGlob_cube,iLoc_cube);
        
        std::array<Real, 3> n0;
        
        switch(_iBoundary)
        {
            case 0: n0[0]= 0.0;n0[1]= 0.0;n0[2]=-1.0;break;
            case 1: n0[0]= 0.0;n0[1]=-1.0;n0[2]= 0.0;break;
            case 2: n0[0]= 1.0;n0[1]= 0.0;n0[2]= 0.0;break;
            case 3: n0[0]= 0.0;n0[1]= 1.0;n0[2]= 0.0;break;
            case 4: n0[0]=-1.0;n0[1]= 0.0;n0[2]= 0.0;break;
            case 5: n0[0]= 0.0;n0[1]= 0.0;n0[2]= 1.0;break;
            default: {assert(false);}
        }
        
        _Cube.MltCoMatGradDef(iEltGlob_cube,iLoc_cube,n0);
        
        Real norm = ArrayAlgebra::Norm(n0);
        
         n[0] = n0[0]/norm;
         n[1] = n0[1]/norm;
         n[2] = n0[2]/norm;
    }
    
    
    
    // Jacobian of the transform from the reference line to the element iEltGlob at the quadrature position iLoc.
    Real getJacobian(Index iEltGlob, Index iLoc) const
    {
        Index iEltGlob_cube;
        Index iLoc_cube;
        
        _loc_bdny2loc_cube(iEltGlob,iLoc,iEltGlob_cube,iLoc_cube);
        
        std::array<Real, 3> n0;
        
        switch(_iBoundary)
        {
            case 0: n0[0]= 0.0;n0[1]= 0.0;n0[2]=-1.0;break;
            case 1: n0[0]= 0.0;n0[1]=-1.0;n0[2]= 0.0;break;
            case 2: n0[0]= 1.0;n0[1]= 0.0;n0[2]= 0.0;break;
            case 3: n0[0]= 0.0;n0[1]= 1.0;n0[2]= 0.0;break;
            case 4: n0[0]=-1.0;n0[1]= 0.0;n0[2]= 0.0;break;
            case 5: n0[0]= 0.0;n0[1]= 0.0;n0[2]= 1.0;break;
            default: {assert(false);}
        }
        
        _Cube.MltCoMatGradDef(iEltGlob_cube,iLoc_cube,n0);
        
        return ArrayAlgebra::Norm(n0);
    }
    
    
    
    void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
    {
        RealVector uv;
        RealVector uvw;
    
        _Square.getBarycenter(iEltGlob,uv);
    
        switch(_iBoundary)
        {
            case 0: uvw[0]=uv[0]; uvw[1]=uv[1]; uvw[2]=  0.0; break;
            case 1: uvw[0]=uv[0]; uvw[1]=  0.0; uvw[2]=uv[1]; break;
            case 2: uvw[0]=  1.0; uvw[1]=uv[0]; uvw[2]=uv[1]; break;
            case 3: uvw[0]=uv[0]; uvw[1]=  1.0; uvw[2]=uv[1]; break;
            case 4: uvw[0]=  0.0; uvw[1]=uv[0]; uvw[2]=uv[1]; break;
            case 5: uvw[0]=uv[0]; uvw[1]=uv[1]; uvw[2]=  1.0; break;
            default: {assert(false);}
        }
        
        // Applying transformation.
        _Cube._GeoTransform.Eval(uvw, xyzCenter);
    }
    
    
    void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector& xyz) const
    {
        RealVector uv;
        RealVector uvw;
    
        _Square.getDoFCoordinate(iEltGlob,iLoc,uv);
    
        switch(_iBoundary)
        {
            case 0: uvw[0]=uv[0]; uvw[1]=uv[1]; uvw[2]=  0.0; break;
            case 1: uvw[0]=uv[0]; uvw[1]=  0.0; uvw[2]=uv[1]; break;
            case 2: uvw[0]=  1.0; uvw[1]=uv[0]; uvw[2]=uv[1]; break;
            case 3: uvw[0]=uv[0]; uvw[1]=  1.0; uvw[2]=uv[1]; break;
            case 4: uvw[0]=  0.0; uvw[1]=uv[0]; uvw[2]=uv[1]; break;
            case 5: uvw[0]=uv[0]; uvw[1]=uv[1]; uvw[2]=  1.0; break;
            default: {assert(false);}
        }
        
        // Applying transformation.
        _Cube._GeoTransform.Eval(uvw, xyz);
    }
    
    
    
    
    //void MltCoMatGradDef(Index iEltGlob, Index iLoc,
   //      std::array<Real, 2>& LocalSolutionGradient) const
   // {
   //
   // }
    
    //void TransposeMltCoMatGradDef(Index iEltGlob, Index iLoc,
  //       std::array<Real, 2>& LocalSolutionGradient) const
  //  {
  //
  //  }
    
    
    Mesh getMesh()
    {
        std::vector<std::shared_ptr<Point> > coords;
        std::vector<std::shared_ptr<const Point> > local_coords(4);
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
            local_coords[2]=coords[numbering[2]];
            local_coords[3]=coords[numbering[3]];
            
            elements.push_back(std::make_shared<Quadrangle4>(local_coords,0));
        }
        
        return Mesh(coords,elements);
    }

    Index getLabel(Index iEltGlob) const 
    {
        return 0;
    }

    bool isAffine(const Index iEltGlob) const
    {
        return false;
    }
};

 

 
} // OndoMathX


// -----------------------------------------------------------------------------------------//
 
