#pragma once

namespace OndoMathX {

 
 
template <typename GeoTransform> class CylinderBndy
{

private:

    GaussLobattoElement _GGLE_lateral_bndy;
    
    Circle<Identity<2>> _Circle;
    Cylinder<GeoTransform> _Cylinder;
    
    //Number of the boundary
    Index _iBoundary;

    //Number of Elements in the Lateral Boundary
    Index _Nlat;
    Index _NDoFlat;
    
    //Boolean indicating if the space should be seen as a trace space
    bool _isTraceSpace;
    
    
public:
    
    typedef GaussLobattoInterpolator Interpolator;
    
    void _loc_bdny2loc_cylinder(Index iElt_bndy,Index iloc_bndy, Index & iElt_cylinder,Index & iloc_cylinder) const
    {
        Index i;
        switch(_iBoundary)
        {
            case 0:
                iElt_cylinder = iElt_bndy;
                iloc_cylinder = iloc_bndy;
            break;
                
            case 1:
            {
                Index iElt_Z = iElt_bndy / _Circle._NTheta;
                Index iLocTheta = iloc_bndy % (_Circle._P + 1);
                Index iLocZ = iloc_bndy / (_Circle._P + 1);
                iElt_cylinder = iElt_bndy + (iElt_Z + 1) * (_Circle._NElem - _Circle._NTheta); 
                iloc_cylinder = _Circle._P * (_Circle._P + 1) + iLocZ * (_Circle._P + 1)*(_Circle._P + 1) + iLocTheta;     
                break;
            }
                
            case 2:
                iElt_cylinder = iElt_bndy + _Cylinder._NElem - _Circle._NElem;
                iloc_cylinder = iloc_bndy + (_Cylinder._NLocDoF - _Circle._NLocDoF);
            break;
                
            default: assert(false && "Should have reached return.");
        }
        
            
    }
    

   
    Index _glob_bdny2glob_cylinder(Index iDof) const
    { 
        switch(_iBoundary)
        {
            case 0: 
                return iDof; 
                break;
            case 1: 
            {
                Index iTheta = iDof % _Circle._NDoFTheta;
                Index iZ = iDof / _Circle._NDoFTheta;
                return iZ*_Circle._NDoF + (_Circle._NDoF - _Circle._NDoFTheta ) + iTheta ; 
                break;
            }    
            case 2:  
                return iDof + _Cylinder._NDoF - _Circle._NDoF  ; 
                break;           
        }
        assert(false && "Should have reached return.");
        return 0;
    }


    //This one is used to get vertex coordinate of the boundary, by using the getvertexcoord implemented with the cylinder 
    Index _glob_bdny2glob_cylinder_vertex(Index iVertex) const
    { 
        Index NVertCircle = _Circle.getNumVertices();
        switch(_iBoundary)
        {
            case 0: 
                return iVertex; 
                break;
            case 1:
            {
                Index iTheta = iVertex % _Circle._NTheta;
                Index iZ = iVertex / _Circle._NTheta;
                return iZ*NVertCircle + (NVertCircle - _Circle._NTheta ) + iTheta ; 
                //Because of periodicity, the number of vertex = number of elements in radial direction
                break;
            }    
            case 2:  
                return iVertex + _Cylinder.getNumVertices() - NVertCircle  ; 
                break;         
        }
        assert(false && "Should have reached return.");
        return 0;
    }
    
    


    static const Index Dim = 2;
   
    CylinderBndy(Index iBoundary, Cylinder<GeoTransform> aCylinder) : _Circle( Circle<Identity<2>>(1,1)), _Cylinder(aCylinder)
    {
        assert(iBoundary<3);
        
        _iBoundary = iBoundary;
        _Circle = Circle<Identity<2>>(_Cylinder._PX,_Cylinder._NX);
 
        _Nlat = _Circle._NTheta * _Cylinder._NZ;
        _NDoFlat = _Circle._NDoFTheta * _Cylinder._NDoFZ ;

        if (_iBoundary == 1)
            _GGLE_lateral_bndy = GaussLobattoElement(_Cylinder._PX+1,_Cylinder._PZ+1,0);
        
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
            return _Cylinder.getNumDoFs();
        else
        {
            if (_iBoundary == 1)
                return _NDoFlat; 
            else
                return _Circle.getNumDoFs();              
        }
        assert(false && "Should have reached return.");
        return 0; 
    }


    Index getNumLocDoFs() const 
    { 
        if (_iBoundary == 1)
            return (_Circle._P + 1) * (_Cylinder._PZ + 1) ;
        else
            return _Circle.getNumLocDoFs();
           
        assert(false && "Should have reached return.");
            return 0; 
    }
    
   
    Index getNumColors() const 
    {
        return _Circle.getNumColors();
    }
    
   
    Index getNumElements(Index iColor) const 
    {
        if (_iBoundary == 1)
            return _Circle._NTheta * _Cylinder._NZ / 4;
        else
            return _Circle.getNumElements(iColor);

        assert(false && "Should have reached return.");
            return 0; 
    }
    
   
    Index getNumElements() const 
    {
       if (_iBoundary == 1)
            return _Circle._NTheta * _Cylinder._NZ;
        else
            return _Circle.getNumElements();

        assert(false && "Should have reached return.");
            return 0; 
    }
    
 
    const GaussLobattoElement & getFE() const 
    { 
        if (_iBoundary == 1)
            return _GGLE_lateral_bndy ; 
        else
            return _Circle.getFE();

        // assert(false && "Should have reached return.");
        //     return 0; 
        
    }


    Index getNumVertices() const 
    {
        return 2 * _Circle.getNumVertices() + (_Cylinder._NZ - 2) * _Circle._NTheta;
    }
        
    
    void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
    {
        if (_iBoundary == 1)
        {
            numbering.resize(4);
            
            //this one treated like a square but considering the periodicity for the last 2 nodes
            // in the radial direction
            Index iEltGlobTheta = iEltGlob % _Circle._NTheta;
            Index iEltGlobZ = iEltGlob / _Circle._NTheta;
            
            numbering[0] = (iEltGlobZ + 0) * _Circle._NTheta + (iEltGlobTheta + 0);
            numbering[1] = (iEltGlobZ + 0) * _Circle._NTheta + (iEltGlobTheta + 1) % _Circle._NTheta;
            numbering[2] = (iEltGlobZ + 1) * _Circle._NTheta + (iEltGlobTheta + 1) % _Circle._NTheta;
            numbering[3] = (iEltGlobZ + 1) * _Circle._NTheta + (iEltGlobTheta + 0);
        }
            
        else
            _Circle.getElementVertices(iEltGlob,numbering);
    }
    
    
    Index EltColorNum2EltGlobalNum(Index c, Index e) const 
    {
        if (_iBoundary == 1)
        {
            
            // As NTheta is always even, we can treat it as a normal rectangle (no prob of periodicity)
            Index eTheta = e%(_Circle._NTheta/2);
            Index eZ = e/(_Circle._NTheta/2);
            Index geTheta{0},geZ{0};

            switch (c)
            {
                case 0: geTheta=2*eTheta;   geZ=2*eZ;   break;
                case 1: geTheta=2*eTheta+1; geZ=2*eZ;   break;
                case 2: geTheta=2*eTheta  ; geZ=2*eZ+1; break;
                case 3: geTheta=2*eTheta+1; geZ=2*eZ+1; break;
                default: assert(false);
            }
            return geZ*_Circle._NTheta+geTheta;
        }

        else
            return _Circle.EltColorNum2EltGlobalNum(c,e);
    }
    // ---------------------------------------------------------------------------------//

  
    void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
    {
        Index  iVertex_Cyl = _glob_bdny2glob_cylinder_vertex(iVertex);
        _Cylinder.getVertexCoordinate(iVertex_Cyl, xyz);  
    }

    Index Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        Index iTmp;
        if (_iBoundary == 1)
        {
            Index iEltGlobTheta = iEltGlob % _Circle._NTheta;
			Index iEltGlobZ = iEltGlob / _Circle._NTheta;

			Index iLocTheta = iLoc % (_Circle._P + 1);
			Index iLocZ = iLoc / (_Circle._P + 1);

			iTmp = (iEltGlobZ*_Cylinder._PZ + iLocZ) * _Circle._NDoFTheta + (iEltGlobTheta * _Circle._P + iLocTheta) % _Circle._NDoFTheta;
        }
        else
            iTmp = _Circle.Loc2Glob(iEltGlob,iLoc);

        if (!_isTraceSpace) 
            return iTmp;
        else 
            return _glob_bdny2glob_cylinder(iTmp);
    }
    
  /*  Index BndyElement2VolElement(Index iEltGlob) const
    {
        if (!_isTraceSpace) return iEltGlob;
        else
        {
            Index iLoc =0;
            Index iEltGlob_cylinder;
            Index iLoc_cylinder;
             
            _loc_bdny2loc_cylinder(iEltGlob,iLoc,iEltGlob_cylinder,iLoc_cylinder);
            
            return iEltGlob_cylinder;
        }
    }*/
    
    void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        Index NLocDoF = getNumLocDoFs();
        loc2glob.resize(NLocDoF);
        
        for (Index iLoc=0; iLoc < NLocDoF; iLoc++)
            loc2glob[iLoc] = Loc2Glob(iEltGlob,iLoc);
    }


    Index Loc2GlobDisc(Index iEltGlob, Index iLoc) const
    {
        Index NLocDoF = _Circle.getNumLocDoFs();
        
        if (!_isTraceSpace)
        {
            return iEltGlob*NLocDoF+iLoc;
        } 
        else 
        {
            Index ve;
            Index viloc;
            _loc_bdny2loc_cylinder(iEltGlob,iLoc,ve,viloc);

            //ve = getVolumeElement(iEltGlob);
            //viloc = Loc2Loc(iLoc);
            Index NLocDoF = _Cylinder.getNumLocDoFs();

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
        Index iEltGlob_cylinder;
        Index iLoc_cylinder;
        
        _loc_bdny2loc_cylinder(iEltGlob,iLoc,iEltGlob_cylinder,iLoc_cylinder);
        
        std::array<Real, 3> n0;
        
        switch(_iBoundary)
        {
            case 0:
                n0[0]= 0.0;n0[1]= 0.0;n0[2]=-1.0;break;
            case 1: 
            {
                // MIGHT BE WRONG:
                // Get the position of the Dof in the cylinder, without considering the third component
                //As it is radial, the norm should correspond to the vector position(?)
                const RealVector & P = _Cylinder._GLE.getPointCoordinate(iLoc_cylinder);
                _Cylinder._TemplateTransform(iEltGlob_cylinder, P, n0);
                n0[2]=0.0;
                break;
            }   
            case 2: 
           
                 n0[0]= 0.0;n0[1]= 0.0;n0[2]= 1.0;break;
            
               
            default: assert(false && "Should have reached break.");
        }
         
        //still, I am not sure of what it is done in this passage, does it work with what is done in case 1?
        _Cylinder.MltCoMatGradDef(iEltGlob_cylinder,iLoc_cylinder,n0);
        
        //Real norm = ArrayAlgebra::Norm(n0);
        
         n[0] = n0[0];// /norm;
         n[1] = n0[1];// /norm;
         n[2] = n0[2];// /norm;
    }
    
    
    
    // Jacobian of the transform from the reference line to the element iEltGlob at the quadrature position iLoc.
    Real getJacobian(Index iEltGlob, Index iLoc) const
    {
        Index iEltGlob_cylinder;
        Index iLoc_cylinder;
        
        _loc_bdny2loc_cylinder(iEltGlob,iLoc,iEltGlob_cylinder,iLoc_cylinder);
        
        std::array<Real, 3> n0;
        
        switch(_iBoundary)
        {
            case 0: n0[0]= 0.0;n0[1]= 0.0;n0[2]=-1.0;break;
            case 1: 
            {
                // same comment done on the norm
                const RealVector & P = _Cylinder._GLE.getPointCoordinate(iLoc_cylinder);
                _Cylinder.TemplateTransform(iEltGlob_cylinder, P, n0);
                n0[2]=0.0;
                break;
            }   
            case 2: n0[0]= 0.0;n0[1]= 0.0;n0[2]= 1.0;break;
            default : assert(false && "Should have reached break.");
        }
        
        _Cylinder.MltCoMatGradDef(iEltGlob_cylinder,iLoc_cylinder,n0);
        
        return ArrayAlgebra::Norm(n0);
    }
    
    
    
    void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
    {
        Index iLoc{0}, iEltGlob_cylinder{0},iLoc_cylinder{0} ;
        _loc_bdny2loc_cylinder(iEltGlob, iLoc, iEltGlob_cylinder,iLoc_cylinder );
     
    
        _Cylinder.getBarycenter(iEltGlob_cylinder,xyzCenter);
    
    }
    
    
    void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector& xyz) const
    {
        Index iElt_cylinder, iloc_cylinder;
        _loc_bdny2loc_cylinder(iEltGlob, iLoc, iElt_cylinder, iloc_cylinder);
    
        _Cylinder.getDoFCoordinate(iElt_cylinder,iloc_cylinder,xyz);
    
    }
    
    Mesh getMesh()
    {
        std::vector<std::shared_ptr<Point> > coords;
        std::vector<std::shared_ptr<Point> > local_coords(4);
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
 
