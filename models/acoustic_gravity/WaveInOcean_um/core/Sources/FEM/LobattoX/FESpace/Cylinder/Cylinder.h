#pragma once

namespace OndoMathX {

 

    template <typename GeoTransform>  class CylinderBndy;
  
    template <
        typename GeoTransform = Identity<3>
    >
    class Cylinder : public FESpaceT<Cylinder<GeoTransform>,3> {
 
    private:

        Circle<Identity<2>> _Circle;
        
        // Associated spectral finite element
        GaussLobattoElement _GLE;
     
        // Number of elements by color (also by directions)
        Index _NElemColor_0,_NElemColor_1;

        //Order
        Index _PX;
        Index _PZ;
        
        // Number of elements
        Index _NX;
        Index _NZ;
    

        // Total number of elements
        Index _NElem;
 
        // Number of local DoF
        Index _NLocDoF;
       
        // Number of DoF
        Index _NDoFZ;
        Index _NDoF;
        
        // Storing mesh steps
        Real _hZ ;

        // Associated transformation of the template finite element space
        GeoTransform _GeoTransform;


                    // ---------------------------------------------------------------------------------//
        //Class methods


        void _TemplateTransform (const Index &iEltGlob, const RealVector &p,
                                        RealVector &uvw) const
        {
            Index iEltGlobZ = iEltGlob / _Circle._NElem;
            Index iEltGlobSurf = iEltGlob - iEltGlobZ * _Circle._NElem;
            
            _Circle._TemplateTransform(iEltGlobSurf,p,uvw);
            
            uvw[2] = (p[2] + iEltGlobZ) * _hZ;
        }


        void _TemplateTransformGradient(const Index &iEltGlob, const RealVector &p,
                                        RealVector &uvw, RealMatrix3x3 &F ) const
        {
            Index iEltGlobZ = iEltGlob / _Circle._NElem;
            Index iEltGlobSurf = iEltGlob - iEltGlobZ * _Circle._NElem;
            
            //2D tansform of the circle
            RealMatrix2x2 F2D;
            _Circle._TemplateTransformGradient(iEltGlobSurf,p,uvw,F2D);
            
            //Correction to take into account 3D aspects
            uvw[2] = (p[2] + iEltGlobZ) * _hZ;
            
            F[0][0] = F2D[0][0]; F[0][1] = F2D[0][1]; F[0][2] = 0;
            F[1][0] = F2D[1][0]; F[1][1] = F2D[1][1]; F[1][2] = 0;
            F[2][0] = 0        ; F[2][1] = 0        ; F[2][2] = _hZ;
        }
        
    public:

        
        // ---------------------------------------------------------------------------------//
        /*! \brief Definition of the dimension as a static variable.*/
        static const Index Dim = 3;
        // ---------------------------------------------------------------------------------//

        
    
        // ---------------------------------------------------------------------------------//
        Cylinder(Index PX, Index PZ,Index NX, Index NZ, GeoTransform aGeoTransform = IdentityMap3)
            : _Circle(PX, NX),_GeoTransform(aGeoTransform)
        {   
            _PX = PX;
            _PZ = PZ;
            _NX = NX;
            _NZ = NZ;

        
            // Number of degrees of freedom (DoF) int the first direction and in total.
            _NDoFZ = (_NZ * PZ + 1);
            
            // Computing space step.
            _hZ = 1.0 / _NZ;

            // ---------------------------------------------------------------------------------//
            _NElem   = _Circle._NElem   *  _NZ;
            _NDoF    = _Circle._NDoF    * _NDoFZ;
            _NLocDoF = _Circle._NLocDoF * (_PZ + 1);
    
        
            _NElemColor_0 = (_NZ / 2 + _NZ % 2);
            _NElemColor_1 = (_NZ / 2          );
            
            // Initializing spectral finite element.
            _GLE = GaussLobattoElement(PX+1,PX+1,PZ+1);
        }
        // ---------------------------------------------------------------------------------//



        // ---------------------------------------------------------------------------------//
        // Extracting number of color in the template geometry, used for parallelization.
        Index getNumColors() const
        {
            return 8;
        }
        
        //  Extracting number of element associated to a given color.
        Index getNumElements(Index iColor) const
        {
            switch (iColor)
            {
                case 0: return _NElemColor_0*_Circle.getNumElements(iColor);
                case 1: return _NElemColor_0*_Circle.getNumElements(iColor);
                case 2: return _NElemColor_0*_Circle.getNumElements(iColor);
                case 3: return _NElemColor_0*_Circle.getNumElements(iColor);
                case 4: return _NElemColor_1*_Circle.getNumElements(iColor-4);
                case 5: return _NElemColor_1*_Circle.getNumElements(iColor-4);
                case 6: return _NElemColor_1*_Circle.getNumElements(iColor-4);
                case 7: return _NElemColor_1*_Circle.getNumElements(iColor-4);
            }
            
            {assert(false);}
            
            return _NElem;
        }
                
        //  Extracting number of global number of element.
        Index getNumElements() const
        {
            return _NElem;
        }
        
        const GaussLobattoElement & getFE() const
        {
            return _GLE;
        }

        Index getNumVertices() const
        {
            return _Circle.getNumVertices() * (_NZ+1);
        }
        
        Index getNumVerticesQ2() const
        {
            Index NVertQ2 = 2;
            return _Circle.getNumVerticesQ2() * (_NZ*NVertQ2 +1);
        }

         
        void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
        {
            RealVector uvw;

            Index NVertSurface = _Circle.getNumVertices();
            
            Index iVertexZ = iVertex / NVertSurface;
            Index iVertexRel = iVertex - iVertexZ * NVertSurface;

            _Circle.getVertexCoordinate(iVertexRel, uvw);

            uvw[2]=iVertexZ*_hZ;
            
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        void getVertexCoordinateQ2(const Index iVertex, RealVector& xyz) const
        {
            RealVector uvw;

            Index NVertSurfaceQ2 = _Circle.getNumVerticesQ2();
            
            Index iVertexZQ2 = iVertex / NVertSurfaceQ2;
            Index iVertexRel = iVertex - iVertexZQ2 * NVertSurfaceQ2;

            _Circle.getVertexCoordinateQ2(iVertexRel, uvw);

            Index NVertQ2 = 2;
            Real  hZQ2 = _hZ/NVertQ2;
            uvw[2]=iVertexZQ2*hZQ2;
            
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

         
        void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            numbering.resize(8);
         
            Index NVertSurface = _Circle.getNumVertices();
            
            Index iEltGlobZ = iEltGlob / _Circle._NElem;
            Index iEltGlobSurf = iEltGlob - iEltGlobZ * _Circle._NElem;
            

            std::vector<Index> numberingSurf ;

            _Circle.getElementVertices(iEltGlobSurf, numberingSurf);

            for (Index i = 0; i < 4; i++)
            {
                numbering[i]   = numberingSurf[i] + (iEltGlobZ + 0) * NVertSurface ;
                numbering[i+4] = numberingSurf[i] + (iEltGlobZ + 1) * NVertSurface ;
            }
        }

        void getElementVerticesQ2(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            numbering.resize(27);

            
         
            Index NVertSurfaceQ2 = _Circle.getNumVerticesQ2();
            
            Index iEltGlobZ = iEltGlob / _Circle._NElem;
            Index iEltGlobSurf = iEltGlob - iEltGlobZ * _Circle._NElem;
            


            std::vector<Index> numberingSurf ;

            _Circle.getElementVerticesQ2(iEltGlobSurf, numberingSurf);

            // std::array<Index,27> temp_numberingQ2;
            // for (Index i = 0; i < 9; i++)
            // {
            //     temp_numberingQ2[i]     = numberingSurf[i] + (iEltGlobZ*2 + 0) * NVertSurfaceQ2 ;
            //     temp_numberingQ2[i+9]   = numberingSurf[i] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2 ;
            //     temp_numberingQ2[i+9*2] = numberingSurf[i] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            // }

            numbering[0]  = numberingSurf[0]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[8]  = numberingSurf[4]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[1]  = numberingSurf[1]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[9]  = numberingSurf[7]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[20] = numberingSurf[8]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[11] = numberingSurf[5]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[3]  = numberingSurf[3]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[13] = numberingSurf[6]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[2]  = numberingSurf[2]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;

            numbering[10] = numberingSurf[0] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[21] = numberingSurf[4] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[12] = numberingSurf[1] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[22] = numberingSurf[7] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[26] = numberingSurf[8] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[23] = numberingSurf[5] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[15] = numberingSurf[3] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[24] = numberingSurf[6] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[14] = numberingSurf[2] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;

            numbering[4]  = numberingSurf[0] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[16] = numberingSurf[4] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[5]  = numberingSurf[1] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[17] = numberingSurf[7] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[25] = numberingSurf[8] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[18] = numberingSurf[5] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[7]  = numberingSurf[3] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[19] = numberingSurf[6] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[6]  = numberingSurf[2] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;


            
            
        }
        // ---------------------------------------------------------------------------------//
 

        // ---------------------------------------------------------------------------------//
        // Extracting the total number of DoF.
		Index getNumDoFs() const
		{
            return _NDoF;
		}
        
        // Extracting the number of local DoF.
        Index getNumLocDoFs() const
        {
            return _NLocDoF;
        }

        bool isAffine(const Index iEltGlob) const
        {
            return false;
        }
        // ---------------------------------------------------------------------------------//


        // Fonction that maps the numbering of the element by color to a global numbering.
        // c is the index of the color.
        // e is the numbering of the element in the color e.
        // It returns the global numbering of the element.
        Index EltColorNum2EltGlobalNum(Index c, Index e) const
        {
            Index NElemColorSurf = _Circle._NElemColor;

            Index ez = e / NElemColorSurf;
            Index exy = e - ez*NElemColorSurf;
    
            if (c < 4)
            {
                Index gez=2*(ez%(_NZ/2+_NZ%2));
        
                return _Circle.EltColorNum2EltGlobalNum(c,exy) + gez*_Circle._NElem;
            }
            else
            {
                Index gez=2*(ez%(_NZ/2      ))+1;
        
                return _Circle.EltColorNum2EltGlobalNum(c-4,exy) + gez*_Circle._NElem;
            }
            return c;
        }
        

        // Extracting global DoF index from local DoF index.
        // iEltGlob is the index of the element in the color group.
        // iLoc is the local DoF index.
        // It returns the global index.
        Index _Loc2Glob(Index iEltGlob, Index iLoc) const
        {
            Index iEltGlobZ = iEltGlob / _Circle._NElem;
            Index iEltGlobSurf = iEltGlob - iEltGlobZ * _Circle._NElem;
                
            Index iLocZ = iLoc / _Circle._NLocDoF;
            Index iLocRel = iLoc - iLocZ * _Circle._NLocDoF;

            return _Circle._Loc2Glob( iEltGlobSurf, iLocRel ) + (iEltGlobZ * _PZ + iLocZ) * _Circle._NDoF;
        }
        // ---------------------------------------------------------------------------------//


    
        // ---------------------------------------------------------------------------------//
        /*! \brief Gives the barycenter of the element
            \param iEltGlob is a global index of an element.
            \param xyzCenter are the coordinates of the barycenter of the element
        */
        void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;
            
            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector& P = {0.5,0.5,0.5};

            _TemplateTransform(iEltGlob, P, uvw);

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyzCenter);
        }


        // Get the label of the element (here 0)
        Index getLabel(const Index iEltGlob) const
        {
            return 0;
        }
        

        /*! \brief Gives a quadrature point associated to a local DoF in an element.
            \param iEltGlob is a global index of an element.
            \param iLoc is a local index in the element iEltGlob.
            \param xyz are the coordinates of the quadrature point associated to iLoc which will be filled.
            */
        void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector& xyz) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector & P = _GLE.getPointCoordinate(iLoc);

            _TemplateTransform(iEltGlob, P, uvw);

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }


        void getDoFCoordinateOnCylinder(const Index iEltGlob, const Index iLoc, RealVector& uvw) const
        {
            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector & P = _GLE.getPointCoordinate(iLoc);

            _TemplateTransform(iEltGlob, P, uvw);
        }



        /*! \brief Gives a degree of freedom coordinate associated to a global DoF.
            \param iGlob is a global index of an DoF.
            \param xyz are the coordinates of the degree of freedo massociated to iLoc which will be filled.
        */
        void getDoFCoordinate(const Index iGlob, RealVector& xyz) const
        {
            RealVector uvw;
            
            // Storing local index of DoF and Global index of Element
            Index iLoc, iEltGlob;
  
            Index iGlobZ = iGlob / _Circle._NDoF;
            Index iGlobSurface = iGlob - iGlobZ * _Circle._NDoF;
            
            _Circle.getDoFCoordinate(iGlobSurface,uvw);
            
            Index iEltGlobZ = iGlobZ / _PZ;
            Index iLocZ = iGlobZ - iEltGlobZ * _PZ;
    
            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector & P = _GLE.getPointCoordinate(iLocZ*_Circle._NLocDoF);

            // Adjusting quadrature point coordinate from the element index.
            uvw[2] = (P[2] + iEltGlobZ) * _hZ;
            
            
            _GeoTransform.Eval(uvw, xyz);
        }
        
        
        void _getGradDef(Index iEltGlob, Index iLoc,
			std::array<std::array<Real, 3>,3>& GradDef) const
        {
            assert(GeoTransform::isContinuouslyDifferentiable);

            RealVector uvw;
            RealMatrix3x3 GradDef_Tmp;
            RealMatrix3x3 GeoGradDef;
  
            const RealVector& p = _GLE.getPointCoordinate(iLoc);
    
            _TemplateTransformGradient(iEltGlob, p, uvw, GradDef_Tmp);
            _GeoTransform.getGradDef(uvw,GeoGradDef);

            ArrayAlgebra::MatMlt(GeoGradDef,GradDef_Tmp,GradDef);       
      
        }

        void getGradDefGeoTransform(Index iEltGlob, Index iLoc,
			std::array<std::array<Real, 3>,3>& GeoGradDef) const
        {
            assert(GeoTransform::isContinuouslyDifferentiable);

            RealVector uvw;

  
            const RealVector& p = _GLE.getPointCoordinate(iLoc);
    
            _TemplateTransform(iEltGlob, p, uvw);
            _GeoTransform.getGradDef(uvw,GeoGradDef);
        }
       
        // ---------------------------------------------------------------------------------//
        
        CylinderBndy<GeoTransform> getBndyFESpace(Index iBoundary)
        {
            assert(iBoundary<3);
            
            return CylinderBndy<GeoTransform>(iBoundary,*this);
        }
        
        friend class CylinderBndy<GeoTransform>;
        
    };





















 

 

  
} // OndoMathX



#include "CylinderBndy.hxx"

// -----------------------------------------------------------------------------------------//
 
