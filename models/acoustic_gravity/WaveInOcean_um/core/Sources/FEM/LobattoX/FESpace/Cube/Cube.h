#pragma once
 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
    template <typename GeoTransform>  class CubeBndy;

    template <
        typename GeoTransform = Identity<3>
    >
	class Cube : public FESpaceT<Cube<GeoTransform>,3> {
 
    private:
        
        Square<Identity<2>> _Square;

        // Associated spectral finite element
        GaussLobattoElement _GLE;
     
        // Number of elements by color (also by directions)
        Index _NElemColor_0, _NElemColor_0_X, _NElemColor_0_Y, _NElemColor_0_XY, _NElemColor_0_Z;
        Index _NElemColor_1, _NElemColor_1_X, _NElemColor_1_Y, _NElemColor_1_XY, _NElemColor_1_Z;
        Index _NElemColor_2, _NElemColor_2_X, _NElemColor_2_Y, _NElemColor_2_XY, _NElemColor_2_Z;
        Index _NElemColor_3, _NElemColor_3_X, _NElemColor_3_Y, _NElemColor_3_XY, _NElemColor_3_Z;
        Index _NElemColor_4, _NElemColor_4_X, _NElemColor_4_Y, _NElemColor_4_XY, _NElemColor_4_Z;
        Index _NElemColor_5, _NElemColor_5_X, _NElemColor_5_Y, _NElemColor_5_XY, _NElemColor_5_Z;
        Index _NElemColor_6, _NElemColor_6_X, _NElemColor_6_Y, _NElemColor_6_XY, _NElemColor_6_Z;
        Index _NElemColor_7, _NElemColor_7_X, _NElemColor_7_Y, _NElemColor_7_XY, _NElemColor_7_Z;

        //Order of finite elements
        Index _PX;
        Index _PY;
        Index _PZ;
   
        
        
        // Number of elements
        Index _NX;
        Index _NY;
        Index _NZ;
        Index _NXY;

        // Total number of elements
        Index _NElem;

        // Number of local DoF
        Index _NLocDoF;
        Index _NLocDoFXY;
        Index _NLocDoFX;

        // Number of DoF
        Index _NDoF;
        Index _NDoFX;
        Index _NDoFY;
        Index _NDoFZ;
        Index _NDoFXY;

        // Storing mesh steps
        Real _hX, _hY, _hZ;

        // Associated transformation of the template finite element space
        GeoTransform _GeoTransform;
        
	public:

        
        // ---------------------------------------------------------------------------------//
        /*! \brief Definition of the dimension as a static variable.*/
        static const Index Dim = 3;
        // ---------------------------------------------------------------------------------//

        
        

		// ---------------------------------------------------------------------------------//
		/*! \brief Constructor of a finite element space with spectral finite element on a regular square.
			\param NX is the number of element in the direction X (uniformly distributed), must be even.
			\param NY is the number of element in the direction Y (uniformly distributed), must be even.
			\param NZ is the number of element in the direction Z (uniformly distributed), must be even.
			\param aGeoTransform is a transformation to apply on the unit square.
		*/
		Cube(Index PX, Index PY, Index PZ, Index NX, Index NY, Index NZ, GeoTransform aGeoTransform = IdentityMap3)
			: _Square(PX, PY, NX, NY), _GeoTransform(aGeoTransform)
		{
			// Number of element in every directions and in total.
			_NX = NX;
			_NY = NY;
			_NZ = NZ;
			_NXY = _NX * _NY;
			_NElem = _NX * _NY * _NZ;
            
            // Order of finite elements
            _PX = PX;
            _PY = PY;
            _PZ = PZ;
         
			 
			// Number of degrees of freedom (DoF) int the first direction and in total.
			_NDoFX = (_NX * _PX + 1);
			_NDoFY = (_NY * _PY + 1);
            _NDoFZ = (_NZ * _PZ + 1);
			_NDoFXY = _NDoFX * _NDoFY;
			_NDoF = _NDoFX * _NDoFY * _NDoFZ;

			// Number of local DoF.
            _NLocDoFX = (_PX + 1);
			_NLocDoFXY = (_PX + 1) * (_PY + 1);
			_NLocDoF = _NLocDoFXY * (_PZ + 1);

			// Number of elements per color.
            _NElemColor_0_X = (_NX / 2 + _NX % 2);
            _NElemColor_1_X = (_NX / 2          );
            _NElemColor_2_X = (_NX / 2 + _NX % 2);
            _NElemColor_3_X = (_NX / 2          );
            _NElemColor_4_X = (_NX / 2 + _NX % 2);
            _NElemColor_5_X = (_NX / 2          );
            _NElemColor_6_X = (_NX / 2 + _NX % 2);
            _NElemColor_7_X = (_NX / 2          );
            
            _NElemColor_0_Y = (_NY / 2 + _NY % 2);
            _NElemColor_1_Y = (_NY / 2 + _NY % 2);
            _NElemColor_2_Y = (_NY / 2)          ;
            _NElemColor_3_Y = (_NY / 2)          ;
            _NElemColor_4_Y = (_NY / 2 + _NY % 2);
            _NElemColor_5_Y = (_NY / 2 + _NY % 2);
            _NElemColor_6_Y = (_NY / 2)          ;
            _NElemColor_7_Y = (_NY / 2)          ;
            
            _NElemColor_0_Z = (_NZ / 2 + _NZ % 2);
            _NElemColor_1_Z = (_NZ / 2 + _NZ % 2);
            _NElemColor_2_Z = (_NZ / 2 + _NZ % 2);
            _NElemColor_3_Z = (_NZ / 2 + _NZ % 2);
            _NElemColor_4_Z = (_NZ / 2)          ;
            _NElemColor_5_Z = (_NZ / 2)          ;
            _NElemColor_6_Z = (_NZ / 2)          ;
            _NElemColor_7_Z = (_NZ / 2)          ;
            
            _NElemColor_0_XY = _NElemColor_0_X * _NElemColor_0_Y;
            _NElemColor_1_XY = _NElemColor_1_X * _NElemColor_1_Y;
            _NElemColor_2_XY = _NElemColor_2_X * _NElemColor_2_Y;
            _NElemColor_3_XY = _NElemColor_3_X * _NElemColor_3_Y;
            _NElemColor_4_XY = _NElemColor_4_X * _NElemColor_4_Y;
            _NElemColor_5_XY = _NElemColor_5_X * _NElemColor_5_Y;
            _NElemColor_6_XY = _NElemColor_6_X * _NElemColor_6_Y;
            _NElemColor_7_XY = _NElemColor_7_X * _NElemColor_7_Y;
            
            _NElemColor_0 = _NElemColor_0_XY * _NElemColor_0_Z;
            _NElemColor_1 = _NElemColor_1_XY * _NElemColor_1_Z;
            _NElemColor_2 = _NElemColor_2_XY * _NElemColor_2_Z;
            _NElemColor_3 = _NElemColor_3_XY * _NElemColor_3_Z;
            _NElemColor_4 = _NElemColor_4_XY * _NElemColor_4_Z;
            _NElemColor_5 = _NElemColor_5_XY * _NElemColor_5_Z;
            _NElemColor_6 = _NElemColor_6_XY * _NElemColor_6_Z;
            _NElemColor_7 = _NElemColor_7_XY * _NElemColor_7_Z;

            
			// Computing space step.
			_hX = 1.0 / _NX;
			_hY = 1.0 / _NY;
			_hZ = 1.0 / _NZ;

            // Initializing spectral finite element.
            _GLE = GaussLobattoElement(_PX+1,_PY+1,_PZ+1);
     
    

		}
		// ---------------------------------------------------------------------------------//



        // ---------------------------------------------------------------------------------//
        // Extracting number of color in the template geometry, used for parallelization.
        template<bool Continuous = true> Index getNumColors() const
        {
            return 8;
        }
        
        //  Extracting number of element associated to a given color.
        Index getNumElements(Index iColor) const
        {
            switch (iColor)
            {
                case 0: return _NElemColor_0;
                case 1: return _NElemColor_1;
                case 2: return _NElemColor_2;
                case 3: return _NElemColor_3;
                case 4: return _NElemColor_4;
                case 5: return _NElemColor_5;
                case 6: return _NElemColor_6;
                case 7: return _NElemColor_7;
            }    
            {assert(false);}
        }
                
        //  Extracting number of global number of element.
        Index getNumElements() const
        {
            return _NElem;
        }

        bool isAffine(const Index iEltGlob)
        {
            return _GeoTransform.isAffine();
        }
        
        const GaussLobattoElement & getFE() const
        {
            return _GLE;
        }
        
        Index getNumVertices() const
        {
            return (_NX+1)*(_NY+1)*(_NZ+1);
        }

        Index getNumVerticesQ2() const
        {
            return (_NX*2+1)*(_NY*2+1)*(_NZ*2+1);
        }
         
        void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
        {
            RealVector uvw;
            
            Index iVertexZ = iVertex / ((_NX+1)*(_NY+1));
            Index iVertexX = (iVertex - iVertexZ * (_NX+1)*(_NY+1) ) % (_NX+1);
            Index iVertexY = (iVertex - iVertexZ * (_NX+1)*(_NY+1) ) / (_NX+1);
          
        
            uvw[0]=iVertexX*_hX;
            uvw[1]=iVertexY*_hY;
            uvw[2]=iVertexZ*_hZ;
            
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        void getVertexCoordinateQ2(const Index iVertex, RealVector& xyz) const
        {
            Index NVX = (_NX * 2 + 1);
			Index NVY = (_NY * 2 + 1);
			Index NVXY = NVX * NVY;

			// Extracting DoF coordinates in the template finite element space.
			Index iGlobZ = iVertex / NVXY;
			Index iGlobX = (iVertex - iGlobZ * NVXY) % NVX;
			Index iGlobY = (iVertex - iGlobZ * NVXY) / NVX;

			Index iEltGlobX = iGlobX / 2;
			Index iEltGlobY = iGlobY / 2;
			Index iEltGlobZ = iGlobZ / 2;

			Index iLocX = iGlobX - iEltGlobX * 2;
			Index iLocY = iGlobY - iEltGlobY * 2;
			Index iLocZ = iGlobZ - iEltGlobZ * 2;

            Real u = iLocX*0.5;
            Real v = iLocY*0.5;
            Real w = iLocZ*0.5;

			// Adjusting quadrature point coordinate from the element index.
            RealVector  uvw;

			uvw[0] = (u + iEltGlobX) * _hX;
			uvw[1] = (v + iEltGlobY) * _hY;
			uvw[2] = (w + iEltGlobZ) * _hZ;

			// Applying transformation.
			_GeoTransform.Eval(uvw, xyz);
        }
         
        void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            numbering.resize(8);
            
            Index iEltGlobZ = iEltGlob / _NXY;
            Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
            Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;
            
            numbering[0] = (iEltGlobZ + 0) * (_NX+1) * (_NY+1) + (iEltGlobY + 0) * (_NX+1) + (iEltGlobX + 0);
            numbering[1] = (iEltGlobZ + 0) * (_NX+1) * (_NY+1) + (iEltGlobY + 0) * (_NX+1) + (iEltGlobX + 1);
            numbering[2] = (iEltGlobZ + 0) * (_NX+1) * (_NY+1) + (iEltGlobY + 1) * (_NX+1) + (iEltGlobX + 1);
            numbering[3] = (iEltGlobZ + 0) * (_NX+1) * (_NY+1) + (iEltGlobY + 1) * (_NX+1) + (iEltGlobX + 0);
            numbering[4] = (iEltGlobZ + 1) * (_NX+1) * (_NY+1) + (iEltGlobY + 0) * (_NX+1) + (iEltGlobX + 0);
            numbering[5] = (iEltGlobZ + 1) * (_NX+1) * (_NY+1) + (iEltGlobY + 0) * (_NX+1) + (iEltGlobX + 1);
            numbering[6] = (iEltGlobZ + 1) * (_NX+1) * (_NY+1) + (iEltGlobY + 1) * (_NX+1) + (iEltGlobX + 1);
            numbering[7] = (iEltGlobZ + 1) * (_NX+1) * (_NY+1) + (iEltGlobY + 1) * (_NX+1) + (iEltGlobX + 0);
        }

        void getElementVerticesQ2(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            numbering.resize(27);
            
            Index iEltGlobZ = iEltGlob / _NXY;
            Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
            Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;

            Index iEltGlobXY = iEltGlobX + iEltGlobY*_NX;

            Index NVertSurfaceQ2 = _Square.getNumVerticesQ2();
            
            std::vector<Index> numberingXY;

            _Square.getElementVerticesQ2(iEltGlobXY, numberingXY);

            numbering[0]  = numberingXY[0]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[8]  = numberingXY[4]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[1]  = numberingXY[1]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[9]  = numberingXY[7]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[20] = numberingXY[8]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[11] = numberingXY[5]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[3]  = numberingXY[3]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[13] = numberingXY[6]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;
            numbering[2]  = numberingXY[2]+ (iEltGlobZ*2 + 0) * NVertSurfaceQ2;

            numbering[10] = numberingXY[0] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[21] = numberingXY[4] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[12] = numberingXY[1] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[22] = numberingXY[7] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[26] = numberingXY[8] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[23] = numberingXY[5] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[15] = numberingXY[3] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[24] = numberingXY[6] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;
            numbering[14] = numberingXY[2] + (iEltGlobZ*2 + 1) * NVertSurfaceQ2;

            numbering[4]  = numberingXY[0] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[16] = numberingXY[4] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[5]  = numberingXY[1] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[17] = numberingXY[7] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[25] = numberingXY[8] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[18] = numberingXY[5] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[7]  = numberingXY[3] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[19] = numberingXY[6] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;
            numbering[6]  = numberingXY[2] + (iEltGlobZ*2 + 2) * NVertSurfaceQ2 ;


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
        // ---------------------------------------------------------------------------------//


        // Fonction that maps the numbering of the element by color to a global numbering.
        // c is the index of the color.
        // e is the numbering of the element in the color e.
        // It returns the global numbering of the element.
        Index EltColorNum2EltGlobalNum(Index c, Index e) const
        {
            Index ex,ey,ez;
            Index gex,gey,gez;

            switch (c)
            {
                case 0: ex=e%_NElemColor_0_X; ey=(e%_NElemColor_0_XY)/_NElemColor_0_X; ez=e/(_NElemColor_0_XY); gex=2*ex;   gey=2*ey;   gez=2*ez;   break;
                case 1: ex=e%_NElemColor_1_X; ey=(e%_NElemColor_1_XY)/_NElemColor_1_X; ez=e/(_NElemColor_1_XY); gex=2*ex+1; gey=2*ey;   gez=2*ez;   break;
                case 2: ex=e%_NElemColor_2_X; ey=(e%_NElemColor_2_XY)/_NElemColor_2_X; ez=e/(_NElemColor_2_XY); gex=2*ex  ; gey=2*ey+1; gez=2*ez;   break;
                case 3: ex=e%_NElemColor_3_X; ey=(e%_NElemColor_3_XY)/_NElemColor_3_X; ez=e/(_NElemColor_3_XY); gex=2*ex+1; gey=2*ey+1; gez=2*ez;   break;
                case 4: ex=e%_NElemColor_4_X; ey=(e%_NElemColor_4_XY)/_NElemColor_4_X; ez=e/(_NElemColor_4_XY); gex=2*ex;   gey=2*ey;   gez=2*ez+1; break;
                case 5: ex=e%_NElemColor_5_X; ey=(e%_NElemColor_5_XY)/_NElemColor_5_X; ez=e/(_NElemColor_5_XY); gex=2*ex+1; gey=2*ey;   gez=2*ez+1; break;
                case 6: ex=e%_NElemColor_6_X; ey=(e%_NElemColor_6_XY)/_NElemColor_6_X; ez=e/(_NElemColor_6_XY); gex=2*ex  ; gey=2*ey+1; gez=2*ez+1; break;
                case 7: ex=e%_NElemColor_7_X; ey=(e%_NElemColor_7_XY)/_NElemColor_7_X; ez=e/(_NElemColor_7_XY); gex=2*ex+1; gey=2*ey+1; gez=2*ez+1; break;
                default: assert(false);
            }
        
            return gez*_NXY+gey*_NX+gex;

        }
        
        // Extracting global DoF index from local DoF index.
        // iEltGlob is the index of the element in the color group.
        // iLoc is the local DoF index.
        // It returns the global index.
        Index _Loc2Glob(Index iEltGlob, Index iLoc) const
		{
			Index iEltGlobZ = iEltGlob / _NXY;
			Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
			Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;

			Index iLocZ = iLoc / _NLocDoFXY;
			Index iLocX = (iLoc - iLocZ * _NLocDoFXY) % (_PX + 1);
			Index iLocY = (iLoc - iLocZ * _NLocDoFXY) / (_PX + 1);

			return (iEltGlobZ * _PZ + iLocZ) * _NDoFXY + (iEltGlobY*_PY + iLocY) * _NDoFX + (iEltGlobX * _PX + iLocX);
		}
        
/*
        void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
        {
            Index iEltGlobZ = iEltGlob / _NXY;
            Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
            Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;
            
            Index glob_ref = iEltGlobZ * _PZ * _NDoFXY + iEltGlobY*_PY * _NDoFX + (iEltGlobX * _PX);
            
            for (Index iLoc=0; iLoc < _NLocDoF; iLoc++)
            {
                Index iLocZ = iLoc / _NLocDoFXY;
                Index iLocX = (iLoc - iLocZ * _NLocDoFXY) % (_PX + 1);
                Index iLocY = (iLoc - iLocZ * _NLocDoFXY) / (_PX + 1);
                
                loc2glob[iLoc] = glob_ref +  iLocZ * _NDoFXY +  iLocY * _NDoFX +  iLocX;
            }
        }
*/

       
		// ---------------------------------------------------------------------------------//


    
        

		// ---------------------------------------------------------------------------------//
		/*! \brief Gives the barycenter of the element
			\param iEltGlob is a global index of an element.
			\param xyzCenter are the coordinates of the barycenter of the element
		*/
		void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
		{
			Index iEltGlobZ = iEltGlob / _NXY;
			Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
			Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;

			// Storing template Quadrature Coordinates.
            RealVector uvw;

			uvw[0] = (0.5 + iEltGlobX) * _hX;
			uvw[1] = (0.5 + iEltGlobY) * _hY;
			uvw[2] = (0.5 + iEltGlobZ) * _hZ;

			// Applying transformation.
			_GeoTransform.Eval(uvw, xyzCenter);
		}


        // Get the label of the element (here 0)
        Index getLabel(const Index iEltGlob) const
        {
            return 0;
        }
        

		/*! \brief Gives a degree of freedom coordinate associated to a global DoF.
			\param iGlob is a global index of an DoF.
			\param xyz are the coordinates of the degree of freedo massociated to iLoc which will be filled.
		*/
		void getDoFCoordinate(const Index iGlob, RealVector& xyz) const
		{
			// Storing template Quadrature Coordinates.
            RealVector  uvw;

			// Extracting DoF coordinates in the template finite element space.
			Index iGlobZ = iGlob / _NDoFXY;
			Index iGlobX = (iGlob - iGlobZ * _NDoFXY) % _NDoFX;
			Index iGlobY = (iGlob - iGlobZ * _NDoFXY) / _NDoFX;

			Index iEltGlobX = iGlobX / _PX;
			Index iEltGlobY = iGlobY / _PY;
			Index iEltGlobZ = iGlobZ / _PZ;

			Index iLocX = iGlobX - iEltGlobX * _PX;
			Index iLocY = iGlobY - iEltGlobY * _PY;
			Index iLocZ = iGlobZ - iEltGlobZ * _PZ;

			Index iLoc = iLocZ * _NLocDoFXY + iLocY * (_PX + 1) + iLocX;

			// Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector & P = _GLE.getPointCoordinate(iLoc);

			// Adjusting quadrature point coordinate from the element index.
			uvw[0] = (P[0] + iEltGlobX) * _hX;
			uvw[1] = (P[1] + iEltGlobY) * _hY;
			uvw[2] = (P[2] + iEltGlobZ) * _hZ;

			// Applying transformation.
			_GeoTransform.Eval(uvw, xyz);
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

			// Extracting quadrature coordinates in the template finite element space.
			Index iEltGlobZ = iEltGlob / _NXY;
			Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
			Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;

			// Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector & P = _GLE.getPointCoordinate(iLoc);

			// Adjusting quadrature point coordinate from the element index.
			uvw[0] = (P[0] + iEltGlobX) * _hX;
			uvw[1] = (P[1] + iEltGlobY) * _hY;
			uvw[2] = (P[2] + iEltGlobZ) * _hZ;

			// Applying transformation.
			_GeoTransform.Eval(uvw, xyz);
		}
        
        

        void _getGradDef(Index iEltGlob, Index iLoc,
			            std::array<std::array<Real, 3>,3>& GradDef) const
            {
                // Storing template Quadrature Coordinates.
                RealVector uvw;
                
                // Extracting quadrature coordinates in the template finite element space.
                Index iEltGlobZ = iEltGlob / _NXY;
                Index iEltGlobX = (iEltGlob - iEltGlobZ * _NXY) % _NX;
                Index iEltGlobY = (iEltGlob - iEltGlobZ * _NXY) / _NX;
                
                // Extracting corresponding quadrature point coordinate in the reference element.
                const RealVector & P = _GLE.getPointCoordinate(iLoc);
                
                // Adjusting quadrature point coordinate from the element index.
                uvw[0] = (P[0] + iEltGlobX) * _hX;
                uvw[1] = (P[1] + iEltGlobY) * _hY;
                uvw[2] = (P[2] + iEltGlobZ) * _hZ;

                // Applying transform co-basis.
                if constexpr(GeoTransform::isContinuouslyDifferentiable)
                {
                    _GeoTransform.getGradDef(uvw, GradDef);
                }
                else
                {
                    RealVector dir;
                    
                    dir[0] = uvw[0] - (0.5 + iEltGlobX) * _hX;
                    dir[1] = uvw[1] - (0.5 + iEltGlobY) * _hY;
                    dir[2] = uvw[2] - (0.5 + iEltGlobZ) * _hZ;
                    
                    _GeoTransform.getGradDef(uvw, dir, GradDef);
                }

                GradDef[0][0] *= _hX; GradDef[0][1] *= _hY; GradDef[0][2] *= _hZ;
                GradDef[1][0] *= _hX; GradDef[1][1] *= _hY; GradDef[1][2] *= _hZ;
                GradDef[2][0] *= _hX; GradDef[2][1] *= _hY; GradDef[2][2] *= _hZ;
            }


		// ---------------------------------------------------------------------------------//
 
 

    


        
        // ---------------------------------------------------------------------------------//
        /*! \brief Gives a copy to the GeoTransform
         */
        GeoTransform getGeoTransform() const
        {
            return  _GeoTransform;
        }
        
        /*! \brief Return the number of element in the direction X
         */
        Index getNx() const
        {
            return  _NX;
        }
        
        /*! \brief Return the number of element in the direction X
         */
        Index getNy() const
        {
            return  _NY;
        }
        
        /*! \brief Return the number of element in the direction X
         */
        Index getNz() const
        {
            return  _NZ;
        }
        // ---------------------------------------------------------------------------------//

        CubeBndy<GeoTransform> getBndyFESpace(Index iBoundary)
        {
            assert(iBoundary<6);
            
            return CubeBndy<GeoTransform>(iBoundary,*this);
        }
        
        friend class CubeBndy<GeoTransform>;

	};

} // OndoMathX


#include "CubeBndy.hxx"

// -----------------------------------------------------------------------------------------//
 
