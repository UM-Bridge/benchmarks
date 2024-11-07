#pragma once
 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 
    template <
        typename GeoTransform = Identity<1>
    >
    class Line : public FESpaceT<Line<GeoTransform>,1>{
        
    
    private:
        
        //Order
        Index _PX;
        
        // Associated spectral finite element
        GaussLobattoElement _GLE;

        // Number of elements by color
        Index _NElemColor_0;
        Index _NElemColor_1;
 
        // Total number of elements
        Index _NX;

        // Number of local DoF
        Index _NLocDoF;

        // Number of DoF (total and the x-y direction)
        Index _NDoF;

        //Mesh step
        Real _hX;

        //Associated transformation of the template finite element space
        GeoTransform _GeoTransform;
        
    public:
        
            // ---------------------------------------------------------------------------------//
            /*! \brief Definition of the dimension as a static variable.*/
            static const Index Dim = 1;
           

      
        Line(Index PX, Index NX, GeoTransform aGeoTransform = IdentityMap1)
            : _GeoTransform(aGeoTransform)
        {
            // Order .
            _PX = PX;
            
            // Number of element .
            _NX = NX;

            // Number of degrees of freedom (DoF)
            _NDoF = (_NX * _PX + 1);

            // Number of local DoF.
            _NLocDoF = (_PX + 1);

            // Number of elements per color.
            _NElemColor_0 = (_NX / 2 + _NX % 2);
            _NElemColor_1 = (_NX / 2          );
   
            // Computing space step.
            _hX = 1.0 / _NX;
         
            // Initializing spectral finite element.
            _GLE = GaussLobattoElement(_PX+1,0,0);
        }
        // ---------------------------------------------------------------------------------//


        // ---------------------------------------------------------------------------------//
        // Extracting number of color in the template geometry, used for parallelization.
        Index getNumColors() const
        {
            return 2;
        }

        //  Extracting number of element associated to a given color.
        Index getNumElements(Index iColor) const
        {
            switch (iColor)
            {
                case 0: return _NElemColor_0;
                case 1: return _NElemColor_1;
            }
            
            {assert(false);}
            
            return _NX;
        }
        
        //  Extracting number of global number of element.
        Index getNumElements() const
        {
            return _NX;
        }
        
        const GaussLobattoElement & getFE() const
        {
            return _GLE;
        }
        
        Index getNumVertices() const
        {
            return (_NX+1);
        }
         
        void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
        {
            RealVector uvw;
            
            uvw[0]=iVertex*_hX;
            uvw[1]=0;
            uvw[2]=0;
            
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }
         
      
         
        void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            numbering.resize(2);
                     
            numbering[0] = (iEltGlob + 0);
            numbering[1] = (iEltGlob + 1);
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
            Index ex;
            Index gex;
 
            switch (c)
            {
                case 0: ex=e%(_NX/2+_NX%2); gex=2*ex;   break;
                case 1: ex=e%(_NX/2      ); gex=2*ex+1; break;
                default: assert(false);
            }
           
            return gex;
        }

        // Extracting global DoF index from local DoF index.
        // iEltGlob is the index of the element in the color group.
        // iLoc is the local DoF index.
        // It returns the global index.
        Index Loc2Glob(Index iEltGlob, Index iLoc) const
        {
            return (iEltGlob * _PX + iLoc);
        }
        
        void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob)
        {
            Index glob_ref =  iEltGlob * _PX;
            
            for (Index iLoc=0; iLoc < _NLocDoF; iLoc++)
                loc2glob[iLoc] = glob_ref + iLoc;
            
        }

 

        // ---------------------------------------------------------------------------------//
        /*! \brief Gives a degree of freedom coordinate associated to a global DoF.
            \param iGlob is a global index of an DoF.
            \param xyz are the coordinates of the degree of freedo massociated to iLoc which will be filled.
        */
        void getDoFCoordinate(const Index iGlob, RealVector& xyz) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;

            // Extracting DoF coordinates in the template finite element space.
            Index iEltGlob = iGlob / _PX;
            Index iLoc = iGlob - iEltGlob * _PX;

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector & P = _GLE.getPointCoordinate(iLoc);

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlob) * _hX;
            uvw[1] = 0.0;
            uvw[2] = 0.0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        

        /*! \brief Gives the barycenter of the element
            \param iEltGlob is a global index of an element.
            \param xyzCenter are the coordinates of the barycenter of the element
        */
        void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;

            uvw[0] = (0.5 + iEltGlob) * _hX;
            uvw[1] = 0.0;
            uvw[2] = 0.0;

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
            const RealVector& P = _GLE.getPointCoordinate(iLoc);
            
            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlob) * _hX;
            uvw[1] = 0;
            uvw[2] = 0;
            
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }
        
    



        // Extracting jacobian of the transformation used to define the template at a local DoF in an element
        Real getJacobian(Index iEltGlob, Index iLoc) const
        {
            // Storing template Quadrature Coordinates.
            std::array<Real, 1> uvw;

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector& P = _GLE.getPointCoordinate(iLoc);

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlob) * _hX;
   
            
            if constexpr(GeoTransform::isContinuouslyDifferentiable)
            {
                return _GeoTransform.getJacobian(uvw) * _hX;
            }
            else
            {
                RealVector dir;
                dir[0] = uvw[0] - (0.5 + iEltGlob) * _hX;
                
                return _GeoTransform.getJacobian(uvw, dir) * _hX;
            }
            
        }



        // Apply to a "gradient-like" finite element vector cobasis of the transformation used
        //  to define any element in the template at a local DoF in an element, i.e.
        //  \mathcal{F}(\mathrm{co}(F)) \nabla u \longrightarrow \nabla u,
        // where  F  is the jacobian matrix of the transformation and  u  is a solution defined
       // at a specific DoF
        
        void MltCoMatGradDef(Index iEltGlob, Index iLoc,
            std::array<Real, 1>& LocalSolutionGradient) const
        {
            // Storing template Quadrature Coordinates.
            std::array<Real, 1> uvw;

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector& P = _GLE.getPointCoordinate(iLoc);

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlob) * _hX;
            
            if constexpr(GeoTransform::isContinuouslyDifferentiable)
            {
                _GeoTransform.MltCoMatGradDef(uvw, LocalSolutionGradient);
            }
            else
            {
                RealVector dir;
                dir[0] = uvw[0] - (0.5 + iEltGlob) * _hX;
                _GeoTransform.MltCoMatGradDef(uvw, dir, LocalSolutionGradient);
            }
        }



        
        // Apply to a "gradient-like" finite element vector cobasis of the transformation used
        //  to define any element in the template at a local DoF in an element, i.e.
        //  \mathcal{F}(\mathrm{co}(F)^T) \nabla u \longrightarrow \nabla u,
        // where  F  is the jacobian matrix of the transformation and  u  is a solution defined
       // at a specific DoF
        
        void TransposeMltCoMatGradDef(Index iEltGlob, Index iLoc,
                          std::array<Real, 1>& LocalSolutionGradient) const
        {
            MltCoMatGradDef(iEltGlob,iLoc,LocalSolutionGradient);
        }
        

        // ---------------------------------------------------------------------------------//
        /*! \brief Gives a reference to the GeoTransform
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
        // ---------------------------------------------------------------------------------//
           
    };
 
}//OndoMathX
 
