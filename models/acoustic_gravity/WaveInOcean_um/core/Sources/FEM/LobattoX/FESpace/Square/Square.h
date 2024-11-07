#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{

    enum PeriodicityType
    {
        None = 0,
        X = 1,
        Y = 2,
        XY = 3
    };

    template <typename GeoTransform>
    class SquareBndy;
    template <typename GeoTransform>
    class SquareSkeleton;

    template <
        typename GeoTransform = Identity<2>>
    class Square : public FESpaceT<Square<GeoTransform>, 2>
    {

    private:
        // Associated spectral finite element
        GaussLobattoElement _GLE;

        // Number of elements by color
        Index _NElemColor_0;
        Index _NElemColor_1;
        Index _NElemColor_2;
        Index _NElemColor_3;

        // Order of the elements in the x-direction and y-direction
        Index _PX;
        Index _PY;

        // Number of element in the x-direction and y-direction
        Index _NX;
        Index _NY;

        // Total number of elements
        Index _NElem;

        // Number of local DoF
        Index _NLocDoFX;
        Index _NLocDoFY;
        Index _NLocDoF;

        // Number of DoF (total and the x-y direction)
        Index _NDoF;
        Index _NDoFX;
        Index _NDoFY;

        // Mesh step
        Real _hX, _hY;

        // Associated transformation of the template finite element space
        GeoTransform _GeoTransform;

    public:
        // ---------------------------------------------------------------------------------//
        /*! \brief Definition of the dimension as a static variable.*/
        static const Index Dim = 2;

        Square(Index PX, Index PY, Index NX, Index NY, GeoTransform aGeoTransform = IdentityMap2) : _GeoTransform(aGeoTransform)
        {
            // Order of the method
            _PX = PX;
            _PY = PY;

            // Number of element in every directions and in total.
            _NX = NX;
            _NY = NY;
            _NElem = _NX * _NY;

            // Number of degrees of freedom (DoF) in the 2 directions and in total.
            _NDoFX = (_NX * _PX + 1);
            _NDoFY = (_NY * _PY + 1);
            _NDoF = _NDoFX * _NDoFY;

            // Number of local DoF.
            _NLocDoFX = (_PX + 1);
            _NLocDoFY = (_PY + 1);
            _NLocDoF = (_PX + 1) * (_PY + 1);

            // Number of elements per color.
            _NElemColor_0 = (_NX / 2 + _NX % 2) * (_NY / 2 + _NY % 2);
            _NElemColor_1 = (_NX / 2) * (_NY / 2 + _NY % 2);
            _NElemColor_2 = (_NX / 2 + _NX % 2) * (_NY / 2);
            _NElemColor_3 = (_NX / 2) * (_NY / 2);

            // Computing space step.
            _hX = 1.0 / _NX;
            _hY = 1.0 / _NY;

            // Initializing spectral finite element.
            _GLE = GaussLobattoElement(_PX + 1, _PY + 1, 0);
        }
        // ---------------------------------------------------------------------------------//

        // ---------------------------------------------------------------------------------//
        // Extracting number of color in the template geometry, used for parallelization.
        Index getNumColors() const
        {
            return 4;
        }

        //  Extracting number of element associated to a given color.
        Index getNumElements(Index iColor) const
        {
            switch (iColor)
            {
            case 0:
                return _NElemColor_0;
            case 1:
                return _NElemColor_1;
            case 2:
                return _NElemColor_2;
            case 3:
                return _NElemColor_3;
            }

            {
                assert(false);
            }

            return _NElem;
        }

        //  Extracting number of global number of element.
        Index getNumElements() const
        {
            return _NElem;
        }

        const GaussLobattoElement &getFE() const
        {
            return _GLE;
        }

        bool isAffine(const Index iEltGlob) const
        {
            return _GeoTransform.isAffine();
        }

        Index getNumVertices() const
        {
            return (_NX + 1) * (_NY + 1);
        }

        Index getNumVerticesQ2() const
        {
            return (_NX * 2 + 1) * (_NY * 2 + 1);
        }

        void getVertexCoordinate(const Index iVertex, RealVector &xyz) const
        {
            RealVector uvw;

            Index iVertexX = iVertex % (_NX + 1);
            Index iVertexY = iVertex / (_NX + 1);

            uvw[0] = iVertexX * _hX;
            uvw[1] = iVertexY * _hY;
            uvw[2] = 0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        void getElementVertices(const Index iEltGlob, std::vector<Index> &numbering) const
        {
            numbering.resize(4);

            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            numbering[0] = (iEltGlobY + 0) * (_NX + 1) + (iEltGlobX + 0);
            numbering[1] = (iEltGlobY + 0) * (_NX + 1) + (iEltGlobX + 1);
            numbering[2] = (iEltGlobY + 1) * (_NX + 1) + (iEltGlobX + 1);
            numbering[3] = (iEltGlobY + 1) * (_NX + 1) + (iEltGlobX + 0);
        }
        // ---------------------------------------------------------------------------------//

        void getVertexCoordinateQ2(const Index iVertex, RealVector &xyz) const
        {
            Index NVX = (_NX * 2 + 1);

            Index iGlobX = iVertex % NVX;
            Index iGlobY = iVertex / NVX;

            Index iEltGlobX = iGlobX / 2;
            Index iEltGlobY = iGlobY / 2;

            Index iLocX = iGlobX - iEltGlobX * 2;
            Index iLocY = iGlobY - iEltGlobY * 2;

            Real u = iLocX * 0.5;
            Real v = iLocY * 0.5;

            // Storing template Quadrature Coordinates.
            RealVector uvw;

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (u + iEltGlobX) * _hX;
            uvw[1] = (v + iEltGlobY) * _hY;
            uvw[2] = 0.0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        void getElementVerticesQ2(const Index iEltGlob, std::vector<Index> &numbering) const
        {

            numbering.resize(9);

            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            numbering[0] = (2 * iEltGlobY + 0) * (2 * _NX + 1) + (2 * iEltGlobX + 0);
            numbering[1] = (2 * iEltGlobY + 0) * (2 * _NX + 1) + (2 * iEltGlobX + 2);
            numbering[2] = (2 * iEltGlobY + 2) * (2 * _NX + 1) + (2 * iEltGlobX + 2);
            numbering[3] = (2 * iEltGlobY + 2) * (2 * _NX + 1) + (2 * iEltGlobX + 0);

            numbering[4] = (2 * iEltGlobY + 0) * (2 * _NX + 1) + (2 * iEltGlobX + 1);
            numbering[5] = (2 * iEltGlobY + 1) * (2 * _NX + 1) + (2 * iEltGlobX + 2);
            numbering[6] = (2 * iEltGlobY + 2) * (2 * _NX + 1) + (2 * iEltGlobX + 1);
            numbering[7] = (2 * iEltGlobY + 1) * (2 * _NX + 1) + (2 * iEltGlobX + 0);

            numbering[8] = (2 * iEltGlobY + 1) * (2 * _NX + 1) + (2 * iEltGlobX + 1);
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

        // Order in the X-direction
        Index getPX() const
        {
            return _PX;
        }

        // Order in the Y-direction
        Index getPY() const
        {
            return _PY;
        }
        // ---------------------------------------------------------------------------------//

        // Fonction that maps the numbering of the element by color to a global numbering.
        // c is the index of the color.
        // e is the numbering of the element in the color e.
        // It returns the global numbering of the element.
        Index EltColorNum2EltGlobalNum(Index c, Index e) const
        {
            Index ex, ey;
            Index gex, gey;

            switch (c)
            {
            case 0:
                ex = e % (_NX / 2 + _NX % 2);
                ey = e / (_NX / 2 + _NX % 2);
                gex = 2 * ex;
                gey = 2 * ey;
                break;
            case 1:
                ex = e % (_NX / 2);
                ey = e / (_NX / 2);
                gex = 2 * ex + 1;
                gey = 2 * ey;
                break;
            case 2:
                ex = e % (_NX / 2 + _NX % 2);
                ey = e / (_NX / 2 + _NX % 2);
                gex = 2 * ex;
                gey = 2 * ey + 1;
                break;
            case 3:
                ex = e % (_NX / 2);
                ey = e / (_NX / 2);
                gex = 2 * ex + 1;
                gey = 2 * ey + 1;
                break;
            default:
                assert(false);
            }

            return gey * _NX + gex;
        }

        Index EltGlobalNum2Color(Index eg) const
        {
            return (eg / _NX) % 2 + 2 * ((eg % _NX) % 2);
        }

        // Extracting global DoF index from local DoF index.
        // iEltGlob is the index of the element in the color group.
        // iLoc is the local DoF index.
        // It returns the global index.
        Index _Loc2Glob(Index iEltGlob, Index iLoc) const
        {
            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            Index iLocX = iLoc % (_PX + 1);
            Index iLocY = iLoc / (_PX + 1);

            return (iEltGlobY * _PY + iLocY) * _NDoFX + (iEltGlobX * _PX + iLocX);
        }

        /*void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
        {
            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            Index glob_ref =  iEltGlobY * _PY * _NDoFX + iEltGlobX * _PX;

            for (Index iLocY=0; iLocY < _NLocDoFY; iLocY++)
            {
                Index tmp = glob_ref +  iLocY * _NDoFX;
                Index tmp_iLoc = iLocY*_NLocDoFX;

                for (Index iLocX=0; iLocX < _NLocDoFX; iLocX++)
                {
                    loc2glob[tmp_iLoc + iLocX] = tmp + iLocX;
                }
            }
        }*/

        // ---------------------------------------------------------------------------------//
        /*! \brief Gives a degree of freedom coordinate associated to a global DoF.
            \param iGlob is a global index of an DoF.
            \param xyz are the coordinates of the degree of freedo massociated to iLoc which will be filled.
        */
        void getDoFCoordinate(const Index iGlob, RealVector &xyz) const
        {
            Index iEltGlobX;
            Index iEltGlobY;
            Index iGlobX;
            Index iGlobY;
            Index iLocX;
            Index iLocY;
            Index iLoc;

            // if (Continuous)
            //{
            //  Extracting DoF coordinates in the template finite element space.
            iGlobX = iGlob % _NDoFX;
            iGlobY = iGlob / _NDoFX;

            iEltGlobX = iGlobX / _PX;
            iEltGlobY = iGlobY / _PY;

            iLocX = iGlobX - iEltGlobX * _PX;
            iLocY = iGlobY - iEltGlobY * _PY;
            // }
            // else
            //  {
            //      iGlobX = iGlob % (iEltGlobX*(_PX+1));
            //      iGlobY = iGlob / (iEltGlobX*(_PX+1));

            //      iEltGlobX = iGlobX / (_PX+1);
            //	    iEltGlobY = iGlobY / (_PY+1);

            //      iLocX = iGlobX - iEltGlobX * (_PX+1);
            //	    iLocY = iGlobY - iEltGlobY * (_PY+1);
            // }

            iLoc = iLocX + iLocY * (_PX + 1);

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector &P = _GLE.getPointCoordinate(iLoc);

            // Storing template Quadrature Coordinates.
            RealVector uvw;

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlobX) * _hX;
            uvw[1] = (P[1] + iEltGlobY) * _hY;
            uvw[2] = 0.0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        /*! \brief Gives the barycenter of the element
            \param iEltGlob is a global index of an element.
            \param xyzCenter are the coordinates of the barycenter of the element
        */
        void getBarycenter(const Index iEltGlob, RealVector &xyzCenter) const
        {
            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            // Storing template Quadrature Coordinates.
            RealVector uvw;

            uvw[0] = (0.5 + iEltGlobX) * _hX;
            uvw[1] = (0.5 + iEltGlobY) * _hY;
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
        void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector &xyz) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;

            // Extracting quadrature coordinates in the template finite element space.
            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector &P = _GLE.getPointCoordinate(iLoc);

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlobX) * _hX;
            uvw[1] = (P[1] + iEltGlobY) * _hY;
            uvw[2] = 0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        void _getGradDef(Index iEltGlob, Index iLoc,
                         RealMatrix2x2 &GradDef) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;

            // Extracting quadrature coordinates in the template finite element space.
            Index iEltGlobX = iEltGlob % _NX;
            Index iEltGlobY = iEltGlob / _NX;

            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector &P = _GLE.getPointCoordinate(iLoc);

            // Adjusting quadrature point coordinate from the element index.
            uvw[0] = (P[0] + iEltGlobX) * _hX;
            uvw[1] = (P[1] + iEltGlobY) * _hY;

            if constexpr (GeoTransform::isContinuouslyDifferentiable)
            {
                _GeoTransform.getGradDef(uvw, GradDef);
            }
            else
            {
                RealVector dir;

                dir[0] = uvw[0] - (0.5 + iEltGlobX) * _hX;
                dir[1] = uvw[1] - (0.5 + iEltGlobY) * _hY;
                dir[2] = 0.0;

                _GeoTransform.getGradDef(uvw, dir, GradDef);
            }

            GradDef[0][0] *= _hX;
            GradDef[0][1] *= _hY;
            GradDef[1][0] *= _hX;
            GradDef[1][1] *= _hY;
        }

        // ---------------------------------------------------------------------------------/

        /*! \brief Gives a copy to the GeoTransform
         */
        GeoTransform getGeoTransform() const
        {
            return _GeoTransform;
        }

        /*! \brief Return the number of element in the direction X
         */
        Index getNX() const
        {
            return _NX;
        }

        /*! \brief Return the number of element in the direction X
         */
        Index getNY() const
        {
            return _NY;
        }

        // ---------------------------------------------------------------------------------//

        // ---------------------------------------------------------------------------------//

        SquareBndy<GeoTransform> getBndyFESpace(Index iBoundary)
        {
            assert(iBoundary < 4);

            return SquareBndy<GeoTransform>(iBoundary, *this);
        }

        SquareSkeleton<GeoTransform> getSkeletonFESpace(PeriodicityType PerType = PeriodicityType::None)
        {
            return SquareSkeleton<GeoTransform>(*this, PerType);
        }

        // ---------------------------------------------------------------------------------//

        friend class SquareBndy<GeoTransform>;
        friend class SquareSkeleton<GeoTransform>;
    };

} // OndoMathX

#include "SquareBndy.hxx"
#include "SquareSkeleton.hxx"
