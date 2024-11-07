#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {


    template <
        typename GeoTransform = Identity<2>
    >
    class Circle  : public FESpaceT<Circle<GeoTransform>,2> {
        
        template <typename GT> friend class Cylinder;
        template <typename GT> friend class CylinderBndy;

    private:
        
        // Associated spectral finite element
        GaussLobattoElement _GLE;

        //Order
        Index _P;

        // Number of elements by color
        Index _NElemColor;

        //Number of element in the x-direction and y-direction
        Index _NX;
        Index _NRho;
        Index _NTheta;

        // Total number of elements
        Index _NElem;
        Index _NSquareElem;
        Index _NCircleElem;

        // Number of local DoF
        Index _NLocDoF;

        // Number of DoF (total and the x-y direction)
        Index _NDoF;
        Index _NDoFSquare;
        Index _NDoFX;
        Index _NDoFRho;
        Index _NDoFTheta;

        //Mesh step
        Real _hX, _hTheta, _L;

        //Associated transformation of the template finite element space
        GeoTransform _GeoTransform;
        
        //buffer for pre-computations of the geometric informations
        std::vector<std::vector<RealMatrix2x2>> _GradDef;
        
        bool _GeomInfoPreComputed;
        
        
        
        
    public:
        

        //Class methods
        Index _getNumVerticesIntSquare() const
        {
            return (_NX+1)*(_NX+1);
        }

        Index _getNumVerticesIntSquareQ2() const
        {
            Index NVertQ2 = 2;
            return (_NX*NVertQ2+1)*(_NX*NVertQ2+1);
        }
        
         Index _getNumVerticesExtCircle() const
        {
            return _NRho * _NTheta ;
        }

        Index _getNumVerticesExtCircleQ2() const
        {
            Index NVertQ2 = 2;
            return (_NRho * NVertQ2) * (_NTheta * NVertQ2);
        }

        Index _getNumElementsinSquare(Index iColor) const
        {
            return _NX * _NX / 4;
        }

        void _Glob2Loc(const Index iGlob, Index & iLoc, Index & iEltGlob ) const
        {
            if (iGlob < _NDoFSquare)
            {
                // Extracting DoF coordinates in the template finite element space.
                Index iGlobX = iGlob % _NDoFX;
                Index iGlobY = iGlob / _NDoFX;
                Index iEltGlobX = iGlobX / _P;
                Index iEltGlobY = iGlobY / _P;
                if(iGlobX == (_NDoFX - 1))
                    iEltGlobX = iGlobX / _P - 1;
                
                if(iGlobY == (_NDoFX - 1))
                     iEltGlobY = iGlobY / _P - 1;
                 
                iEltGlob = iEltGlobX + iEltGlobY * _NX;
            
                Index iLocX = iGlobX - iEltGlobX * _P;
                Index iLocY = iGlobY - iEltGlobY * _P;

                iLoc = iLocX + iLocY * (_P + 1);
            }
            else
            {
                Index iGlobRel = iGlob - _NDoFSquare;
                
                // Extracting DoF coordinates in the template finite element space.
                Index iGlobTheta = iGlobRel % _NDoFTheta;
                Index iGlobRho = iGlobRel / _NDoFTheta;

                //Element coordinate along x
                Index iEltGlobTheta = iGlobTheta / _P; 
                //Element coordinate along y
                Index iEltGlobRho = (iGlobRho + 1) / _P; 
                
                iEltGlob = (iEltGlobTheta + iEltGlobRho * _NTheta) + _NSquareElem;
                
                Index iLocTheta = iGlobTheta - iEltGlobTheta * _P;
                Index iLocRho = iGlobRho - iEltGlobRho * _P + 1;

                iLoc = iLocRho + iLocTheta * (_P + 1);
             
            }
            
        }
        
        Index _Loc2GlobCircle(const Index iEltRel, const Index iLocX, const Index iLocY) const
        {
            Index iEltTheta = iEltRel % _NTheta;
            Index iEltRho = iEltRel / _NTheta;
            
            if (iLocX==0 && iEltRho==0) //CRITICAL LAYER OF DOF
            {
                
                Index iEltSquareEdge = iEltTheta % _NX;

                Index iQuarter = iEltRel / _NX;
                
                //first edge of the square (south)
                if (iQuarter == 0)
                    return iEltTheta*_P + iLocY;
    
                //Second edge of the square  (east)
                else if (iQuarter == 1)
                    return (iLocY + 1 + iEltSquareEdge*_P) * _NDoFX -1 ;
                
                // Third edge (north)
                else if (iQuarter == 2)
                    return _NDoFSquare - (iLocY + 1 + iEltSquareEdge*_P);
                
                // Fourth edge (west...)
                else if (iQuarter == 3)
                {
                    if (iEltTheta == _NTheta-1 && iLocY == _P)
                        return 0;
                    else
                        return (_NDoFX - (iLocY + 1 + iEltSquareEdge*_P))*_NDoFX;
                }
            }
            else if (iEltTheta == _NTheta-1 && iLocY == _P)
                return (iEltRho*_P - 1 + iLocX)*_NDoFTheta + _NDoFSquare;
            else
                return (iEltRho*_P - 1 + iLocX)*_NDoFTheta + (iEltTheta*_P + iLocY) + _NDoFSquare;

            assert(false && "Should have reached return.");
            return 0;
        }
        
 
        void _Rectangle2QuarterCircle(const RealVector &xyz, Index iQuarter, RealVector &uvw) const
        {
                //First Transfomration xyz -> bar xyz
                Real alpha = (sqrt(2)*exp(_L)-1.0)/_L;
                Real xBar1 = _L + (log((alpha*xyz[0]+1.0)/sqrt(2.0))-_L)/(log(sqrt(2.0))/_L+1.0);
                Real xBar2 = atan(xyz[1]);

                //Second Transfomration bar xyz -> hat xyz
                Real xHat1 = (xBar1 - _L) * ( log(sqrt(2) * cos(xBar2)) / _L + 1 ) + _L;
                Real xHat2 = xBar2;

                //Third Transformation hat xyz -> uvw
                uvw[0] = sqrt(2) * exp(xHat1) * cos(xHat2 + (M_PI/2) * iQuarter - M_PI/2);
                uvw[1] = sqrt(2) * exp(xHat1) * sin(xHat2 + (M_PI/2) * iQuarter - M_PI/2);
                uvw[2] = 0.0;
        }

                
        void _TemplateTransform (const Index &iEltGlob, const RealVector &p,
                                        RealVector &uvw) const
        {
            if (iEltGlob < _NSquareElem)
            {
                // Extracting quadrature coordinates in the template finite element space.
                Index iEltGlobX = iEltGlob % _NX;
                Index iEltGlobY = iEltGlob / _NX;

                // Adjusting quadrature point coordinate from the element index.
                uvw[0] = (p[0] + iEltGlobX) * _hX - 1;
                uvw[1] = (p[1] + iEltGlobY) * _hX - 1;
                uvw[2] = 0.0;
            }
            else
            {
                Index iEltRel = iEltGlob - _NSquareElem;
                Index iEltTheta = iEltRel % _NTheta;
                Index iEltRho = iEltRel / _NTheta;
                
                // Quadrature coordinates in the template finite element space.
                Index iQuarter = iEltTheta / _NX; //Between 0 and 3, number of the quarter circle
                Index iEltTilde1 = iEltRho;
                Index iEltTilde2 = iEltTheta % _NX;
                
                //space step in (xTilde1 , xTilde2)
                Real hTilde1 = _L / _NRho;
                Real hTilde2 = _hX;
                
                RealVector xTilde;
                //First Transfomration p->tilde : [0 , L]x[-1 , 1]
                xTilde[0] = (p[0] + iEltTilde1) * hTilde1;
                xTilde[1] = (p[1] + iEltTilde2) * hTilde2 - 1;
                xTilde[2] = 0.0;
                
                _Rectangle2QuarterCircle(xTilde, iQuarter, uvw);
            }

            uvw[0] = uvw[0] / 2.0;
            uvw[1] = uvw[1] / 2.0;
            //Circle of Unitary Radius!!!

        }

                /*! \brief Gives the coordinate of a point associated to a local DoF in an element.
            \param iEltGlob is a global index of an element.
            \param iLoc is a local index in the element iEltGlob.
            \param uvw are the coordinates of the point associated to iLoc which will be filled.
            \param F is the Jacobian matrix of the transform
        */
        void _TemplateTransformGradient(const Index &iEltGlob, const RealVector &p,
                                        RealVector &uvw, RealMatrix2x2 &F ) const
        {
            _TemplateTransform(iEltGlob,p,uvw);
            
            if (iEltGlob < _NSquareElem)
            {
                F[0][0] = _hX/2; F[0][1] = 0;
                F[1][0] = 0  ; F[1][1] = _hX/2;
            }
            else
            {
                Index iEltRel = iEltGlob - _NSquareElem;
                Index iEltTheta = iEltRel % _NTheta;
                Index iEltRho = iEltRel / _NTheta;

                Index iQuarter = (iEltTheta / _NX);
                Index iEltTilde1 = iEltRho;
                Index iEltTilde2 = iEltTheta % _NX;
                
                RealMatrix2x2 F_10;
                RealMatrix2x2 F_21;
                RealMatrix2x2 F_20;
                RealMatrix2x2 F_32;
                RealMatrix2x2 F_30;
                RealMatrix2x2 F_43;
                
                //space step in (xTilde1 , xTilde2)
                Real hTilde1 = _L / _NRho;
                Real hTilde2 = _hX;
                
                //First transformation p->tilde :  [0 , L]x[-1 , 1]
                Real xTilde1 = (p[0] + iEltTilde1) * hTilde1;
                Real xTilde2 = (p[1] + iEltTilde2) * hTilde2 - 1.0;
               
                F_10[0][0] = hTilde1; F_10[0][1] = 0.0;
                F_10[1][0] = 0.0;     F_10[1][1] = hTilde2;

                //Second transformation tilde->bar: [-pi/4 , pi/4]x[0 , L]
                Real alpha = (sqrt(2)*exp(_L)-1.0)/_L;
                Real xBar1 = _L + (log((alpha * xTilde1 + 1.0)/sqrt(2.0))-_L)/(log(sqrt(2.0))/_L+1.0);
                Real xBar2 = atan(xTilde2);
               
                F_21[0][0] = (alpha/(alpha * xTilde1 + 1.0))/(log(sqrt(2.0))/_L+1.0); F_21[0][1] = 0.0;
                F_21[1][0] = 0.0;                                                     F_21[1][1] = 1.0/(1.0 + (xTilde2)*(xTilde2));

                //Third transformation bar->hat: [-pi/4 , pi/4]x[log(1/sqrt(2)) , L]
                Real xHat1 = (xBar1 - _L) * ( (log(sqrt(2.0) * cos(xBar2)) / _L) + 1.0 ) + _L;
                Real xHat2 = xBar2;
             
                F_32[0][0] = (log(sqrt(2.0) * cos(xBar2)) / _L) + 1; F_32[0][1] = (xBar1 - _L) * (-1.0 / _L) * tan(xBar2);
                F_32[1][0] = 0.0;                                    F_32[1][1] = 1.0;
                
                // Fourth (conform) transformation
                F_43[0][0] = uvw[0]; F_43[0][1] = -uvw[1];
                F_43[1][0] = uvw[1]; F_43[1][1] = uvw[0];

                ArrayAlgebra::MatMlt(F_21 , F_10, F_20);
                ArrayAlgebra::MatMlt(F_32 , F_20, F_30);
                ArrayAlgebra::MatMlt(F_43 , F_30, F);
            }
        }


            // ---------------------------------------------------------------------------------//
            /*! \brief Definition of the dimension as a static variable.*/
            static const Index Dim = 2;
        
        
        Circle(Index P, Index N, GeoTransform aGeoTransform = IdentityMap2)
        {
            //Obvious initialisation
            _GeoTransform = aGeoTransform;
            _P = P;
            _GeomInfoPreComputed = false;
            
            //Internal Square
            // ---------------------------------------------------------------------------------//
            // Number of element in every directions and in total.
            _NX = 4*N;
            _NSquareElem = _NX * _NX;
            
            // Number of degrees of freedom (DoF) in the 2 directions and in total.
            _NDoFX = (_NX * _P + 1);
            _NDoFSquare = _NDoFX * _NDoFX;

            // Computing space step. 
            _hX = 2.0 / _NX; //NB: at the beginning the internal square is [-1,1]^2 - then I divide by 2 to get a unitary circle
            
            //External Circle
            // ---------------------------------------------------------------------------------//
            //The radius of the Circle is given by sqrt(2)*exp(L)
            _L = log(sqrt(2)); //Radius = 2; -> Be careful!!! In template transform I divide by 2 to get a unitary circle!
        
            // Number of element in every directions and in total.
            _NTheta = 4*_NX;
            _NRho = round(1.0 / _hX);
            
            _NCircleElem = _NRho * _NTheta;
            
            // Number of degrees of freedom (DoF) in the 2 directions and in total.
            _NDoFRho = (_NRho * _P);
            _NDoFTheta = (_NTheta * _P);
            
            // Computing space step.
            _hTheta = (2 * M_PI) / _NTheta;
            

            // ---------------------------------------------------------------------------------//
            //Total Number of Elements and of DoFs
            _NElem = _NSquareElem + _NCircleElem;
            _NDoF = (_NDoFX * _NDoFX) + (_NDoFRho * _NDoFTheta);

            // Number of local DoF.
            _NLocDoF = (_P + 1) * (_P + 1);

            // Number of elements per color.
            // Considering a restriction of N=Nx=Ny=Nrho and N%4=0 + the coloring that I have in mind
            _NElemColor = (_NX / 2 ) * (_NX / 2 ) + (_NTheta * _NRho / 4);


            // Initializing spectral finite element.
            _GLE = GaussLobattoElement(_P+1,_P+1,0);
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
            if (iColor<4) return _NElemColor;
    
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
            return _getNumVerticesIntSquare() + _getNumVerticesExtCircle();
        }

        Index getNumVerticesQ2() const
        {
            return _getNumVerticesIntSquareQ2() + _getNumVerticesExtCircleQ2();
        }

        bool isAffine(const Index iEltGlob) const
        {
            return false;
        }

        void getVertexCoordinate(const Index iVertex, RealVector& xyz) const
        {
            RealVector abc;
            RealVector uvw;
            
            Index NVertSquare = _getNumVerticesIntSquare();
            
            if (iVertex < NVertSquare)
            {
                Index iVertexX = iVertex % (_NX+1);
                Index iVertexY = iVertex / (_NX+1);
            
                uvw[0]=(iVertexX*_hX) - 1.0;
                uvw[1]=(iVertexY*_hX) - 1.0;
                uvw[2]=0;
            }
            else
            {
                
                Index iVertexTheta = (iVertex - NVertSquare) % _NTheta;
                Index iVertexRho   = (iVertex - NVertSquare) / _NTheta + 1;
                Index iQuarter     = iVertexTheta / _NX;

                Index iVertexX = iVertexTheta % _NX;
            
                abc[0]= iVertexRho* _L / _NRho;
                abc[1]= iVertexX* _hX -1.0;
                abc[2]=0;
                
                _Rectangle2QuarterCircle(abc, iQuarter, uvw);
                
            }

            uvw[0] = uvw[0] / 2.0;
            uvw[1] = uvw[1] / 2.0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }

        void getVertexCoordinateQ2(const Index iVertex, RealVector& xyz) const
        {
            
            RealVector abc;
            RealVector uvw;

            Index NVertQ2 = 2; 

            Index NVertSquareQ2 = _getNumVerticesIntSquareQ2();

            Real hXQ2 = _hX / NVertQ2;
            
            Index NVertThetaQ2 = _NTheta * NVertQ2;
            Index NVertRhoQ2 = _NRho * NVertQ2;
           
            if (iVertex < NVertSquareQ2)
            {
                Index iVertexX = iVertex % (_NX*NVertQ2 +1);
                Index iVertexY = iVertex / (_NX*NVertQ2 +1);
            
                uvw[0]=(iVertexX*hXQ2) - 1.0;
                uvw[1]=(iVertexY*hXQ2) - 1.0;
                uvw[2]=0;
            }
            else
            {
                
                Index iVertexTheta = (iVertex - NVertSquareQ2) % NVertThetaQ2;
                Index iVertexRho   = (iVertex - NVertSquareQ2) / NVertThetaQ2 + 1;
                Index iQuarter     = iVertexTheta / (_NX * NVertQ2 );

                Index iVertexX = iVertexTheta % (_NX * NVertQ2 );
            
                abc[0]= iVertexRho* _L / NVertRhoQ2;
                abc[1]= iVertexX* hXQ2 -1.0;
                abc[2]=0;
                
                _Rectangle2QuarterCircle(abc, iQuarter, uvw);
                
            }

            // Here as well it is to have a unitary circle
            uvw[0] = uvw[0] / 2.0; 
            uvw[1] = uvw[1] / 2.0;

            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
          
        }
         
        void getElementVertices(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            numbering.resize(4);
            Index NVertSquare = _getNumVerticesIntSquare();

            //Numbering of elements within the Central Square
            if (iEltGlob < _NSquareElem)
            {
                Index iEltGlobX = iEltGlob % _NX;
                Index iEltGlobY = iEltGlob / _NX;
                
                numbering[0] = (iEltGlobY + 0) * (_NX+1) + (iEltGlobX + 0);
                numbering[1] = (iEltGlobY + 0) * (_NX+1) + (iEltGlobX + 1);
                numbering[3] = (iEltGlobY + 1) * (_NX+1) + (iEltGlobX + 0);
                numbering[2] = (iEltGlobY + 1) * (_NX+1) + (iEltGlobX + 1);
                
            }

            //Critical Circular layer with internal border as a square
            else if ( (iEltGlob >= _NSquareElem)  &&  (iEltGlob < (_NSquareElem + _NTheta)) )
            {
                Index iEltRel = iEltGlob - _NSquareElem;
                Index iEltGlobTheta = iEltRel % _NTheta;
                Index temp = 0;

               
                numbering[1] = NVertSquare + (iEltGlobTheta + 0);
                numbering[2] = NVertSquare + (iEltGlobTheta + 1)%_NTheta;
                
                // NOTE: This SHOULD work ONLY IF _NX = _NX = N
                if ((iEltRel / _NX) == 0)
                {
                    numbering[0] = iEltGlobTheta;
                    numbering[3] = iEltGlobTheta + 1;
                }
                else if ((iEltRel / _NX) == 1)
                {
                    temp = iEltGlobTheta % _NX ;
                    numbering[0] = (_NX + 1)*(temp + 0) + _NX;
                    numbering[3] = (_NX + 1)*(temp + 1) + _NX;
                }
                else if ( (iEltRel / _NX) == 2 )
                {
                    temp = iEltGlobTheta % (2*_NX);
                    numbering[0] = NVertSquare - 1 - (temp + 0);
                    numbering[3] = NVertSquare - 1 - (temp + 1);
                }
                else if ( (iEltRel / _NX) == 3 )
                {
                    temp = iEltGlobTheta % (3*_NX);

                    numbering[0] = _NX*(_NX + 1) - (temp + 0)*(_NX + 1);
                    numbering[3] = _NX*(_NX + 1) - (temp + 1)*(_NX + 1);
                }
            }
            
            else if (iEltGlob >= (_NSquareElem + _NTheta) )
            {
                Index iEltRel = iEltGlob - _NSquareElem;
                Index iEltTheta = iEltRel % _NTheta;
                Index iEltRho = (iEltRel / _NTheta) - 1;

                numbering[0] = (iEltRho + 0) * _NTheta + (iEltTheta + 0)           + NVertSquare;
                numbering[3] = (iEltRho + 0) * _NTheta + (iEltTheta + 1) % _NTheta + NVertSquare;
                numbering[1] = (iEltRho + 1) * _NTheta + (iEltTheta + 0)           + NVertSquare;
                numbering[2] = (iEltRho + 1) * _NTheta + (iEltTheta + 1) % _NTheta + NVertSquare;
               
            }
        }

        void getElementVerticesQ2(const Index iEltGlob, std::vector<Index> & numbering) const
        {
            
            numbering.resize(9);
            Index NVertSquareQ2 = _getNumVerticesIntSquareQ2();
            Index NVertQ2 = 2;
            //Number of Vert in edge of internal square square 
            Index NVertXQ2 = (NVertQ2 * _NX) + 1;
            Index NVertThetaQ2 = _NTheta * NVertQ2 ;

            //Numbering of elements within the Central Square
            if (iEltGlob < _NSquareElem)
            {
                Index iEltGlobX = iEltGlob % _NX;
                Index iEltGlobY = iEltGlob / _NX;
                
                numbering[0] = (2*iEltGlobY + 0) * NVertXQ2 + (2*iEltGlobX + 0);
                numbering[4] = (2*iEltGlobY + 0) * NVertXQ2 + (2*iEltGlobX + 1);
                numbering[1] = (2*iEltGlobY + 0) * NVertXQ2 + (2*iEltGlobX + 2);
                numbering[7] = (2*iEltGlobY + 1) * NVertXQ2 + (2*iEltGlobX + 0);
                numbering[8] = (2*iEltGlobY + 1) * NVertXQ2 + (2*iEltGlobX + 1);
                numbering[5] = (2*iEltGlobY + 1) * NVertXQ2 + (2*iEltGlobX + 2);
                numbering[3] = (2*iEltGlobY + 2) * NVertXQ2 + (2*iEltGlobX + 0);
                numbering[6] = (2*iEltGlobY + 2) * NVertXQ2 + (2*iEltGlobX + 1);
                numbering[2] = (2*iEltGlobY + 2) * NVertXQ2 + (2*iEltGlobX + 2);
            }

            //Critical Circular layer with internal border as a square
            else if ( (iEltGlob >= _NSquareElem)  &&  (iEltGlob < (_NSquareElem + _NTheta)) )
            {
                Index iEltRel = iEltGlob - _NSquareElem;
                Index iEltGlobTheta = iEltGlob % _NTheta;
                Index temp = 0;
                
                numbering[4] = NVertSquareQ2 + (2*iEltGlobTheta + 0);
                numbering[8] = NVertSquareQ2 + (2*iEltGlobTheta + 1);
                numbering[6] = NVertSquareQ2 + (2*iEltGlobTheta + 2)%NVertThetaQ2;

                numbering[1] = NVertSquareQ2 + NVertThetaQ2 + (2*iEltGlobTheta + 0);
                numbering[5] = NVertSquareQ2 + NVertThetaQ2 + (2*iEltGlobTheta + 1);
                numbering[2] = NVertSquareQ2 + NVertThetaQ2 + (2*iEltGlobTheta + 2)%NVertThetaQ2;
                
                // NOTE: This SHOULD work ONLY IF _NX = _NX = N
                if ((iEltRel / _NX) == 0)
                {
                    numbering[0] = 2*iEltGlobTheta;
                    numbering[7] = 2*iEltGlobTheta + 1;
                    numbering[3] = 2*iEltGlobTheta + 2;
                }
                else if ((iEltRel / _NX) == 1)
                {
                    temp = iEltGlobTheta % _NX ;
                    numbering[0] = NVertXQ2 * (2*temp + 0) + (NVertXQ2 - 1);
                    numbering[7] = NVertXQ2 * (2*temp + 1) + (NVertXQ2 - 1);
                    numbering[3] = NVertXQ2 * (2*temp + 2) + (NVertXQ2 - 1);
                }
                else if ( (iEltRel / _NX) == 2 )
                {
                    temp = iEltGlobTheta % (_NX*2);
                    numbering[0] = NVertSquareQ2 - 1 - (temp*NVertQ2 + 0);
                    numbering[7] = NVertSquareQ2 - 1 - (temp*NVertQ2 + 1);
                    numbering[3] = NVertSquareQ2 - 1 - (temp*NVertQ2 + 2);
                }
                else if ( (iEltRel / _NX) == 3 )
                {
                    temp = iEltGlobTheta % (3*_NX);
                    numbering[0] = (NVertXQ2-1) * NVertXQ2 - (2*temp + 0)*(NVertXQ2);
                    numbering[7] = (NVertXQ2-1) * NVertXQ2 - (2*temp + 1)*(NVertXQ2);
                    numbering[3] = (NVertXQ2-1) * NVertXQ2 - (2*temp + 2)*(NVertXQ2);
                }
            }
            
            else if (iEltGlob >= (_NSquareElem + _NTheta) )
            {
                Index iEltRel = iEltGlob - _NSquareElem;
                Index iEltTheta = iEltRel % _NTheta;
                Index iEltRho = (iEltRel / _NTheta) -1 ;

                numbering[0] = (2*iEltRho + 1) * NVertThetaQ2 + (2*iEltTheta + 0)                + NVertSquareQ2;
                numbering[7] = (2*iEltRho + 1) * NVertThetaQ2 + (2*iEltTheta + 1)                + NVertSquareQ2;
                numbering[3] = (2*iEltRho + 1) * NVertThetaQ2 + (2*iEltTheta + 2) % NVertThetaQ2 + NVertSquareQ2;

                numbering[4] = (2*iEltRho + 2) * NVertThetaQ2 + (2*iEltTheta + 0)                + NVertSquareQ2;
                numbering[8] = (2*iEltRho + 2) * NVertThetaQ2 + (2*iEltTheta + 1)                + NVertSquareQ2;
                numbering[6] = (2*iEltRho + 2) * NVertThetaQ2 + (2*iEltTheta + 2) % NVertThetaQ2 + NVertSquareQ2;

                numbering[1] = (2*iEltRho + 3) * NVertThetaQ2 + (2*iEltTheta + 0)                + NVertSquareQ2;
                numbering[5] = (2*iEltRho + 3) * NVertThetaQ2 + (2*iEltTheta + 1)                + NVertSquareQ2;
                numbering[2] = (2*iEltRho + 3) * NVertThetaQ2 + (2*iEltTheta + 2) % NVertThetaQ2 + NVertSquareQ2;
            }
            
            
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


        

        // Function that maps the numbering of the element by color to a global numbering.
        // c is the index of the color.
        // e is the numbering of the element in the color c.
        // It returns the global numbering of the element.
        Index EltColorNum2EltGlobalNum(Index c, Index e) const
        {
            Index gex{0},gey{0};
            Index N_c_Square = _NX*_NX / 4 ; // Num of elements of color c in the square
            Index N_c_Quarter = _NX*_NRho / 4;
            
            if (e < N_c_Square)
            {
                Index ex = e%(_NX/2);
                Index ey = e/(_NX/2);

                switch (c)
                {
                    case 0: gex=2*ex;   gey=2*ey;   break;
                    case 1: gex=2*ex+1; gey=2*ey;   break;
                    case 2: gex=2*ex  ; gey=2*ey+1; break;
                    case 3: gex=2*ex+1; gey=2*ey+1; break;
                    default: assert(false);
                }
                return gey*_NX+gex;
            }

            else
            {
                Index eRel = e - N_c_Square;
                Index iQuarter = eRel / N_c_Quarter;
                Index eBar = eRel - (iQuarter * N_c_Quarter);
                Index ex = eBar%(_NX/2);
                Index ey = eBar/(_NX/2);

                if (iQuarter==0)
                {
                    switch (c)
                    {
                        case 0: gex=2*ex  ; gey=2*ey+1; break;
                        case 1: gex=2*ex+1; gey=2*ey+1; break;
                        case 2: gex=2*ex+1; gey=2*ey;   break;
                        case 3: gex=2*ex;   gey=2*ey;   break;
                        default: assert(false);
                    }
                }
                else if (iQuarter==1)
                {
                    switch (c)
                    {
                        case 0: gex=2*ex;   gey=2*ey;   break;
                        case 1: gex=2*ex+1; gey=2*ey+1; break;
                        case 2: gex=2*ex+1; gey=2*ey;   break;
                        case 3: gex=2*ex  ; gey=2*ey+1; break;
                        default: assert(false);
                    }
                }
                else if (iQuarter==2)
                {
                    switch (c)
                    {
                        case 0: gex=2*ex;   gey=2*ey;   break;
                        case 1: gex=2*ex+1; gey=2*ey;   break;
                        case 2: gex=2*ex+1; gey=2*ey+1;   break;
                        case 3: gex=2*ex  ; gey=2*ey+1; break;
                        default: assert(false);
                    }
                }
                else if (iQuarter==3)
                {
                    switch (c)
                    {
                        case 0: gex=2*ex;   gey=2*ey+1;   break;
                        case 1: gex=2*ex+1; gey=2*ey;   break;
                        case 2: gex=2*ex+1; gey=2*ey+1;   break;
                        case 3: gex=2*ex  ; gey=2*ey; break;
                        default: assert(false);
                    }
                }

                return gey*_NTheta + gex + iQuarter*_NX + _NSquareElem;
            }
        }
        

        Index _Loc2Glob(Index iEltGlob, Index iLoc) const
        {
            Index iLocX = iLoc % (_P + 1);
            Index iLocY = iLoc / (_P + 1);
                
            if (iEltGlob < _NSquareElem)
            {
                Index iEltGlobX = iEltGlob % _NX;
                Index iEltGlobY = iEltGlob / _NX;

                return (iEltGlobY*_P + iLocY) * _NDoFX + (iEltGlobX * _P + iLocX);

            }
            else return _Loc2GlobCircle(iEltGlob - _NSquareElem, iLocX, iLocY);
 
        }

       

    // ---------------------------------------------------------------------------------//
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
            const RealVector& p = _GLE.getPointCoordinate(iLoc);

            _TemplateTransform(iEltGlob, p, uvw);
 
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyz);
        }



      

        /*! \brief Gives a degree of freedom coordinate associated to a global DoF.
            \param iGlob is a global index of an DoF.
            \param xyz are the coordinates of the degree of freedom associated to iLoc which will be filled.
        */

        void getDoFCoordinate(const Index iGlob, RealVector& xyz) const
        {
            Index iLoc, iEltGlob;

            _Glob2Loc(iGlob, iLoc, iEltGlob);

            getDoFCoordinate(iEltGlob,iLoc,xyz);
        }
        

        /*! \brief Gives the barycenter of the element
            \param iEltGlob is a global index of an element.
            \param xyzCenter are the coordinates of the barycenter of the element
        */
        void getBarycenter(const Index iEltGlob, RealVector& xyzCenter) const
        {
            // Storing template Quadrature Coordinates.
            RealVector uvw;
            
            // Extracting corresponding quadrature point coordinate in the reference element.
            const RealVector& p = {0.5,0.5,0.0};

            _TemplateTransform(iEltGlob, p, uvw);
 
            // Applying transformation.
            _GeoTransform.Eval(uvw, xyzCenter);
        }
        
  
 
        // Get the label of the element (here 0)
        Index getLabel(const Index iEltGlob) const
        {
            return 0;
        }


        
        void _getGradDef(Index iEltGlob, Index iLoc,
			RealMatrix2x2& GradDef) const
        {
            assert(GeoTransform::isContinuouslyDifferentiable);

            RealVector uvw;
            RealMatrix2x2 GradDef_Tmp;
            RealMatrix2x2 GeoGradDef;
  
            const RealVector& p = _GLE.getPointCoordinate(iLoc);
    
            _TemplateTransformGradient(iEltGlob, p, uvw, GradDef_Tmp);
            _GeoTransform.getGradDef(uvw,GeoGradDef);

             ArrayAlgebra::MatMlt(GeoGradDef,GradDef_Tmp,GradDef);
        }


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
            return  _NX;
        }
        // ---------------------------------------------------------------------------------//
    };

 
}//OndoMathX
 
