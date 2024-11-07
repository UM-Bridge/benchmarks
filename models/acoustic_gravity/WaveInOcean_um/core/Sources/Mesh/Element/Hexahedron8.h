#pragma once

#include "ElementT.h"
#include "../ReferenceElement/RHexahedron.h"

namespace OndoMathX {
        
        class Hexahedron8 : public ElementT<class Hexahedron8,RHexahedron>
        {
        private:

            static constexpr Index _NVertexQ1 = 8;

        protected:

            std::array< std::shared_ptr<Point>, 8> _v;
            
            public :
            
            explicit Hexahedron8(std::vector< std::shared_ptr<Point> > &v, Index label=0)
            : ElementT<class Hexahedron8,RHexahedron>(label)
            {
                for(Index i = 0; i < _NVertexQ1; i++)
                     _v[i] = v[i];
            }
            
            
            Hexahedron8() = delete;
            explicit Hexahedron8(const Hexahedron8& rhs) = default;
            explicit Hexahedron8(Hexahedron8&& rhs) = default;
            Hexahedron8& operator=(const Hexahedron8& rhs) = default;
            virtual ~Hexahedron8();
            
            
            inline GeoElement getType() const { return GeoHexahedron8; }
            
            inline Index getNumVertices() const { return _NVertexQ1; }
            
            inline const Point & getVertex(Index num) const
            {
                assert(num < _NVertexQ1);
                return *_v[num];
            }
/*
            inline std::shared_ptr<Point> getVertexPtr(Index num)
            {
                assert(num < _NVertexQ1);
                return _v[num];
            }
     */       
            inline Index getNumNodes() const { return getNumVertices(); };
            inline const Point & getNode(Index num) const
            {
                return getVertex(num);
            };
            
            
            
            inline void getShapeFunction(Index num,const RealVector& uvw, Real &r) const
            {
                
                assert(num < _NVertexQ1);
                
                static std::array< std::function<void(const RealVector& pt, Real &s)> , 8> functions
                {
                    {
                        [](const RealVector& pt,Real &s) { s = (1.0 - pt[0]) * (1.0 - pt[1]) * (1. - pt[2] );},
                        [](const RealVector& pt,Real &s) { s = pt[0]         * (1. - pt[1])  * (1. - pt[2] );},
                        [](const RealVector& pt,Real &s) { s = pt[0]         * pt[1]         * (1. - pt[2] );},
                        [](const RealVector& pt,Real &s) { s = (1. - pt[0])  * pt[1]         * (1. - pt[2] );},
                        [](const RealVector& pt,Real &s) { s = (1. - pt[0])  * (1. - pt[1])  * pt[2];},
                        [](const RealVector& pt,Real &s) { s = pt[0] 	      * (1. - pt[1])  * pt[2];},
                        [](const RealVector& pt,Real &s) { s = pt[0] 	      * pt[1] 		  * pt[2];},
                        [](const RealVector& pt,Real &s) { s = (1. - pt[0])  * pt[1] 	 	  * pt[2];}
                    }
                };
                
                functions[num](uvw,r);
                
            }
            
            inline void getGradShapeFunction(Index num,const RealVector& uvw, RealVector &r) const
            {
                assert(num < _NVertexQ1);
                
                
                static std::array< std::function<void(const RealVector&, RealVector&)> , _NVertexQ1> functions
                {
                    {
                        [](const RealVector& pt,RealVector &s)
                        {   s[0]= - (1. - pt[1]) * (1. - pt[2]);
                            s[1] = - (1. - pt[0]) * (1. - pt[2]);
                            s[2] = - (1. - pt[0]) * (1. - pt[1]);},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = (1. - pt[1]) * (1. - pt[2]);
                            s[1] = - pt[0] * (1. - pt[2]);
                            s[2] = - pt[0] * (1. - pt[1]);},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = pt[1] * (1. - pt[2]);
                            s[1] = pt[0] * (1. - pt[2]);
                            s[2] = - pt[0] * pt[1];},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = - pt[1] * (1. - pt[2]);
                            s[1] = (1. - pt[0]) * (1. - pt[2]);
                            s[2] = - (1. - pt[0]) * pt[1];},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = - (1. - pt[1]) * pt[2];
                            s[1] = - (1. - pt[0]) * pt[2];
                            s[2] = (1. - pt[0]) * (1. - pt[1]);},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = (1. - pt[1]) * pt[2];
                            s[1] = - pt[0] * pt[2];
                            s[2] = pt[0] * (1. - pt[1]);},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = pt[1] * pt[2];
                            s[1] = pt[0] * pt[2];
                            s[2] = pt[0] * pt[1];},
                        [](const RealVector& pt,RealVector &s)
                        {   s[0] = - pt[1] * pt[2];
                            s[1] = (1. - pt[0]) * pt[2];
                            s[2] = (1. - pt[0]) * pt[1];},
                    }
                };
                
                functions[num](uvw,r);
                
            }
            
            inline Edge getEdge(Index num) const
            {
                assert(num < 12);
                return Edge(_v[RHexahedron::edges(num,0)],
                            _v[RHexahedron::edges(num,1)]);
            }
            
            
            
            inline Face getFace(Index num) const
            {
                assert(num < 6);
                
                return Face(_v[RHexahedron::faces(num,0)],
                            _v[RHexahedron::faces(num,1)],
                            _v[RHexahedron::faces(num,2)],
                            _v[RHexahedron::faces(num,3)]);
            }
            
            inline Index getOrientationEdge(Index e) const
            {
                assert(e < 12);
                
                if ( _v[RHexahedron::edges(e,0)]->getIndex()
                    >_v[RHexahedron::edges(e,1)]->getIndex())
                    return 1;
                else  return 0;
            }
            
            inline Index getOrientationFace(Index f) const
            {
                assert(f < 6);
                
                Index numOrientation = 0;
                
                std::array< std::shared_ptr<const Point>, 4> v;
                
                v[0]=_v[RHexahedron::faces(f,0)];
                v[1]=_v[RHexahedron::faces(f,1)];
                v[2]=_v[RHexahedron::faces(f,2)];
                v[3]=_v[RHexahedron::faces(f,3)];
                
                if (v[1]->getIndex()<v[numOrientation]->getIndex()) numOrientation = 1;
                if (v[2]->getIndex()<v[numOrientation]->getIndex()) numOrientation = 2;
                if (v[3]->getIndex()<v[numOrientation]->getIndex()) numOrientation = 3;
                
                if (v[(numOrientation+3)%4]->getIndex() < v[(numOrientation+1)%4]->getIndex()) numOrientation+=4;
                
                return numOrientation;
            }
            
            inline bool isAffine() const
            {
                return false;
            }
            
            //virtual void xyz2uvw(const Real *xyz, Real *uvw) const;
            
        };

        inline Hexahedron8::~Hexahedron8() = default;

} // OndoMathX
