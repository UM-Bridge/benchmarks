#pragma once

#include "ElementT.h"
#include "../ReferenceElement/RTetrahedron.h"

namespace OndoMathX {
    
        class Tetrahedron4 : public ElementT<class Tetrahedron4,RTetrahedron>
        {
            
        protected:
            
            std::array< std::shared_ptr<Point>, 4> _v;
            
            public :
            
            
            explicit Tetrahedron4(std::vector< std::shared_ptr<Point> > &v,Index label=0)
            : ElementT<class Tetrahedron4,RTetrahedron>(label)
            {
                for(Index i = 0; i < 4; i++) _v[i] = v[i];
            }
            
            Tetrahedron4() = delete;
            explicit Tetrahedron4(const Tetrahedron4& rhs) = default;
            explicit Tetrahedron4(Tetrahedron4&& rhs) = default;
            Tetrahedron4& operator=(const Tetrahedron4& rhs) = default;
            virtual ~Tetrahedron4();
            
            inline GeoElement getType() const { return GeoTetrahedron4; }
            
            inline Index getNumVertices() const { return 4; }
            inline const Point & getVertex(Index num) const
            {
                assert(num < 4);
                return *_v[num];
            }
/*
            inline std::shared_ptr<Point> getVertexPtr(Index num)
            {
                assert(num < 4);
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
                assert(num < 4);
                
                static std::array< std::function<void(const RealVector& pt, Real &s)> , 4> functions
                {
                    {
                        [](const RealVector& pt,Real &s) { s = 1. - pt[0] - pt[1] - pt[2] ;},
                        [](const RealVector& pt,Real &s) { s =      pt[0]                 ;},
                        [](const RealVector& pt,Real &s) { s =          	 pt[1]    	   ;},
                        [](const RealVector& pt,Real &s) { s =              		  pt[2];}
                    }
                };
                
                functions[num](uvw,r);
            }
            
            inline void getGradShapeFunction(Index num,const RealVector& /*uvw*/, RealVector &r) const
            {
                assert(num < 4);
                
                
                
                static std::array< std::function<void(RealVector&)> , 4> functions
                {
                    {
                        [](RealVector &s) {s[0] = -1.; s[1] = -1.; s[2] = -1.;},
                        [](RealVector &s) {s[0] =  1.; s[1] =  0.; s[2] =  0.;},
                        [](RealVector &s) {s[0] =  0.; s[1] =  1.; s[2] =  0.;},
                        [](RealVector &s) {s[0] =  0.; s[1] =  0.; s[2] =  1.;}
                    }
                };
                
                functions[num](r);
            }
            
            inline Edge getEdge(Index num) const
            {
                assert(num < 6);
                return Edge(_v[RTetrahedron::edges(num,0)],
                            _v[RTetrahedron::edges(num,1)]);
            }
            
            
            inline Face getFace(Index num) const
            {
                assert((num >= 0) && (num < 4));
                return Face(_v[RTetrahedron::faces(num,0)],
                            _v[RTetrahedron::faces(num,1)],
                            _v[RTetrahedron::faces(num,2)]);
            }
            
            inline Index getOrientationEdge(Index e) const
            {
                assert(e < 6);
                
                if ( _v[RTetrahedron::edges(e,0)]->getIndex()
                    >_v[RTetrahedron::edges(e,1)]->getIndex())
                    return 1;
                
                else return 0;
            }
            
            inline Index getOrientationFace(Index f) const
            {
                assert(f < 4);
                
                Index numOrientation = 0;
                std::array< std::shared_ptr<const Point>, 3> v;
                
                v[0]=_v[RTetrahedron::faces(f,0)];
                v[1]=_v[RTetrahedron::faces(f,1)];
                v[2]=_v[RTetrahedron::faces(f,2)];
                
                if (v[1]->getIndex()<v[numOrientation]->getIndex()) numOrientation = 1;
                if (v[2]->getIndex()<v[numOrientation]->getIndex()) numOrientation = 2;
                
                if (v[(numOrientation+2)%3]->getIndex() < v[(numOrientation+1)%3]->getIndex()) numOrientation += 3;
                
                return numOrientation;
            }
            
            inline bool isAffine() const
            {
                return true;
            }
            
        };
        
        inline Tetrahedron4::~Tetrahedron4() = default;
 
} // OndoMathX

 
