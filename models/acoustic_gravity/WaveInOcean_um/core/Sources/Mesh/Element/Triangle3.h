#pragma once

#include "ElementT.h"
#include "../ReferenceElement/RTriangle.h"


namespace OndoMathX {
    
        /************************************************************************
         *  Triangle3
         *
         *     2
         *     | \
         *     |   \
         *     |     \
         *     |       \
         *     |         \
         *     0-----------1
         *
         *************************************************************************/
        
        class Triangle3: public ElementT<class Triangle3,RTriangle>
        {
            
            
        protected:
            
            std::array< std::shared_ptr<Point>, 3> _v;
            
        public:
            
            explicit Triangle3(std::vector< std::shared_ptr<Point> > &v,Index label=0)
            : ElementT<class Triangle3,RTriangle>(label)
            {
                for(Index i = 0; i < 3; i++) _v[i] = v[i];
            }
            
            Triangle3() = delete;
            explicit Triangle3(const Triangle3& rhs) = default;
            explicit Triangle3(Triangle3&& rhs) = default;
            Triangle3& operator=(const Triangle3& rhs) = default;
            virtual ~Triangle3();
            
            inline GeoElement getType() const { return GeoTriangle3; }
            
            inline Index getNumVertices() const { return 3; }
            inline const Point & getVertex(Index num) const
            {
                assert(num < 3);
                return *_v[num];
            }
            /*
            inline std::shared_ptr<Point> getVertexPtr(Index num)
            {
                assert(num < 3);
                return _v[num];
            }*/
            
            inline Index getNumNodes() const { return getNumVertices(); };
            inline const Point & getNode(Index num) const
            {
                return getVertex(num);
            };
     
            
            inline void getShapeFunction(Index num,const RealVector& uvw, Real &r) const
            {
                assert(num < 3);
                
                static std::array< std::function<void(const RealVector& pt, Real &s)> , 3> functions
                {
                    {
                        [](const RealVector& pt,Real &s) { s = (1.0 - pt[0] - pt[1]);},
                        [](const RealVector& pt,Real &s) { s = pt[0]                ;},
                        [](const RealVector& pt,Real &s) { s = pt[1]                ;},
                    }
                };
                
                functions[num](uvw,r);
                
            }
            
            inline void getGradShapeFunction(Index num, const RealVector& /*uvw*/, RealVector &r) const
            {
                assert(num < 3);
                
                r[2]=0.0;
                
                static std::array< std::function<void(RealVector&)> , 3> functions
                {
                    {
                        [](RealVector &s) {s[0]=-1.0;s[1]=-1.0;},
                        [](RealVector &s) {s[0]=1.0; s[1]=0.0; },
                        [](RealVector &s) {s[0]=0.0; s[1]=1.0; }
                    }
                };
                
                functions[num](r);
                
            }
            
            
            inline Edge getEdge(Index num) const
            {
                assert(num < 3);
                
                return Edge(_v[RTriangle::edges(num,0)],
                            _v[RTriangle::edges(num,1)]);
            }
            
            inline Face getFace(Index num) const
            {
                assert(num==0);
                
                return Face(_v[0], _v[1], _v[2]);
            }
            
            
            inline Index getOrientationEdge(Index e) const
            {
                assert(e < 3);
                if ( _v[RTriangle::edges(e,0)]->getIndex()
                    >_v[RTriangle::edges(e,1)]->getIndex())
                    return 1;
                else return 0;
            }
            
            
            inline Index getOrientationFace(Index f) const
            {
                assert(f==0);
                
                Index numOrientation = 0;
                
                if (_v[1]->getIndex()<_v[numOrientation]->getIndex()) numOrientation = 1;
                if (_v[2]->getIndex()<_v[numOrientation]->getIndex()) numOrientation = 2;
                
                if (_v[(numOrientation+2)%3]->getIndex() < _v[(numOrientation+1)%3]->getIndex()) numOrientation += 3;
                
                return numOrientation;
            }
            
            
            inline bool isAffine() const
            {
                return true;
            }
            
            
            bool xyz2uvw(const RealVector & xyz,  RealVector & uvw) const
            {
                if (getDim()==2)
                {
                    std::array<Real, 2> xyz_tmp = {xyz[0]-(*_v[0])[0], xyz[1]-(*_v[0])[1]};
                    std::array<Real, 2> uv;
                    
                    RealVector zero = {0, 0 ,0};
                    RealMatrix2x2 graddef,inv;
                    
                    getGradDef(zero, graddef);
                    ArrayAlgebra::Inv(graddef, inv);
                    ArrayAlgebra::MatMlt(inv, xyz_tmp, uv);
                    
                    uvw[0]=uv[0];
                    uvw[1]=uv[1];
                    uvw[2]=0.0;
                    
                    return true;
                }
                else
                {
                    return ElementT<class Triangle3,RTriangle>::xyz2uvw(xyz, uvw);
                }
            }
          
        };


        inline Triangle3::~Triangle3() = default;

} // OndoMathX
