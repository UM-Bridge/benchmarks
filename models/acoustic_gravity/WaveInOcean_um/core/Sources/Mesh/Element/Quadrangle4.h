#pragma once

#include "ElementT.h"
#include "../ReferenceElement/RQuadrangle.h"


namespace OndoMathX {

        class Quadrangle4: public ElementT<class Quadrangle4,RQuadrangle>
        {
        protected:
            
            std::array< std::shared_ptr<Point>, 4> _v;
            
        public:
            
            explicit Quadrangle4(std::vector< std::shared_ptr<Point> > &v,Index label=0)
            : ElementT<class Quadrangle4,RQuadrangle>(label)
            {
                for(Index i = 0; i < 4; i++) _v[i] = v[i];
            }
            
            Quadrangle4() = delete;
            explicit Quadrangle4(const Quadrangle4& rhs) = default;
            explicit Quadrangle4(Quadrangle4&& rhs) = default;
            Quadrangle4& operator=(const Quadrangle4& rhs) = default;
            virtual ~Quadrangle4();
            
            inline GeoElement getType() const { return GeoQuadrangle4; }
            
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
                        [](const RealVector& pt,Real &s) { s = (1.0 - pt[0]) * (1.0 - pt[1]);},
                        [](const RealVector& pt,Real &s) { s = pt[0] * (1.0 - pt[1])        ;},
                        [](const RealVector& pt,Real &s) { s = pt[0] * pt[1]                ;},
                        [](const RealVector& pt,Real &s) { s = (1.0 - pt[0]) * pt[1]        ;}
                    }
                };
                
                functions[num](uvw,r);
            }
            
            inline void getGradShapeFunction(Index num, const RealVector& uvw, RealVector &r) const
            {
                assert(num < 4);
                
                r[2]=0.0;
                
                static std::array< std::function<void(const RealVector&, RealVector&)> , 4> functions
                {
                    {
                        [](const RealVector& pt,std::array<Real, 3> &s) {s[0]=-(1.0-pt[1]); s[1]=-(1.0-pt[0]);},
                        [](const RealVector& pt,std::array<Real, 3> &s) {s[0]=(1.0 - pt[1]);s[1]=-pt[0];      },
                        [](const RealVector& pt,std::array<Real, 3> &s) {s[0]=pt[1]; 		 s[1]=pt[0]; 	   },
                        [](const RealVector& pt,std::array<Real, 3> &s) {s[0]=-pt[1];       s[1]=(1.0- pt[0]);}
                    }
                };
                
                functions[num](uvw,r);
            }
            
            
            
            inline Edge getEdge(Index num) const
            {
                assert(num < 4);
                return Edge(_v[RQuadrangle::edges(num,0)],
                            _v[RQuadrangle::edges(num,1)]);
            }
            
            
            
            inline Face getFace(Index num) const
            {
                assert(num==0);
                return Face(_v[0], _v[1], _v[2], _v[3]);
            }
            
            inline Index getOrientationEdge(Index e) const
            {
                assert(e < 4);
                if ( _v[RQuadrangle::edges(e,0)]->getIndex()
                    >_v[RQuadrangle::edges(e,1)]->getIndex())
                    return 1;
                else  return 0;
            }
            
            inline Index getOrientationFace(Index f) const
            {
                assert(f==0);
                
                Index numOrientation = 0;
                
                if (_v[1]->getIndex()<_v[numOrientation]->getIndex()) numOrientation = 1;
                if (_v[2]->getIndex()<_v[numOrientation]->getIndex()) numOrientation = 2;
                if (_v[3]->getIndex()<_v[numOrientation]->getIndex()) numOrientation = 3;
                
                if (_v[(numOrientation+3)%4]->getIndex() < _v[(numOrientation+1)%4]->getIndex()) numOrientation+=4;
                
                return numOrientation;
            }
            
            
            bool isAffine() const;
            
        };



inline Quadrangle4::~Quadrangle4() = default;


inline bool Quadrangle4::isAffine() const
{
    RealVector tmp = getBarycenter();
    
    Real e1x = 0.5*(_v[0]->x()+_v[2]->x())-tmp[0];
    Real e1y = 0.5*(_v[0]->y()+_v[2]->y())-tmp[1];
    Real e1z = 0.5*(_v[0]->z()+_v[2]->z())-tmp[2];
    
    Real e2x = 0.5*(_v[1]->x()+_v[3]->x())-tmp[0];
    Real e2y = 0.5*(_v[1]->y()+_v[3]->y())-tmp[1];
    Real e2z = 0.5*(_v[1]->z()+_v[3]->z())-tmp[2];
    
    Real length1 = fabs(_v[2]->x()-_v[0]->x())
    +fabs(_v[2]->y()-_v[0]->y())
    +fabs(_v[2]->z()-_v[0]->z());
    
    Real length2 = fabs(_v[3]->x()-_v[1]->x())
    +fabs(_v[3]->y()-_v[1]->y())
    +fabs(_v[3]->z()-_v[1]->z());
    
    if (e1x+e1y+e1z < REF_COORD_TOL*length1 && e2x+e2y+e2z < REF_COORD_TOL*length2)
        return true;
    else
        return false;
}

} // OndoMathX
