#pragma once

#include "ElementT.h"
#include "../ReferenceElement/RSegment.h"
#include "Segment2.h"

namespace OndoMathX {
        
        /*
         * Segment3
         *
         *   0-----2-----1
         *
         *   0 -- x
         *
         */
        
        class Segment3 : public Segment2
        {
            
            std::shared_ptr<Point> _geoV;
            
        public:
            
            static const std::array< std::function<void(const Real&,Real &)> , 3> shapeFunction;
            
            static const std::array< std::function<void(const Real&,Real &)> , 3> dShapeFunction;
            
            
            explicit Segment3(std::vector< std::shared_ptr<Point> > &v, Index label=0)
            :  Segment2(v,label)
            {
                _geoV = v[2];
            }
            
            Segment3() = delete;
            explicit Segment3(const Segment3& rhs) = default;
            explicit Segment3(Segment3&& rhs) = default;
            Segment3& operator=(const Segment3& rhs) = default;
            virtual ~Segment3();
   
            inline GeoElement getType() const { return GeoSegment3; }
            
            inline Index getNumVertices() const { return 2; }
            inline const Point & getVertex(Index num) const
            {
                assert(num<2);
                return Segment2::getVertex(num);
            }
            
            inline Index getNumNodes() const { return 3; };
            inline const Point & getNode(Index num) const
            {
                assert(num<3);
                
                if (num<2) return Segment2::getVertex(num);
                else return *_geoV;
            };
    
            
            
            inline void getShapeFunction(Index num,const RealVector& uvw, Real &s) const
            {
                assert(num < 3);
                
                
                shapeFunction[num](uvw[0],s);
                
            }
            
            inline void getGradShapeFunction(Index num,const RealVector& uvw,  RealVector& s) const
            {
                assert(num < 3);
                
                s[1]=0.0;
                s[2]=0.0;
                
                dShapeFunction[num](uvw[0],s[0]);
            }
            
            inline bool isAffine() const
            {
                Real length = fabs(_v[1]->x()-_v[0]->x())
                +fabs(_v[1]->y()-_v[0]->y())
                +fabs(_v[1]->z()-_v[0]->z());
                
                Real ex = fabs(_geoV->x() -  0.5*(_v[0]->x()+_v[1]->x()));
                Real ey = fabs(_geoV->y() -  0.5*(_v[0]->y()+_v[1]->y()));
                Real ez = fabs(_geoV->z() -  0.5*(_v[0]->z()+_v[1]->z()));
                
                if (ex+ey+ez < REF_COORD_TOL*length)
                    return true;
                else
                    return false;
            } 
        };


inline Segment3::~Segment3() = default;

inline const std::array< std::function<void(const Real&, Real &)> , 3> Segment3::shapeFunction =
{
    {
        [](const Real& u,Real &s) { s=2*(0.5-u)*(1.0-u);},
        [](const Real& u,Real &s) { s=2*(u-0.5)*u;},
        [](const Real& u,Real &s) { s=4*u*(1.0-u);},
    }
};

inline const std::array< std::function<void(const Real&,Real &)> , 3> Segment3::dShapeFunction =
{
    {
        [](const Real& u,Real &s) {s=4*u-3.0;},
        [](const Real& u,Real &s) {s=4*u-1.0;},
        [](const Real& u,Real &s) {s=4-8.0*u;}
    }
};


} // OndoMathX
