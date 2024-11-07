#pragma once

#include "ElementT.h"
#include "Quadrangle4.h"
#include "Segment3.h"

namespace OndoMathX {
        
        /*
         * Quadrangle9
         *
         *   3-----6-----2
         *   |           |
         *   |           |
         *   7     8     5
         *   |           |
         *   |           |
         *   0-----4-----1
         *
         *
         *   y
         *   |
         *   0 -- x
         *
         *
         */
        
        
        class Quadrangle9: public Quadrangle4
        {
        private:
            
            static const std::array<std::function<void(std::array<Index,2> &)> , 9> _tensorialNum;
            std::array<std::shared_ptr<const Point>, (9-4) >_geoV;
            
        public:
            
            explicit Quadrangle9(std::vector< std::shared_ptr<Point> > &v, Index label=0)
            :  Quadrangle4(v,label)
            {
                for(Index i = 0; i < 5; i++) _geoV[i] = v[i+4];
            }
            
            Quadrangle9() = delete;
            explicit Quadrangle9(const Quadrangle9& rhs) = default;
            explicit Quadrangle9(Quadrangle9&& rhs) = default;
            Quadrangle9& operator=(const Quadrangle9& rhs) = default;
            virtual ~Quadrangle9();
            
            inline GeoElement getType() const { return GeoQuadrangle9; }
            
            inline Index getNumVertices() const { return 4; }
            inline const Point & getVertex(Index num) const
            {
                assert(num<4);
                
                return Quadrangle4::getVertex(num);
            }
            
            inline Index getNumNodes() const { return 9; };
            inline const Point & getNode(Index num) const
            {
                assert(num<9);
                
                if (num<4) return Quadrangle4::getVertex(num);
                else return *_geoV[num-4];
            };
            
            
            inline void getShapeFunction(Index num,const RealVector& uvw, Real &s) const
            {
                Real sx,sy;
                std::array<Index,2> numT;
                
                _tensorialNum[num](numT);
                Segment3::shapeFunction[numT[0]](uvw[0],sx);
                Segment3::shapeFunction[numT[1]](uvw[1],sy);
                s = sx*sy;
            }
            
            inline void getGradShapeFunction(Index num,const RealVector& uvw, RealVector &s) const
            {
                Real sx,sy,dsx,dsy;
                std::array<Index,2> numT;
                
                _tensorialNum[num](numT);
                
                Segment3::dShapeFunction[numT[0]](uvw[0],dsx);
                Segment3::dShapeFunction[numT[1]](uvw[1],dsy);
                
                Segment3::shapeFunction [numT[0]](uvw[0],sx);
                Segment3::shapeFunction [numT[1]](uvw[1],sy);
                
                s[0]=dsx*sy;
                s[1]=sx*dsy;
                s[2]=0.0;
            }
            
            inline bool isAffine() const
            {
                return false;
            }
            
        };
    

inline Quadrangle9::~Quadrangle9() = default;

inline const std::array< std::function<void(std::array<Index,2> &)> , 9> Quadrangle9::_tensorialNum =
{
    {
        [](std::array<Index,2> & numT){numT[0]=0;numT[1]=0;},
        [](std::array<Index,2> & numT){numT[0]=1;numT[1]=0;},
        [](std::array<Index,2> & numT){numT[0]=1;numT[1]=1;},
        [](std::array<Index,2> & numT){numT[0]=0;numT[1]=1;},
        [](std::array<Index,2> & numT){numT[0]=2;numT[1]=0;},
        [](std::array<Index,2> & numT){numT[0]=1;numT[1]=2;},
        [](std::array<Index,2> & numT){numT[0]=2;numT[1]=1;},
        [](std::array<Index,2> & numT){numT[0]=0;numT[1]=2;},
        [](std::array<Index,2> & numT){numT[0]=2;numT[1]=2;}
    }
};

} // OndoMathX
