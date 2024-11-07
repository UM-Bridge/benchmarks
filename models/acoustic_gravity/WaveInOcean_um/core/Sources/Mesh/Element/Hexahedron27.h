#pragma once

#include "ElementT.h"
#include "Hexahedron8.h"
#include "Segment3.h"

namespace OndoMathX {
 
        
        /*
         * Hexahedron27
         *
         *   3----13----2
         *   |\         |\
         *   |15    24  | 14
         *   9  \ 20    11 \
         *   |   7----19+---6
         *   |22 |  26  | 23|
         *   0---+-8----1   |
         *    \ 17    25 \  18
         *    10 |  21    12|
         *      \|         \|
         *       4----16----5
         *
         *   y
         *   |
         *   0 -- x
         *    \
         *     z
         */
        
        class Hexahedron27: public Hexahedron8
        {
        private:
            static constexpr Index _NVertexQ1 = 8;
            static constexpr Index _NNodeQ2 = 27;
            static constexpr Index _Node_Not_Vertex = _NNodeQ2 - _NVertexQ1;
            
            static const std::array< std::function<void(std::array<Index,3> &)> , _NNodeQ2> _tensorialNum;
            
            std::array<std::shared_ptr<Point>, (_Node_Not_Vertex) >_geoV;
            
        public:
            
            explicit Hexahedron27(std::vector< std::shared_ptr<Point> > &v, Index label=0)
            :  Hexahedron8(v,label)
            {
                // v contains all the nodes of Hexa27 (element Q2)
                // The first 8 nodes correspond to the vertex of Hexa27 
                // The vertex of Hexa27 are initiated with Hexa8 
                // The rest of nodes are put within _geoV 
                // _geoV = v[_NVertexQ1:end]
                
                for(Index i = 0; i < _Node_Not_Vertex; i++) 
                    _geoV[i] = v[i+_NVertexQ1];
            }
            
            Hexahedron27() = delete;
            explicit Hexahedron27(const Hexahedron27& rhs) = default;
            explicit Hexahedron27(Hexahedron27&& rhs) = default;
            Hexahedron27& operator=(const Hexahedron27& rhs) = default;
            virtual ~Hexahedron27();
            
            inline GeoElement getType() const { return GeoHexahedron27; }
            
            inline Index getNumVertices() const { return _NVertexQ1; }
            inline const Point & getVertex(Index num) const
            {
                assert(num<_NVertexQ1);
                return Hexahedron8::getVertex(num);
            }
            
            inline Index getNumNodes() const { return _NNodeQ2; }
            inline const Point & getNode(Index num) const
            {
                assert(num<_NNodeQ2);
                
                if (num<_NVertexQ1) return Hexahedron8::getVertex(num);
                else return *_geoV[num-_NVertexQ1];
            }
            
            
            
            
            
            inline void getShapeFunction(Index num,const RealVector& uvw, Real &s) const
            {
                Real sx,sy,sz;
                std::array<Index,3> numT;

                _tensorialNum[num](numT);
                Segment3::shapeFunction[numT[0]](uvw[0],sx);
                Segment3::shapeFunction[numT[1]](uvw[1],sy);
                Segment3::shapeFunction[numT[2]](uvw[2],sz);
                
                s = sx*sy*sz;
            }
            
            inline void getGradShapeFunction(Index num,const RealVector& uvw, RealVector &s) const
            {
                Real sx,sy,sz,dsx,dsy,dsz;
                std::array<Index,3> numT;
                
                _tensorialNum[num](numT);


                Segment3::dShapeFunction[numT[0]](uvw[0],dsx);
                Segment3::dShapeFunction[numT[1]](uvw[1],dsy);
                Segment3::dShapeFunction[numT[2]](uvw[2],dsz);
                
                Segment3::shapeFunction[numT[0]](uvw[0],sx);
                Segment3::shapeFunction[numT[1]](uvw[1],sy);
                Segment3::shapeFunction[numT[2]](uvw[2],sz);
                
                s[0]=dsx*sy*sz;
                s[1]=sx*dsy*sz;
                s[2]=sx*sy*dsz;

           
            }
            
            inline bool isAffine() const
            {
                return false;
            }
            
            
        };

inline Hexahedron27::~Hexahedron27() = default;

inline const std::array< std::function<void(std::array<Index,3> &)> , 27> Hexahedron27::_tensorialNum =
{
    /**/
    {
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=0;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=0;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=1;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=1;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=0;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=0;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=1;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=1;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=0;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=2;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=0;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=2;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=0;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=1;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=1;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=1;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=0;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=2;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=2;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=1;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=2;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=0;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=2;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=2;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=1;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=2;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=2;numT[2]=2;}
    }
    /*
    {
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=0;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=0;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=2;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=2;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=0;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=0;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=2;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=2;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=0;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=1;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=0;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=1;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=0;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=2;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=2;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=2;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=0;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=1;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=1;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=2;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=1;numT[2]=0;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=0;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=0;numT[1]=1;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=2;numT[1]=1;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=2;numT[2]=1;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=1;numT[2]=2;},
        [](std::array<Index,3> & numT){numT[0]=1;numT[1]=1;numT[2]=1;}
    }*/
};

} // OndoMathX
