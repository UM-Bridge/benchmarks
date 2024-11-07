#pragma once

#include <memory>
#include <array>
#include <cassert>
#include "Point.h"

namespace OndoMathX {
             
            // A mesh edge.
            class Edge
            {
            private:
                
                std::array< std::shared_ptr<const Point>, 2> _v;
                std::array<Index, 2> _si;
                
            public:
                
                //! \name Constructors and destructor.
                //@{
                
                // Default construcotr  
                //explicit Edge() {_v[0]=_v[1]=nullptr;_si[0]=_si[1]=0;};
                Edge() = delete;
                
                //! From two vertices
                explicit Edge(std::shared_ptr<const Point> v0, std::shared_ptr<const Point> v1)
                {
                    _v[0] = v0;
                    _v[1] = v1;
                    
                    if(_v[1]->getIndex() < _v[0]->getIndex() )
                    {
                        _si[0] = 1;
                        _si[1] = 0;
                    }
                    else
                    {
                        _si[0] = 0;
                        _si[1] = 1;
                    }
                }
                
                //! Copy constructor.
                Edge(const Edge& rhs) = default;
                
                //! Move constructor.
                explicit Edge(Edge&& rhs) = default;
                
                //! Destructor.
                ~Edge() = default;
            

                //@}
                
                //! Assignment
                Edge& operator=(const Edge& rhs) = default;
                
                inline Index getNumVertices() const {return 2;}
                
                inline const Point& getVertex(Index i)        const
                {assert(i<2);assert(_v[i]!=nullptr); return *_v[i];}
                
                inline const Point& getSortedVertex(Index i)  const
                {assert(i<2);assert(_v[_si[i]]!=nullptr); return *_v[_si[i]];}
            };
            
            inline bool operator<(const Edge &e1, const Edge &e2)
            {
                if(e1.getSortedVertex(0).getIndex()  < e2.getSortedVertex(0).getIndex() ) return true;
                if(e1.getSortedVertex(0).getIndex()  > e2.getSortedVertex(0).getIndex() ) return false;
                if(e1.getSortedVertex(1).getIndex()  < e2.getSortedVertex(1).getIndex() ) return true;
                
                return false;
            }
            
            inline bool getOrientation(const Edge &e1, const Edge &e2, Index & t)
            {
                if (e1.getVertex(0).getIndex() == e2.getVertex(0).getIndex()
                 && e1.getVertex(1).getIndex() == e2.getVertex(1).getIndex())
                {
                    t = 0;
                    return true;
                }
 
                if (e1.getVertex(0).getIndex() == e2.getVertex(1).getIndex()
                    && e1.getVertex(1).getIndex() == e2.getVertex(0).getIndex())
                {
                    t = 1;
                    return true;
                }
                
                t=0;
                return false;
            }
            
        
        
    }  //OndoMathX

 
