#pragma once

#include <memory>
#include <array>
#include <cassert>
#include "Point.h"

namespace OndoMathX {
            
    
            // A mesh Face.
            class Face
            {
            private:
                
                std::array< std::shared_ptr<const Point>, 4> _v;
                std::array<Index, 4> _si;
                
            public:
                
                //! \name Constructors and destructor.
                //@{
                
                // Default construcotr
                //explicit Face() {_v[0]=_v[1]=_v[2]=_v[3]=nullptr;_si[0]=_si[1]=_si[2]=_si[3]=0;};
                Face() = delete;
                
                //! From 4 vertices
                explicit Face(std::shared_ptr<const Point> v0,
                              std::shared_ptr<const Point> v1,
                              std::shared_ptr<const Point> v2,
                              std::shared_ptr<const Point> v3 = nullptr);
                
                //! Copy constructor.
                Face(const Face& rhs) = default;
                
                //! Move constructor.
                explicit Face(Face&& rhs) = default;
                
                //! Destructor.
                ~Face() = default;
                
                //@}
 
                //! Assignment
                Face& operator=(const Face& rhs) = default;
                
                inline Index getNumVertices() const { return (nullptr!=_v[3]) ? 4 : 3; }
                
                inline const Point& getVertex(Index i) const
                {
                    assert(i<getNumVertices());
                    assert(_v[i]!=nullptr);
                    return *_v[i];
                }
                
                inline const Point& getSortedVertex(Index i) const
                {
                    assert(i<getNumVertices());
                    assert(_v[_si[i]]!=nullptr);
                    return *_v[_si[i]];
                }
            };
            
            inline bool operator<(const Face &f1, const Face &f2)
            {
                
                if(f1.getSortedVertex(0).getIndex() < f2.getSortedVertex(0).getIndex()) return true;
                if(f1.getSortedVertex(0).getIndex() > f2.getSortedVertex(0).getIndex()) return false;
                if(f1.getSortedVertex(1).getIndex() < f2.getSortedVertex(1).getIndex()) return true;
                if(f1.getSortedVertex(1).getIndex() > f2.getSortedVertex(1).getIndex()) return false;
                if(f1.getSortedVertex(2).getIndex() < f2.getSortedVertex(2).getIndex()) return true;
                
                if(f1.getNumVertices()==4 && f2.getNumVertices()==4)
                {
                    if(f1.getSortedVertex(2).getIndex() > f2.getSortedVertex(2).getIndex()) return false;
                    if (f1.getSortedVertex(3).getIndex() < f2.getSortedVertex(3).getIndex()) return true;
                }
                else if (f1.getNumVertices()< f2.getNumVertices()) return true;
                
                return false;
            }
         
            
            
            inline Face::Face(std::shared_ptr<const Point> v0,std::shared_ptr<const Point> v1,
                       std::shared_ptr<const Point> v2,std::shared_ptr<const Point> v3)
            {
                _v[0] = v0; _v[1] = v1; _v[2] = v2; _v[3] = v3;

                // This is simply an unrolled insertion sort (hopefully fast).  Note that if
                // _v[3] == 0, _v[3] is not sorted.

                if(_v[1]->getIndex() < _v[0]->getIndex())
                {
                    _si[0] = 1;
                    _si[1] = 0;
                }
                else
                {
                    _si[0] = 0;
                    _si[1] = 1;
                }
                if(_v[2]->getIndex() < _v[_si[1]]->getIndex())
                {
                    _si[2] = _si[1];
                    
                    if(_v[2]->getIndex() < _v[_si[0]]->getIndex())
                    {
                        _si[1] = _si[0];
                        _si[0] = 2;
                    }
                    else _si[1] = 2;
                }
                else _si[2] = 2;
                
                if( _v[3]!=nullptr && (_v[3]->getIndex() < _v[_si[2]]->getIndex()) )
                {
                    _si[3] = _si[2];
                    
                    if(_v[3]->getIndex() < _v[_si[1]]->getIndex())
                    {
                        _si[2] = _si[1];
                        
                        if(_v[3]->getIndex() < _v[_si[0]]->getIndex())
                        {
                            _si[1] = _si[0];
                            _si[0] = 3;
                        } else _si[1] = 3;
                    }
                    else _si[2] = 3;
                }
                else _si[3] = 3;
            }


 
            inline bool getOrientation(const Face &f1, const Face &f2, Index & t)
            {
                assert(f1.getNumVertices()==f2.getNumVertices());
                
                Index nbVertices = f1.getNumVertices();
                
                //We assume that the first face is the reference
                if (f1.getVertex(0).getIndex()==f2.getVertex(0).getIndex())
                {
                    if (f1.getVertex(1).getIndex()==f2.getVertex(1).getIndex()) {t = 0;return true;}
                    else {t = nbVertices;return true;}
                }
                
                if (f1.getVertex(0).getIndex()==f2.getVertex(1).getIndex())
                {
                    if (f1.getVertex(1).getIndex()==f2.getVertex(2).getIndex()) {t = 1;return true;}
                    else {t = nbVertices+1;return true;}
                }
                    
                if (f1.getVertex(0).getIndex()==f2.getVertex(2).getIndex())
                {
                    if (f1.getVertex(1).getIndex()==f2.getVertex(3 % nbVertices).getIndex()) {t = 2;return true;}
                    else {t = nbVertices+2;return true;}
                }
                
                if (nbVertices==4 && (f1.getVertex(0).getIndex()==f2.getVertex(3).getIndex()))
                {
                    if (f1.getVertex(1).getIndex()==f2.getVertex(0).getIndex()) {t = 3;return true;}
                    else {t = nbVertices+3;return true;}
                }
                
                
                 t = 0;
                return false;
            }
            
        
 
        
    } //OndoMathX
  
