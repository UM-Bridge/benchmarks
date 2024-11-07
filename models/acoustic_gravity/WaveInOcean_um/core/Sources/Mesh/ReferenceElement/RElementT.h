#pragma once

#include <cassert>
#include "../../Utility/Defines.h"
#include "../Topology/Point.h"



        namespace OndoMathX {
            
            
            template<class DerivedT>
            class RElementT
            {
            public:
                
                //!Returns the number of interfaces of the element
                /*!
                 *If element dimension is 3, returns the number of faces,
                 *if element dimension is 2, returns the number of edges,
                 *if element dimension is 1, returns 2.
                 * @return the number of interfaces of the element
                 */
                static inline Index getNumInterface()
                {
                    if (DerivedT::getDim()==2) return DerivedT::getNumEdges();
                    else if (DerivedT::getDim()==1) return 2;
                    else if (DerivedT::getDim()==0) return 0;
                    else return DerivedT::getNumFaces();
                }
                
                //!Tests if a position is on a given interface of the element
                /**
                 * @param  pos a position
                 * @return true if \a pos belongs the interface \a numInterface, false otherwise.
                 * @see isInside
                 * @see isOnEdge
                 * @see isOnCoord
                 * @see isOnFace
                 */
                static inline bool isOnInterface(const RealVector &pos,int numInterface=0)
                {
                    if (DerivedT::getDim()==2) return  DerivedT::isOnEdge(pos,numInterface);
                    else return DerivedT::isOnFace(pos,numInterface);
                }
                
                //!Tests if a position is on one of the vertices of the element
                /*!
                 * @param  pos a position
                 * @return true if \a pos is one of the vertices of the element, false otherwise.
                 */
                static inline bool isOnVertices(const RealVector &pos)
                {
                    bool result = false;
                    for(int i=0;i<DerivedT::getNumVertices();++i)
                        result=result|DerivedT::isOnVertex(pos,i);
                    return result;
                }
                
         
                
                //!Tests if position is on one of the edges of the element
                /*!
                 * @param  pos a position
                 * @return true if \a pos is on one of the edges of the element, false otherwise.
                 */
                static inline bool isOnEdges(const RealVector &pos)
                {
                    bool result = false;
                    for(int i=0;i<DerivedT::getNumEdges();++i)
                        result=result|DerivedT::isOnEdge(pos,i);
                    return result;
                }
                
                //!Tests if position is on one of the faces of the element
                /*!
                 * @param pos a position
                 * @return true if \a pos is on one of the faces of the element, false otherwise.
                 */
                static inline bool isOnFaces(const RealVector &pos)
                {
                    bool result = false;
                    for(int i=0;i<DerivedT::getNumFaces();++i)
                        result=result|DerivedT::isOnFace(pos,i);
                    return result;
                }
            };
             
        
    } // OndoMathX

