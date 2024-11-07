#pragma once

#include "RElementT.h"


        namespace OndoMathX {
           
            class RSegment : public RElementT<RSegment>
            {
            public:
              
                
                static inline Index getNumFaces()                   {return 0;}
                static inline Index getNumEdges()                   {return 1;}
                static inline Index getNumVertices()                {return 2;}
                static inline Index getDim()                        {return 1;}
                static inline RefElement getRefType()    {return Segment;}
                static inline Index getNumInterfaces()              {return 2;}
                
                static inline RefElement getInterfaceRefType(int numInterface)
                {
                    assert((numInterface >= 0) && (numInterface < 2));
                    return RefPoint;
                }

                static inline RefElement getFaceRefType(Index /*numFace*/)
                {
                    return RefNone ;
                }
                
                static inline void getNormal(Index numInterface,RealVector &normal)
                {
                    assert((numInterface >= 0) && (numInterface < 2));
                    
                    normal[1]=0.0;
                    normal[2]=0.0;
                    switch(numInterface)
                    {
                        case 0  : normal[0]=-1.0; break;
                        case 1  : normal[0]=1.0; break;
                            
                        default : break;
                    }
                }
                
                static inline void getPosVertices(RealVector &pos,Index numRealVector)
                {
                    assert(numRealVector < 3);
                    
                    switch(numRealVector)
                    {
                        case 0: pos[0] = 0.0; break;
                        case 1: pos[0] = 1.0; break;
                    }
                }
                
 
              
                
                static inline bool isInside(const RealVector &pos)
                {
                    if (pos[2] < -REF_COORD_TOL || pos[2] > REF_COORD_TOL) return false;
                    if (pos[1] < -REF_COORD_TOL || pos[1] > REF_COORD_TOL) return false;
                    
                    if( pos[0] < 1.0+REF_COORD_TOL && pos[0] > -REF_COORD_TOL) return true;
                    return false;
                }
                
                
                static inline bool isOnEdge(const RealVector &pos,Index numInterface)
                {
                    assert(numInterface==0);
                    if (!isInside(pos)) return false;
                    return true;
                }
                
                static inline bool isOnVertex(const RealVector &pos,Index numInterface)
                {
                    if (!isInside(pos)) return false;
                    
                    if(numInterface == 0 && pos[0] < REF_COORD_TOL) return true;
                    if(numInterface == 1 && pos[0] > 1.0-REF_COORD_TOL) return true;
                    
                    return false;
                }
                
                static inline bool isOnFace(const RealVector &/*pos*/,Index /*numInterface*/)
                {
                    assert(false);
                    return false;
                }
                
                static inline Index edges(Index edge, Index vert)
                {
                    assert(edge < 1);
                    assert(vert == 1);
                    
                    return vert;
                }
                
                static inline Index faces(Index /*face*/, Index /*vert*/)
                {
                    assert(false);
                    return 0;
                }

                static inline void posEdgeTransformed(const RealVector &posRes,RealVector &pos,
                                                      Index numEdge, Index numTransfo)
                {
                    assert(isInside(posRes));
                    //assert(numTransfo<1);
                    assert(numEdge==0);
                    
                    pos = posRes;
                    
                    if (numTransfo == 1) pos[0] = 1-posRes[0];
                
                    
                    
                }
                
                static inline void posFaceTransformed(const RealVector &/*posRes*/,RealVector &/*pos*/,
                                                                  Index /*numFace*/, Index /*numTransfo*/)
                {
                    assert(false);
                }
       
            };
      
    } // OndoMathX
