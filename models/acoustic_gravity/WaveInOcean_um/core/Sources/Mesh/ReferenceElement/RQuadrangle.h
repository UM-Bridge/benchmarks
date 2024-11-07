#pragma once

#include "RElementT.h"



        namespace OndoMathX {
 
            
            class RQuadrangle : public RElementT<RQuadrangle>
            {
            public:
                                
                static inline Index getNumFaces()      {return 1;}
                static inline Index getNumEdges()      {return 4;}
                static inline Index getNumVertices()   {return 4;}
                static inline Index getDim()           {return 2;}
                static inline RefElement getRefType() {return Quadrangle;}
                
                static inline RefElement getInterfaceRefType(int numInterface)
                {
                    assert((numInterface >= 0) && (numInterface < 4));
                    return Segment;
                }
  
                static inline RefElement getFaceRefType(Index numFace)
                {
                    assert(numFace < 1);
                    return Quadrangle ;
                }
                
                static inline void getNormal(Index numInterface,RealVector &normal)
                {
                    assert(numInterface < 4);
                    
                    normal[2]=0.0;
                    switch(numInterface)
                    {
                        case 0  : normal[0]=0.0;  normal[1]=-1.0;break;
                        case 1  : normal[0]=1.0;  normal[1]=0.0; break;
                        case 2  : normal[0]=0.0;  normal[1]=1.0; break;
                        case 3  : normal[0]=-1.0; normal[1]=0.0;break;
                        default : break;
                    }
                }
                
                static inline void getPosVertices(RealVector &pos,Index numRealVector)
                {
                    assert(numRealVector < 4);
                    
                    switch(numRealVector)
                    {
                        case 0: pos[0] = 0.0; pos[1]= 0.0; break;
                        case 1: pos[0] = 1.0; pos[1]= 0.0; break;
                        case 2: pos[0] = 1.0; pos[1]= 1.0; break;
                        case 3: pos[0] = 0.0; pos[1]= 1.0; break;
                    }
                }
                
                
                static inline bool isInside(const RealVector &pos)
                {
                    if (pos[2] < -REF_COORD_TOL || pos[2] > REF_COORD_TOL) return false;
                    
                    if( pos[0] < 1.0+REF_COORD_TOL && pos[0] > -REF_COORD_TOL
                       && pos[1] < 1.0+REF_COORD_TOL && pos[1] > -REF_COORD_TOL ) return true;
                    
                    return false;
                }
                
                static inline bool isOnEdge(const RealVector &pos,
                                     Index numInterface)
                {
                    if (!isInside(pos)) return false;
                    
                    if(numInterface == 0 && pos[1] < REF_COORD_TOL) return true;
                    if(numInterface == 1 && pos[0] > 1.0-REF_COORD_TOL) return true;
                    if(numInterface == 2 && pos[1] > 1.0-REF_COORD_TOL) return true;
                    if(numInterface == 3 && pos[0] < REF_COORD_TOL) return true;
                    
                    return false;
                }
                
                static inline bool isOnVertex(const RealVector &pos,
                                       Index numInterface)
                {
                    if (!isInside(pos)) return false;
                    
                    if(numInterface == 0 && pos[1] < REF_COORD_TOL && pos[0] < REF_COORD_TOL)
                        return true;
                    if(numInterface == 1 && pos[0] > 1.0-REF_COORD_TOL && pos[1] < REF_COORD_TOL)
                        return true;
                    if(numInterface == 2 && pos[1] > 1.0-REF_COORD_TOL && pos[0] > 1.0-REF_COORD_TOL)
                        return true;
                    if(numInterface == 3 && pos[0] < REF_COORD_TOL && pos[1] > 1.0-REF_COORD_TOL)
                        return true;
                    
                    return false;
                }
                
                static inline bool isOnFace(const RealVector &pos,
                                     Index numInterface)
                {
                    if (numInterface == 0 && isInside(pos)) return true;
                    
                    return false;
                }
                
                static inline Index edges(Index edge, Index vert)
                {
                    assert(edge < 4);
                    assert((vert == 0) || (vert == 1));
                    
                    static const Index e[4][2] = {
                        {0, 1},
                        {1, 2},
                        {3, 2},
                        {0, 3}
                    };
                    return e[edge][vert];
                }
                
                static inline Index faces(Index face, Index vert)
                {
                    assert(face == 0);
                    assert(vert < 4);
                    
                    return vert;
                }
                
                
                static inline void posEdgeTransformed(const RealVector &posRef,RealVector &pos,
                                                      Index numEdge, Index numOrientation)
                {
                    assert(isOnEdge(posRef,numEdge));
                    assert(numOrientation<2);
                    
                    pos = posRef;
                    
                    if (numOrientation == 1)
                        switch(numEdge)
                        {
                            case 0  : pos[0]=1-pos[0];  break;
                            case 1  : pos[1]=1-pos[1];  break;
                            case 2  : pos[0]=1-pos[0];  break;
                            case 3  : pos[1]=1-pos[1];  break;
                        }
                }

                static inline void posFaceTransformed(const RealVector &posRef1,RealVector &posRef2,
                                                      Index numFace, Index numOrientation)
                {
                    assert(isInside(posRef1));
                    assert(numFace==0);
                    assert((numOrientation <8));
                    
                    posRef2[2]=0.0;
  
                    switch(numOrientation)
                    {
                            
                        case 0: posRef2[0] = posRef1[0];      posRef2[1] = posRef1[1];    break;
                            
                        case 1: posRef2[0]= 1.0-posRef1[1];   posRef2[1] = posRef1[0]; 	  break;
                            
                        case 2: posRef2[0] = 1.0-posRef1[0];  posRef2[1] = 1.0-posRef1[1];break;
                            
                        case 3: posRef2[0] = posRef1[1];      posRef2[1] = 1.0-posRef1[0];break;
                            
                        case 4: posRef2[0] = posRef1[1];      posRef2[1] = posRef1[0];    break;
                    
                        case 5: posRef2[0] = 1.0-posRef1[0];  posRef2[1] = posRef1[1]; 	  break;
                            
                        case 6: posRef2[0] = 1.0-posRef1[1];  posRef2[1] = 1-posRef1[0];  break;
                            
                        case 7: posRef2[0] = posRef1[0];      posRef2[1] = 1- posRef1[1]; break;
                    }
                }
                
            };
       
        
    } // OndoMathX
     
