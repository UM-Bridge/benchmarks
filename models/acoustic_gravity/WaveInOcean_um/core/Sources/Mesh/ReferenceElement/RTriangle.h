#pragma once

#include "RElementT.h"


        namespace OndoMathX {
             
            class RTriangle : public RElementT<RTriangle>
            {
            public:
                
                
                static inline Index getNumFaces()       {return 1;}
                static inline Index getNumEdges()       {return 3;}
                static inline Index getNumVertices() 	{return 3;}
                static inline Index getDim() 			{return 2;}
                static inline RefElement getRefType() {return  Triangle;}
                static inline Index getNumInterfaces()  {return 3;}
                
                
                static inline RefElement getInterfaceRefType(Index numInterface )
                {
                    assert(numInterface < 3);
                    return Segment;
                }
                
                static inline RefElement getFaceRefType(Index numFace)
                {
                    assert(numFace < 1);
                    return Triangle ;
                }
                
                static inline void getNormal(Index numInterface,RealVector &normal)
                {
                    assert((numInterface >= 0) && (numInterface < 3));
                    
                    const Real sqrt2Div2 = std::sqrt(2.)/2.;
                    
                    normal[2]=0.0;
                    switch(numInterface)
                    {
                        case 0  : normal[0]=0.0;normal[1]=-1.0;      break;
                        case 1  : normal[0]=sqrt2Div2;
                                  normal[1]=sqrt2Div2;               break;
                        case 2  : normal[0]=-1.0;normal[1]=0.0;      break;
                        default : break;
                    }
                }
                
                static inline void getPosVertices(RealVector &pos,Index numRealVector)
                {
                    assert((numRealVector >= 0) && (numRealVector < 3));
                    
                    switch(numRealVector)
                    {
                        case 0: pos[0] = 0.0; pos[1]= 0.0; break;
                        case 1: pos[0] = 1.0; pos[1]= 0.0; break;
                        case 2: pos[0] = 0.0; pos[1]= 1.0; break;
                    }
                }
                
                static inline bool isInside(const RealVector &pos)
                {
                    if (pos[2] < -REF_COORD_TOL || pos[2] > REF_COORD_TOL) return false;
                    
                    if( pos[0] > -REF_COORD_TOL && pos[1] > -REF_COORD_TOL
                       && pos[1] < 1.0-pos[0]+REF_COORD_TOL ) return true;
                    
                    return false;
                }
                
                static inline bool isOnEdge(const RealVector &pos,Index numInterface)
                {
                    if (!isInside(pos)) return false;
                    
                    if(numInterface == 0 && pos[1] < REF_COORD_TOL) return true;
                    if(numInterface == 1 && pos[1] > 1.0-pos[0]-REF_COORD_TOL) return true;
                    if(numInterface == 2 && pos[0] < REF_COORD_TOL) return true;
                    
                    return false;
                }
                
                static inline bool isOnVertex(const RealVector &pos,Index numInterface)
                {
                    if (!isInside(pos)) return false;
                    
                    if(numInterface == 0 && pos[1] < REF_COORD_TOL && pos[0] < REF_COORD_TOL) return true;
                    if(numInterface == 1 && pos[0] > 1.0-REF_COORD_TOL && pos[1] < REF_COORD_TOL) return true;
                    if(numInterface == 2 && pos[0] < REF_COORD_TOL && pos[1] > 1.0-REF_COORD_TOL) return true;
                    
                    return false;
                }
                
                static inline bool isOnFace(const RealVector &pos,Index numInterface)
                {
                    if (numInterface == 0 && isInside(pos)) return  true;
                    
                    return false;
                }
                
                static inline Index edges(Index edge, Index vert)
                {
                    assert((edge >= 0) && ((edge < 3)));
                    assert((vert == 0) || (vert == 1));
                    
                   const Index e[3][2] = {
                        {0, 1},
                        {1, 2},
                        {0, 2}
                    };
                    
                    return e[edge][vert];
                }
                
                static inline Index faces(Index face, Index vert)
                {
                    assert(face == 0);
                    assert(vert < 3);
                    
                    return vert;
                }
         
                static inline void posEdgeTransformed(const RealVector &posRes,RealVector &pos,
                                                      Index numEdge, Index numOrientation)
                {
                    assert(isOnEdge(posRes,numEdge));
                    assert(numOrientation<2);
                    
                    pos = posRes;
                    
                    if (numOrientation == 1)
                        switch(numEdge)
                    {
                        case 0  : pos[0]=1-pos[0];                  break;
                        case 1  : pos[0]=1-pos[0];pos[1]=1-pos[1];  break;
                        case 2  : pos[0]=1-pos[0];                  break;

                    }
                }
                
                static inline void posFaceTransformed(const RealVector &posRef1,RealVector &posRef2,
                                                      Index numFace, Index numOrientation)
                {
                    assert(isInside(posRef1));
                    assert(numFace==0);
                    assert((numOrientation < 6));
                    
                    posRef2[2]=0.0;
                    
                    switch(numOrientation)
                    {
                        case 0: posRef2[0] = posRef1[0];               posRef2[1] = posRef1[1]; 			break;
                        case 1: posRef2[0] = 1-posRef1[0]-posRef1[1];  posRef2[1] = posRef1[0]; 			break;
                        case 2: posRef2[0] = posRef1[1];               posRef2[1] = 1-posRef1[0]-posRef1[1];break;
                            
                        case 3: posRef2[0] = posRef1[1];               posRef2[1] = posRef1[0]; 			break;
                        case 4: posRef2[0] = 1-posRef1[0]-posRef1[1];  posRef2[1] = posRef1[1]; 			break;
                        case 5: posRef2[0] = posRef1[0];               posRef2[1] = 1-posRef1[0]-posRef1[1];break;
                    }
                    
                    
                    
                }
                
            };
            
 
    
}  //OndoMathX

 
