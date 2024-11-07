#pragma once

#include <cassert>

#include "../../Utility/Defines.h"

#include <vector>

#include "Segment2.h"
#include "Segment3.h"
#include "Triangle3.h"
#include "Quadrangle4.h"
#include "Tetrahedron4.h"
#include "Hexahedron8.h"
#include "Quadrangle9.h"
#include "Hexahedron27.h"

#include "../Topology/Point.h"


namespace OndoMathX {
    
 
        //!Just a switch that creates an element of the desired type
        inline std::shared_ptr<Element> Factory(GeoElement type,
                                         std::vector<std::shared_ptr<Point> > &v,
                                         Index label)
        {
            switch (type)
            {
                case GeoNone:         {assert(false);}
                case GeoPoint:        {assert(false);}
                case GeoTriangle6:    {assert(false);} break;
                case GeoSegment2:     return std::make_shared<Segment2>       (v,label);break;
                case GeoTriangle3:    return std::make_shared<Triangle3>      (v,label);break;
                case GeoQuadrangle4:  return std::make_shared<Quadrangle4>    (v,label);break;
                case GeoTetrahedron4: return std::make_shared<Tetrahedron4>   (v,label);break;
                case GeoHexahedron8:  return std::make_shared<Hexahedron8>    (v,label);break;
                case GeoSegment3:     return std::make_shared<Segment3>       (v,label);break;
                case GeoQuadrangle9:  return std::make_shared<Quadrangle9>    (v,label);break;
                case GeoHexahedron27: return std::make_shared<Hexahedron27>   (v,label);break;
            }
            
            return nullptr;
        }
        
    
    
} //OndoMathX


 
