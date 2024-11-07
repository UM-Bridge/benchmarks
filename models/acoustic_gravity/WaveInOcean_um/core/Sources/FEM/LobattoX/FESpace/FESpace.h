#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

// -----------------------------------------------------------------------------------------//
// Including transformations finite element space.
#include "GeoTransform/Identity.h"
#include "GeoTransform/Translation.h"
#include "GeoTransform/Scaling.h"
#include "GeoTransform/Affine.h"
#include "GeoTransform/Diffeomorphism.h"
#include "GeoTransform/Homeomorphism.h"
// -----------------------------------------------------------------------------------------//



namespace OndoMathX
{
    inline Identity<3> IdentityMap3;
    inline Identity<2> IdentityMap2;
    inline Identity<1> IdentityMap1;


}

// -----------------------------------------------------------------------------------------//
// Generic finite element space.
#include "FESpaceT.h"
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
// Including 1D finite element space.
#include "Line/Line.h"
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
// Including 2D finite element space.
#include "Square/Square.h"
#include "Circle/Circle.h"
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
// Including 3D finite element space.
#include "Cube/Cube.h"
#include "Cylinder/Cylinder.h"
// -----------------------------------------------------------------------------------------//


// -----------------------------------------------------------------------------------------//
// Finite element space from q mesh.
#include "FEMesh.h"
// -----------------------------------------------------------------------------------------//


