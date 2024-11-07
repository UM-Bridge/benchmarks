#pragma once

namespace OndoMathX {

    namespace Field {

        enum Regularity
        {
            Constant,
            PiecewiseConstant,
            C0,
            PiecewiseC0,
            PointData
        };

        class Identity { public : static const Field::Regularity Regularity = Constant;};
        inline Identity IdentityField ;
       
    }
}

#include "Scalar.h"
#include "Vector.h"
#include "Matrix.h"
#include "ThirdOrderTensor.h"
#include "IsotropicTensor.h"
#include "FourthOrderTensor.h"
#include "NeoHookeanTensor.h"

