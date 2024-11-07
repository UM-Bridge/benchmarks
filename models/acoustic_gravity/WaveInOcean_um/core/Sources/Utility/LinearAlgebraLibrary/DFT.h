#pragma once
 

// -----------------------------------------------------------------------------------------//
#include <cmath>
// -----------------------------------------------------------------------------------------//


// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

namespace LAL {

	
// DFT of a Complex vector 2D.
// Vector is the input Vector
// NDoFx*NDoFy is the dimension of the input Vector (rows/columns)

inline static void DFT2D(
                         ComplexVector& vector,
                         Index NDoFx,
                         Index NDoFy
                         )

{
    ComplexVector Temp;
    LAL::Allocate(Temp, NDoFx*NDoFy);
    
    Complex I(0,1)        ;
    Complex exp_value(0,0);
    Complex temp(0,0)     ;
 
    Real arg = 0;
    Index ind_res = 0;
    Index ind_sum = 0;
   
    //Recovers the data
    Complex * Data = getData(vector);
    Complex * DataTemp = getData(Temp);
    
    
    for (Index jDoF = 0; jDoF < NDoFy; jDoF++)
        for (Index iDoF = 0; iDoF < NDoFx; iDoF++)
        {
            ind_res = jDoF * NDoFx +iDoF;
            DataTemp[ind_res]= 0;
        }
    
    Real twoPi = 4 * asin(1);
    //DFT by column
    for (Index yDoF = 0; yDoF < NDoFy; yDoF++)
        for (Index tDoF = 0; tDoF < NDoFy; tDoF++)
        {
            arg = twoPi*yDoF*tDoF/NDoFy;
            exp_value = exp(-I * arg);
            for (Index xDoF = 0; xDoF < NDoFx; xDoF++)
            {
                temp = Data[tDoF * NDoFx + xDoF];
                ind_res = yDoF * NDoFx + xDoF;
                DataTemp[ind_res] += temp * exp_value;
            }
        }
    
    for (Index jDoF = 0; jDoF < NDoFy; jDoF++)
        for (Index iDoF = 0; iDoF < NDoFx; iDoF++)
        {
            ind_res = jDoF * NDoFx +iDoF;
            Data[ind_res]= 0;
        }
    
    //DFT by row
    for (Index xDoF = 0; xDoF < NDoFx; xDoF++)
        for (Index sDoF = 0; sDoF < NDoFx; sDoF++)
        {
            arg = twoPi*sDoF*xDoF/NDoFx;
            exp_value = exp(-I * arg);
            for (Index yDoF = 0; yDoF < NDoFy; yDoF++)
            {
                ind_sum = yDoF * NDoFx + sDoF;
                ind_res = yDoF * NDoFx + xDoF;
                temp = DataTemp[ind_sum];
                Data[ind_res] += temp * exp_value;
            }
        }
}




// DFT of a Complex vector 3D.
// Vector is the input Vector
// NDoFx*NDoFy*NDoFz is the dimension of the input Vector (rows/columns)

inline static void DFT3D(
                     ComplexVector& vector,
                     Index NDoFx,
                     Index NDoFy,
                     Index NDoFz
                     )
{
    
}



    } //LAL

} // OndoMathX
 
 
