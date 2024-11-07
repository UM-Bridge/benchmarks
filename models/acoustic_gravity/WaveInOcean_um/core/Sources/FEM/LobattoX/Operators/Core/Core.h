#pragma once


namespace OndoMathX {
 
namespace Core {

enum DiffOp
{
    Zero,
    Identity,
    Gradient,
    Divergence,
    Trace,
    NormalTrace,
    TraceNormal
};

enum QuadratureOpt
{
    Default,
    Trial,
    Other
};

 

template<Index DimA, Index DimB> auto InitBuffer()
{
    if constexpr(DimA == 1 && DimB == 1)
    {
        Real t;
        return t;
    }
    else if constexpr(DimA == 1)
    {
        std::array<Real, DimB> t;
        return t;
    }
    else if constexpr(DimB == 1)
    {
        std::array<Real, DimA> t;
        return t;
    }
    else
    {
        std::array<std::array<Real,DimB>, DimA> t;
        return t;
    }

}


class DefaultFE
{

};

DefaultFE _DefaultFE;



} //Core
 
}//OndoMathX
