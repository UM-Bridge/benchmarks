#pragma once

// -----------------------------------------------------------------------------------------//
#include <vector>
#include <cmath>
#include <cassert>
#include <array>

#include "GaussLobattoElement.h"
 
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
namespace OndoMathX
{


class GaussLobattoInterpolator {

private:
   
   std::vector<std::vector<Real>>  _ValuesI;
   std::vector<std::vector<Real>>  _ValuesJ;
   std::vector<std::vector<Real>>  _ValuesK;

   Index _Dim;
   
   Index _NI;
   Index _NJ;
   Index _NK;
   
   Index _MI;
   Index _MJ;
   Index _MK;
   
public:

   GaussLobattoInterpolator() {};
   
   //Interpolation from N_GL to M_GL
   GaussLobattoInterpolator(
        const GaussLobattoElement &N_GL,
        const GaussLobattoElement &M_GL
       )
   {
        _NI = N_GL.getNI();
        _NJ = N_GL.getNJ();
        _NK = N_GL.getNK();

        _MI = M_GL.getNI();
        _MJ = M_GL.getNJ();
        _MK = M_GL.getNK();
        
        assert(N_GL.getDim()==M_GL.getDim());
       
        _Dim = N_GL.getDim();
   
   
        std::vector<Real> PointsNI = N_GL.getPointsI();
        std::vector<Real> PointsNJ = N_GL.getPointsJ();
        std::vector<Real> PointsNK = N_GL.getPointsK();
    
        std::vector<Real> PointsMI = M_GL.getPointsI();
        std::vector<Real> PointsMJ = M_GL.getPointsJ();
        std::vector<Real> PointsMK = M_GL.getPointsK();
       
        _ValuesI.resize(_NI);
        for (Index i1 = 0; i1 < _NI; ++i1)
            _ValuesI[i1].resize(_MI);
        
        _ValuesJ.resize(_NJ);
        for (Index j1 = 0; j1 < _NJ; ++j1)
            _ValuesJ[j1].resize(_MJ);
        
        _ValuesK.resize(_NK);
        for (Index k1 = 0; k1 < _NK; ++k1)
            _ValuesK[k1].resize(_MK);
    
        if (_NI > 0 && _MI > 0)
        {
       
            for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index i2 = 0; i2 < _MI; ++i2)
                    _ValuesI[i1][i2] = _Interpolation(PointsNI, _NI, PointsMI[i2], i1);
       
        }
         
       
        if (_NJ > 0 && _MJ > 0)
        {
       
            for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index j2 = 0; j2 < _MJ; ++j2)
                    _ValuesJ[j1][j2] = _Interpolation(PointsNJ, _NJ, PointsMJ[j2], j1);
       
        }
         
       
        if (_NK > 0 && _MK > 0)
        {
            for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index k2 = 0; k2 < _MK; ++k2)
                    _ValuesK[k1][k2] = _Interpolation(PointsNK, _NK, PointsMK[k2], k1);
        }
         
    }
    
    Index getNumPoints() const
    {
        switch(_Dim)
        {
            case 1: return _NI; break;
            case 2: return _NI*_NJ; break;
            case 3: return _NI*_NJ*_NK; break;
        }
        
        return 0;
    }
    
    Index getNumPointsInterpolation() const
    {
        switch(_Dim)
        {
            case 1: return _MI; break;
            case 2: return _MI*_MJ; break;
            case 3: return _MI*_MJ*_MK; break;
        }
        
        return 0;
    }
    
    
    template<Index Dim=1, bool erase_output = true> void Interpolation(const std::vector<double>& Vec, std::vector<double>& InterpolationVec) const
    {

        if (_Dim == 1)
        {
            assert(InterpolationVec.size() == Dim*_MI);
            assert(Vec.size() == Dim*_NI);

            if (erase_output == true)
                for (Index i = 0; i < Dim*_MI; ++i) InterpolationVec[i] = 0.0;

            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index i2 = 0; i2 < _MI; ++i2)
                    for (Index d = 0 ; d<Dim; ++d)
                        InterpolationVec[i2*Dim+d] += Vec[i1*Dim+d] * _ValuesI[i1][i2];
        }
        if (_Dim == 2)
        {
            assert(InterpolationVec.size() == Dim*_MI*_MJ);
            assert(Vec.size() == Dim*_NI*_NJ);

            std::vector<double> TmpInterpolationVecI(Dim*_MI*_NJ,0.0);

            if (erase_output == true)
                for (Index i = 0; i < Dim*_MI*_MJ; ++i) InterpolationVec[i] = 0.0;

            // First dimension.
            for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index i1 = 0; i1 < _NI; ++i1)
                    for (Index i2 = 0; i2 < _MI; ++i2)
                        for (Index d = 0 ; d<Dim; ++d) TmpInterpolationVecI[(i2 + _MI*j1)*Dim+d] += Vec[(i1 + _NI*j1)*Dim+d] * _ValuesI[i1][i2];


            // Second dimension.
            for (Index d = 0 ; d<Dim; ++d)
                for (Index i1 = 0; i1 < _MI; ++i1)
                    for (Index j1 = 0; j1 < _NJ; ++j1)
                        for (Index j2 = 0; j2 < _MJ; ++j2)
                            for (Index d = 0 ; d<Dim; ++d) InterpolationVec[(i1 + _MI*j2)*Dim+d] += TmpInterpolationVecI[(i1 + _MI*j1)*Dim+d] * _ValuesJ[j1][j2];

        }
        else if (_Dim == 3)
        {
            assert(InterpolationVec.size() == Dim*_MI*_MJ*_MK);
            assert(Vec.size() == Dim*_NI*_NJ*_NK);

            std::vector<double> TmpInterpolationVecI(Dim*_MI*_NJ*_NK,0.0);
            std::vector<double> TmpInterpolationVecJ(Dim*_MI*_MJ*_NK,0.0);

            if (erase_output == true)
                for (Index i = 0; i < Dim*_MI*_MJ*_MK; ++i) InterpolationVec[i] = 0.0;

            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                    for (Index k1 = 0; k1 < _NK; ++k1)
                        for (Index i2 = 0; i2 < _MI; ++i2)
                            for (Index d = 0 ; d<Dim; ++d) TmpInterpolationVecI[(i2 + _MI*j1 + _MI*_NJ*k1)*Dim+d] += Vec[(i1 + _NI*j1 + _NI*_NJ*k1)*Dim+d] * _ValuesI[i1][i2];

            // Second dimension.
            for (Index i1 = 0; i1 < _MI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                    for (Index k1 = 0; k1 < _NK; ++k1)
                        for (Index j2 = 0; j2 < _MJ; ++j2)
                            for (Index d = 0 ; d<Dim; ++d) TmpInterpolationVecJ[(i1 + _MI*j2 + _MI*_MJ*k1)*Dim+d] += TmpInterpolationVecI[(i1 + _MI*j1 + _MI*_NJ*k1)*Dim+d] * _ValuesJ[j1][j2];

            // Third dimension.
            for (Index i1 = 0; i1 < _MI; ++i1)
                for (Index j1 = 0; j1 < _MJ; ++j1)
                    for (Index k1 = 0; k1 < _NK; ++k1)
                        for (Index k2 = 0; k2 < _MK; ++k2)
                            for (Index d = 0 ; d<Dim; ++d) InterpolationVec[(i1 + _MI*j1 + _MI*_MJ*k2)*Dim+d] += TmpInterpolationVecJ[(i1 + _MI*j1 + _MI*_MJ*k1)*Dim+d] * _ValuesK[k1][k2];
        }
    }


    template<Index Dim=1, bool erase_output = true> void TransposeInterpolation(const std::vector<double>& InterpolationVec, std::vector<double>& Vec) const
    {


        if (_Dim == 1)
        {
            assert(InterpolationVec.size() == Dim*_MI);
            assert(Vec.size() == Dim*_NI);

            if (erase_output == true)
                for (Index i = 0; i < Dim*_NI; ++i) Vec[i] = 0.0;

            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index i2 = 0; i2 < _MI; ++i2)
                    for (Index d = 0 ; d<Dim; ++d) Vec[i1*Dim+d] += InterpolationVec[i2*Dim+d] * _ValuesI[i1][i2];
        }
        if (_Dim == 2)
        {
            assert(InterpolationVec.size() == Dim*_MI*_MJ);
            assert(Vec.size() == Dim*_NI*_NJ);

            std::vector<double> TmpInterpolationVecI(Dim*_MI*_NJ,0.0);

            if (erase_output == true)
                for (Index i = 0; i < Dim*_NI*_NJ; ++i) Vec[i] = 0.0;

            // Second dimension.
            for (Index i1 = 0; i1 < _MI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                    for (Index j2 = 0; j2 < _MJ; ++j2)
                        for (Index d = 0 ; d<Dim; ++d) TmpInterpolationVecI[(i1 + _MI*j1)*Dim+d] += InterpolationVec[(i1 + _MI*j2)*Dim+d]  * _ValuesJ[j1][j2];


            // First dimension.
            for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index i1 = 0; i1 < _NI; ++i1)
                    for (Index i2 = 0; i2 < _MI; ++i2)
                        for (Index d = 0 ; d<Dim; ++d) Vec[(i1 + _NI*j1)*Dim+d] += TmpInterpolationVecI[(i2 + _MI*j1)*Dim+d] * _ValuesI[i1][i2];

        }
        else if (_Dim == 3)
        {
            assert(InterpolationVec.size() == Dim*_MI*_MJ*_MK);
            assert(Vec.size() == Dim*_NI*_NJ*_NK);

            std::vector<double> TmpInterpolationVecI(Dim*_MI*_NJ*_NK,0.0);
            std::vector<double> TmpInterpolationVecJ(Dim*_MI*_MJ*_NK,0.0);

            if (erase_output == true)
                for (Index i = 0; i < Dim*_NI*_NJ*_NK; ++i) Vec[i] = 0.0;

            // Third dimension.
            for (Index i1 = 0; i1 < _MI; ++i1)
                for (Index j1 = 0; j1 < _MJ; ++j1)
                    for (Index k1 = 0; k1 < _NK; ++k1)
                        for (Index k2 = 0; k2 < _MK; ++k2)
                            for (Index d = 0 ; d<Dim; ++d) TmpInterpolationVecJ[(i1 + _MI*j1 + _MI*_MJ*k1)*Dim+d] += InterpolationVec[(i1 + _MI*j1 + _MI*_MJ*k2)*Dim+d] * _ValuesK[k1][k2];

            // Second dimension.
            for (Index i1 = 0; i1 < _MI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                    for (Index k1 = 0; k1 < _NK; ++k1)
                        for (Index j2 = 0; j2 < _MJ; ++j2)
                            for (Index d = 0 ; d<Dim; ++d) TmpInterpolationVecI[(i1 + _MI*j1 + _MI*_NJ*k1)*Dim+d] += TmpInterpolationVecJ[(i1 + _MI*j2 + _MI*_MJ*k1)*Dim+d]  * _ValuesJ[j1][j2];

            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                    for (Index k1 = 0; k1 < _NK; ++k1)
                        for (Index i2 = 0; i2 < _MI; ++i2)
                            for (Index d = 0 ; d<Dim; ++d) Vec[(i1 + _NI*j1 + _NI*_NJ*k1)*Dim+d] += TmpInterpolationVecI[(i2 + _MI*j1 + _MI*_NJ*k1)*Dim+d]  * _ValuesI[i1][i2];
        }
    }
};

  
GaussLobattoInterpolator Interpolator(const GaussLobattoElement &N_GL,
                                       const GaussLobattoElement &M_GL)
{
    return GaussLobattoInterpolator(N_GL,M_GL);
}
 
 

} // OndoMathX





