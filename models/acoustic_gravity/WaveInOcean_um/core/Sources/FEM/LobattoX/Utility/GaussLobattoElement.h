#pragma once

// -----------------------------------------------------------------------------------------//
#include <vector>
#include <cmath>
#include <cassert>
#include <array>

#include "GaussLobattoElement.hxx"
// -----------------------------------------------------------------------------------------//

// -----------------------------------------------------------------------------------------//
namespace OndoMathX
{


class GaussLobattoElement
{
    
    
private:
    
    Index _NI;
    Index _NJ;
    Index _NK;
    Index _Dim;
    
    /*! \brief Quadrature weights of a spectral element.*/
    std::vector<Real> _WeightsList;
    std::vector<Real> _WeightsListI;
    std::vector<Real> _WeightsListJ;
    std::vector<Real> _WeightsListK;
    
    /*! \brief Points of a spectral element.*/
    std::vector<RealVector> _PointsList;
    std::vector<Real> _PointsI;
    std::vector<Real> _PointsJ;
    std::vector<Real> _PointsK;

    std::vector<Real> _InterpolationWeightI;
    std::vector<Real> _InterpolationWeightJ;
    std::vector<Real> _InterpolationWeightK;
    
    /*! \brief Projected gradient of basis function in first direction.*/
    std::vector<std::vector<Real>>  _DerivativeValuesI;
  
    /*! \brief Projected gradient of basis function in second direction.*/
    std::vector<std::vector<Real>>  _DerivativeValuesJ;
    
    /*! \brief Projected gradient of basis function in third direction.*/
    std::vector<std::vector<Real>>  _DerivativeValuesK;
    
public:
    
    GaussLobattoElement() {};
    
    
    GaussLobattoElement(Index NI, Index NJ, Index NK) : _NI(NI), _NJ(NJ),_NK(NK) 
    {
        assert(NI>0);
        if (NJ==0) {_Dim=1; assert(NK==0);}
        else if (NK==0) _Dim=2;
        else _Dim=3;
        
        bool Lobatto = true;
        
        switch(_Dim)
        {
            case 1:
                _ComputeGaussFormulas(NI, _PointsI, _WeightsList, Lobatto);
                _PointsList.resize(NI);
                
                for (Index i = 0; i < NI; ++i)
            {
                _PointsList[i][0] = _PointsI[i];
                _PointsList[i][1] = 0.0;
                _PointsList[i][2] = 0.0;
            }
                break;
            case 2:
                _ComputeGaussFormulas(NI, _PointsI, _WeightsListI, Lobatto);
                _ComputeGaussFormulas(NJ, _PointsJ, _WeightsListJ, Lobatto);
                
                _WeightsList.resize(NI*NJ);
                _PointsList.resize(NI*NJ);
                
                for (Index j = 0; j < NJ; ++j)
                for (Index i = 0; i < NI; ++i)
            {
                Index index =   j*NI + i;
                
                _WeightsList[index] = _WeightsListI[i] * _WeightsListJ[j];
                
                _PointsList[index][0] = _PointsI[i];
                _PointsList[index][1] = _PointsJ[j];
                _PointsList[index][2] = 0.0;
            }
                break;
            case 3:
                _ComputeGaussFormulas(NI, _PointsI, _WeightsListI, Lobatto);
                _ComputeGaussFormulas(NJ, _PointsJ, _WeightsListJ, Lobatto);
                _ComputeGaussFormulas(NK, _PointsK, _WeightsListK, Lobatto);
                
                _WeightsList.resize(NI*NJ*NK);
                _PointsList.resize(NI*NJ*NK);
                
                for (Index k = 0; k < NK; ++k)
                for (Index j = 0; j < NJ; ++j)
                for (Index i = 0; i < NI; ++i)
            {
                Index index = k*NI*NJ + j*NI + i;
                
                _WeightsList[index] = _WeightsListI[i] * _WeightsListJ[j] * _WeightsListK[k];
                
                _PointsList[index][0] = _PointsI[i];
                _PointsList[index][1] = _PointsJ[j];
                _PointsList[index][2] = _PointsK[k];
            }
                break;
        }
        
        
        
        if (NI > 0)
        {
            _DerivativeValuesI.resize(NI);
            _InterpolationWeightI.resize(NI);

            for (Index i1 = 0; i1 < NI; ++i1)
            {
                _DerivativeValuesI[i1].resize(NI);

                _InterpolationWeightI[i1] = _InterpolationWeight(_PointsI, NI, i1);
                
                for (Index i2 = 0; i2 < NI; ++i2)
                {
                    Real value =  _DerivativeInterpolation(_PointsI, NI, _PointsI[i2], i1);

                    _DerivativeValuesI[i1][i2] = value;
                }
            }

            
        }
        
        if (NJ > 0)
        {
            _DerivativeValuesJ.resize(NJ);
            _InterpolationWeightJ.resize(NJ);
           
            for (Index j1 = 0; j1 < NJ; ++j1)
            {
                _DerivativeValuesJ[j1].resize(NJ);

                _InterpolationWeightJ[j1] = _InterpolationWeight(_PointsJ, NJ, j1);
                
                for (Index j2 = 0; j2 < (NJ); ++j2)
                _DerivativeValuesJ[j1][j2] = _DerivativeInterpolation(_PointsJ, NJ, _PointsJ[j2], j1);
            }
        }
        
        if (NK > 0)
        {
            _DerivativeValuesK.resize(NK);
            _InterpolationWeightK.resize(NK);
            
            for (Index k1 = 0; k1 < (NK); ++k1)
            {
                _DerivativeValuesK[k1].resize(NK);

                _InterpolationWeightK[k1] = _InterpolationWeight(_PointsK, NK, k1);
                
                for (Index k2 = 0; k2 < (NK); ++k2)
                _DerivativeValuesK[k1][k2] = _DerivativeInterpolation(_PointsK, NK, _PointsK[k2], k1);
            }
        }



        


    }
    // ---------------------------------------------------------------------------------//
    
    // ---------------------------------------------------------------------------------//
    
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
    
    Index getNI() const {return _NI;}
    Index getNJ() const {return _NJ;}
    Index getNK() const {return _NK;}
   
    Index getDim() const {return _Dim;}
    
    Real getQuadratureWeightI(Index iLoc) const
    {
        return _WeightsListI[iLoc];
    }

    Real getQuadratureWeightJ(Index iLoc) const
    {
        return _WeightsListJ[iLoc];
    }

    Real getQuadratureWeightK(Index iLoc) const
    {
        return _WeightsListK[iLoc];
    }
    
    Real getQuadratureWeight(Index iLoc) const
    {
        return _WeightsList[iLoc];
    }

    Real getPointI(Index iLoc) const
    {
        return _PointsI[iLoc];
    }

    Real getPointJ(Index iLoc) const
    {
        return _PointsJ[iLoc];
    }

    Real getPointK(Index iLoc) const
    {
        return _PointsK[iLoc];
    }

    const std::vector<Real> & getPointsI() const
    {
        return _PointsI;
    }

    const std::vector<Real> & getPointsJ() const
    {
        return _PointsJ;
    }

    const std::vector<Real> & getPointsK() const
    {
        return _PointsK;
    }

    
    /*! \brief Gives a quadrature point associated to a local DoF.
     \param iLoc is a local index.
     \return The quadrature point associated to iLoc.
     */
    const RealVector & getPointCoordinate(Index iLoc) const
    {
        return _PointsList[iLoc];
    }
    // ---------------------------------------------------------------------------------//
    
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Return the value of the basis function number n
     \param n is a local index.
     \param uvw are the coordinates of the evaluation point.
     \return The value of the basis function.
     */
    
    Real getValueBasisFunc(Index n, const RealVector & uvw) const
    {
        Real res;
        
        if (_Dim == 1)
        {
            Index i = n % _NI;
            res = _Interpolation(_PointsI, _NI, uvw[0], i);
        }
        else if (_Dim == 2)
        {
            Index i = n % _NI;
            Index j = n / _NJ;
            res =  _Interpolation(_PointsI, _NI, uvw[0], i);
            res *= _Interpolation(_PointsJ, _NJ, uvw[1], j);
        }
        else
        {
            Index m = n % (_NI*_NJ);
            Index i = m % _NI;
            Index j = m / _NJ;
            Index k = n / (_NI*_NJ);
            
            res =  _Interpolation(_PointsI, _NI, uvw[0], i);
            res *= _Interpolation(_PointsJ, _NJ, uvw[1], j);
            res *= _Interpolation(_PointsK, _NK, uvw[2], k);
        }
        
        return res;
    }
    
    
    /*! \brief Return the value of the derivative of the basis function number n
     \param n is a local index.
     \param d the dimension of differentation
     \param uvw are the coordinates of the evaluation point.
     \return The value of the basis function.
     */
    
    Real getValueDerivativeBasisFunc(Index n, Index d, const RealVector & uvw) const
    {
        Real res;
        
        if (_Dim == 1)
        {
            Index i = n % _NI;
            res = _DerivativeInterpolation(_PointsI, _NI, uvw[0], i);
        }
        else if (_Dim == 2)
        {
            Index i = n % _NI;
            Index j = n / _NJ;
            
            if (d==1)
            {
                res =  _DerivativeInterpolation(_PointsI, _NI, uvw[0], i);
                res *= _Interpolation(_PointsJ, _NJ, uvw[1], j);
            }
            else
            {
                res = _Interpolation(_PointsI, _NI, uvw[0], i);
                res *= _DerivativeInterpolation(_PointsJ, _NJ, uvw[1], j);
            }
        }
        else
        {
            Index m = n % (_NI*_NJ);
            Index i = m % _NI;
            Index j = m / _NJ;
            Index k = n / (_NI*_NJ);
            
            
            if (d==1)
            {
                res =  _DerivativeInterpolation(_PointsI, _NI, uvw[0], i);
                res *= _Interpolation(_PointsJ, _NJ, uvw[1], j);
                res *= _Interpolation(_PointsK, _NK, uvw[2], k);
            }
            else if (d==2)
            {
                res = _Interpolation(_PointsI, _NI, uvw[0], i);
                res *= _DerivativeInterpolation(_PointsJ, _NJ, uvw[1], j);
                res *= _Interpolation(_PointsK, _NK, uvw[2], k);
            }
            else
            {
                res = _Interpolation(_PointsI, _NI, uvw[0], i);
                res *= _Interpolation(_PointsJ, _NJ, uvw[1], j);
                res *= _DerivativeInterpolation(_PointsK, _NK, uvw[2], k);
            }
        }
        
        return res;
        
    }
    // ---------------------------------------------------------------------------------//
    
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Extracting from a local solution vector at each quadrature point
     the projection of the solution gradient at each quadrature point.
     \param LocalSolution is the solution at each quadrature point of the element.
     \param LocalSolutionGradient is the solution gradient at each quadrature point of the element.
     This last input argument will be filled during this procedure.
     */
    void GradientInterpolation(const std::vector<Real>& LocalSolution, std::vector<Real>& LocalSolutionGradient) const
    {
        if (_Dim == 1)
        {
            for (Index i = 0; i < _NI; ++i) LocalSolutionGradient[i] = 0.0;
            
            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
            for (Index i2 = 0; i2 < _NI; ++i2)
            LocalSolutionGradient[i2] += LocalSolution[i1] * _DerivativeValuesI[i1][i2];
        }
        if (_Dim == 2)
        {
            if (_NI != _NJ)
            {
                for (Index i = 0; i < 2 * _NI * _NJ; ++i) LocalSolutionGradient[i] = 0.0;
                
                
                // First dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolutionGradient[0 + 2 * i2 + 2 * _NI*j1] += LocalSolution[i1 + _NI*j1] * _DerivativeValuesI[i1][i2];
                
                // Second dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolutionGradient[1 + 2 * i1 + 2 * _NI*j2] += LocalSolution[i1 + _NI*j1] * _DerivativeValuesJ[j1][j2];
                
            }
            else
            {
                const Index M = _NI;
                
                for (Index i = 0; i < 2 * M * M; ++i) LocalSolutionGradient[i] = 0.0;
                
                const Real * p_LocalSolution = LocalSolution.data();
                
                for (Index j1 = 0; j1 < M; ++j1)
                {
                    Real * p_LocalSolutionGradient_I = LocalSolutionGradient.data()  + 2 * M*j1;
                    
                    for (Index i1 = 0; i1 < M; ++i1)
                    {
                        Real * p_LocalSolutionGradient_J = LocalSolutionGradient.data() + 1 + 2 * i1;
                        
                        const Real * p_der_valuesI = _DerivativeValuesI[i1].data();
                        const Real * p_der_valuesJ = _DerivativeValuesJ[j1].data();
                        
                        Real * p_tmp_local_sol_I = p_LocalSolutionGradient_I;
                        Real * p_tmp_local_sol_J = p_LocalSolutionGradient_J;
                        
                        const Real * p_stop = p_der_valuesI+M;
                        Real loc_sol = *p_LocalSolution;
                        
                        while (p_der_valuesI != p_stop)
                        {
                            
                            *p_tmp_local_sol_I += loc_sol * (*p_der_valuesI);
                            *p_tmp_local_sol_J += loc_sol * (*p_der_valuesJ);
                            
                            p_tmp_local_sol_I+=2;
                            p_tmp_local_sol_J+=2*M;
                            
                            ++p_der_valuesI;
                            ++p_der_valuesJ;
                        }
                        
                        ++p_LocalSolution;
                    }
                }
            }
        }
        else if (_Dim == 3)
        {
            for (Index i = 0; i < 3 * _NI*_NJ*_NK; ++i) LocalSolutionGradient[i] = 0.0;
        
            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
            for (Index j1 = 0; j1 < _NJ; ++j1)
            for (Index k1 = 0; k1 < _NK; ++k1)
            for (Index i2 = 0; i2 < _NI; ++i2)
            LocalSolutionGradient[0 + 3 * (i2 + _NI*j1 + _NI*_NJ*k1)] += LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] * _DerivativeValuesI[i1][i2];
   
            // Second dimension.
            for (Index j1 = 0; j1 < _NJ; ++j1)
            for (Index k1 = 0; k1 < _NK; ++k1)
            for (Index i1 = 0; i1 < _NI; ++i1)
            for (Index j2 = 0; j2 < _NJ; ++j2)
            LocalSolutionGradient[1 + 3 * (i1 + _NI*j2 + _NI*_NJ*k1)] += LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] * _DerivativeValuesJ[j1][j2];
 
            // Third dimension.
            for (Index k1 = 0; k1 < _NK; ++k1)
            for (Index i1 = 0; i1 < _NI; ++i1)
            for (Index j1 = 0; j1 < _NJ; ++j1)
            for (Index k2 = 0; k2 < _NK; ++k2)
            LocalSolutionGradient[2 + 3 * (i1 + _NI*j1 + _NI*_NJ*k2)] += LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] * _DerivativeValuesK[k1][k2];
 
        }
    }
    
    
    /*! \brief Compute the inner product between the interpolation of the gradient of
     LocalSolution and the vectorial unknown LocalSolutionGradient
     \param LocalSolutionGradient is homogeneous to the local solution gradient.
     \param LocalSolution will be used to reconstruct the other solution gradient
     appearing in the inner-product and will bear the result of the operation.
     */
    template<bool erase_output = true> void TransposeGradientInterpolation(const std::vector<Real>& LocalSolutionGradient, std::vector<Real>& LocalSolution) const
    {
        if (_Dim == 1)
        {
            if constexpr(erase_output)
                for (Index i = 0; i < _NI; ++i) LocalSolution[i] = 0.0;
            
            // First dimension.
            for (Index i1 = 0; i1 < _NI; ++i1)
            for (Index i2 = 0; i2 < _NI; ++i2)
            LocalSolution[i1] +=  _DerivativeValuesI[i1][i2] *  LocalSolutionGradient[i2];
        }
        if (_Dim == 2)
        {
            //if (_NI != _NJ)
            //{
                if constexpr(erase_output)
                    for (Index i = 0; i < _NI*_NJ; ++i) LocalSolution[i] = 0.0;
                
                //First dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolution[i1 + _NI*j1] += _DerivativeValuesI[i1][i2] * LocalSolutionGradient[0 + 2 * i2 + 2 * _NI*j1];
                
                
                //Second dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolution[i1 + _NI*j1] += _DerivativeValuesJ[j1][j2] * LocalSolutionGradient[1 + 2 * i1 + 2 * _NI*j2];
            //}
            //else
            // {
            //     const Index M = _NI;
                
            //     if constexpr(erase_output)
            //         for (Index i = 0; i < M * M; ++i) LocalSolution[i] = 0.0;
                
            //     Real * p_LocalSolution = LocalSolution.data();
                
            //     for (Index j1 = 0; j1 < M; ++j1)
            //     {
            //         const Real * p_LocalSolutionGradient_I = LocalSolutionGradient.data()  + 2 * M*j1;
                    
            //         for (Index i1 = 0; i1 < M; ++i1)
            //         {
            //             const Real * p_LocalSolutionGradient_J = LocalSolutionGradient.data() + 1 + 2 * i1;
                        
            //             const Real * p_der_valuesI = _DerivativeValuesI[i1].data();
            //             const Real * p_der_valuesJ = _DerivativeValuesJ[j1].data();
                        
            //             const Real * p_tmp_local_sol_I = p_LocalSolutionGradient_I;
            //             const Real * p_tmp_local_sol_J = p_LocalSolutionGradient_J;
                        
            //             const Real * p_stop = p_der_valuesI+M;
                        
            //             while (p_der_valuesI != p_stop)
            //             {
            //                 *p_LocalSolution += (*p_tmp_local_sol_I) * (*p_der_valuesI)
            //                 + (*p_tmp_local_sol_J) * (*p_der_valuesJ);
                            
            //                 p_tmp_local_sol_I+=2;
            //                 p_tmp_local_sol_J+=2*M;
                            
            //                 ++p_der_valuesI;
            //                 ++p_der_valuesJ;
            //             }
                        
            //             ++p_LocalSolution;
            //         }
            //     }
            // }
            
        }
        else if (_Dim == 3)
        {
            if (_NI == _NJ && _NI == _NK)
            {
                const Index M = _NI;
                
                if constexpr(erase_output)
                    for (Index i = 0; i < M*M*M; ++i) LocalSolution[i] = 0.0;
                
                //First dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesI[i1][i2] * LocalSolutionGradient[0 + 3 * (i2 + _NI*j1 + _NI*_NJ*k1)];
                
                //Second dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesJ[j1][j2] * LocalSolutionGradient[1 + 3 * (i1 + _NI*j2 + _NI*_NJ*k1)];
                
                //Third dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index k2 = 0; k2 < _NK; ++k2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesK[k1][k2] * LocalSolutionGradient[2 + 3 * (i1 + _NI*j1 + _NI*_NJ*k2)];
            }
            else
            {
                if constexpr(erase_output)
                    for (Index i = 0; i < _NI*_NJ*_NK; ++i) LocalSolution[i] = 0.0;
                
                //First dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesI[i1][i2] * LocalSolutionGradient[0 + 3 * (i2 +  _NI*j1 +  _NI*_NJ*k1)];
                
                //Second dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesJ[j1][j2] * LocalSolutionGradient[1 + 3 * (i1 + _NI*j2 + _NI*_NJ*k1)];
                
                //Third dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index k2 = 0; k2 < _NK; ++k2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesK[k1][k2] * LocalSolutionGradient[2 + 3 * (i1 + _NI*j1 + _NI*_NJ*k2)];
            }
        }
    }


    void MltAddInterpolationWeightI(const std::vector<Real>& LocalSolutionIn, std::vector<Real>& LocalSolutionOut) const
    {
        if (_Dim == 2)
        {
          

            for (Index j1 = 0; j1 < _NJ; ++j1)
            {
                Real sp = 0;

                for (Index i1 = 0; i1 < _NI; ++i1)
                {
                    sp+=LocalSolutionIn[i1 + _NI*j1]*_InterpolationWeightI[i1];
                }

                for (Index i1 = 0; i1 < _NI; ++i1)
                {
                    LocalSolutionOut[i1 + _NI*j1] += sp*_InterpolationWeightI[i1];
                }
            }
        }
    }


        void MltAddInterpolationWeightJ(const std::vector<Real>& LocalSolutionIn, std::vector<Real>& LocalSolutionOut) const
    {
        if (_Dim == 2)
        {
            

            for (Index i1 = 0; i1 < _NI; ++i1)
            {
                Real sp = 0;

                for (Index j1 = 0; j1 < _NJ; ++j1)
                {
                    sp+=LocalSolutionIn[i1 + _NI*j1]*_InterpolationWeightJ[j1];
                }

                for (Index j1 = 0; j1 < _NJ; ++j1)
                {
                    LocalSolutionOut[i1 + _NI*j1] += sp*_InterpolationWeightJ[j1];
                }
            }
        }
    }

    
    void GradientInterpolationMatrix(std::vector< std::vector<Real> >& BasisFunctionGradient) const
    {
        BasisFunctionGradient.resize(getNumPoints());
        
        for (Index iBasisFunc = 0; iBasisFunc < getNumPoints(); iBasisFunc++)
        {
            BasisFunctionGradient[iBasisFunc].resize(_Dim*getNumPoints(), 0);
            std::vector<Real> LocalSolution(getNumPoints(), 0);
            LocalSolution[iBasisFunc] = 1.0;
            
            GradientInterpolation(LocalSolution, BasisFunctionGradient[iBasisFunc]);
        }
    }
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Extracting from a local solution vector at each quadrature point
     the projection of the solution derivative at each quadrature point.
     \param DerivativeDim dimension along which the derivative is taken
     \param LocalSolution is the solution at each quadrature point of the element.
     \param LocalSolutionGradient is the solution derivative at each quadrature point of the element.
     This last input argument will be filled during this procedure.
     */
    void DerivativeInterpolation(Index DerivativeDim, const std::vector<Real>& LocalSolution, std::vector<Real>& LocalSolutionGradient)
    {
        if (_Dim == 1)
        {
            GradientInterpolation(LocalSolution, LocalSolutionGradient);
        }
        if (_Dim == 2)
        {
            for (Index i = 0; i < _NI*_NJ; ++i) LocalSolutionGradient[i] = 0.0;
            
            if (DerivativeDim == 0)
            {
                // First dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolutionGradient[i2 + _NI*j1] += LocalSolution[i1 + _NI*j1] * _DerivativeValuesI[i1][i2];
            }
            
            else if (DerivativeDim == 1)
            {
                // Second dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolutionGradient[i1 + _NI*j2] += LocalSolution[i1 + _NI*j1] * _DerivativeValuesJ[j1][j2];
                
            }
            
        }
        else if (_Dim == 3)
        {
            for (Index i = 0; i < _NI*_NJ*_NK; ++i) LocalSolutionGradient[i] = 0.0;
            
            if (DerivativeDim == 0)
            {
                
                // First dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolutionGradient[i2 + _NI*j1 + _NI*_NJ*k1] += LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] * _DerivativeValuesI[i1][i2];
                
            }
            
            else if (DerivativeDim == 1)
            {
                
                // Second dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolutionGradient[i1 + _NI*j2 + _NI*_NJ*k1] += LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] * _DerivativeValuesJ[j1][j2];
            }
            
            else if (DerivativeDim == 2)
            {
                
                // Third dimension.
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index k2 = 0; k2 < _NK; ++k2)
                LocalSolutionGradient[i1 + _NI*j1 + _NI*_NJ*k2] += LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] * _DerivativeValuesK[k1][k2];
            }
        }
    }
    
    
    
    
    
    /*! \brief Compute the inner product between the interpolation of the derivative of
     LocalSolution and the vectorial unknown LocalSolutionGradient
     \param DerivativeDim dimension along which the derivative is taken
     \param LocalSolutionGradient is homogeneous to the local solution derivative.
     \param LocalSolution will be used to reconstruct the other solution derivative
     appearing in the inner-product and will bear the result of the operation.
     */
    void TransposeDerivativeInterpolation(Index DerivativeDim, const std::vector<Real>& LocalSolutionGradient, std::vector<Real>& LocalSolution)
    {
        if (_Dim == 1)
        {
            TransposeGradientInterpolation(LocalSolutionGradient, LocalSolution);
        }
        if (_Dim == 2)
        {
            for (Index i = 0; i < _NI*_NJ; ++i) LocalSolution[i] = 0.0;
            
            if (DerivativeDim == 0)
            {
                
                //First dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolution[i1 + _NI*j1] += _DerivativeValuesI[i1][i2] * LocalSolutionGradient[i2 + _NI*j1];
            }
            else if (DerivativeDim == 1)
            {
                //Second dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolution[i1 + _NI*j1] += _DerivativeValuesJ[j1][j2] * LocalSolutionGradient[i1 + _NI*j2];
            }
            
        }
        else if (_Dim == 3)
        {
            for (Index i = 0; i < _NI*_NJ*_NK; ++i) LocalSolution[i] = 0.0;
            
            if (DerivativeDim == 0)
            {
                //First dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index i2 = 0; i2 < _NI; ++i2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesI[i1][i2] * LocalSolutionGradient[i2 + _NI*j1 + _NI*_NJ*k1];
            }
            else if (DerivativeDim == 1)
            {
                //Second dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index j2 = 0; j2 < _NJ; ++j2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesJ[j1][j2] * LocalSolutionGradient[i1 + _NI*j2 + _NI*_NJ*k1];
            }
            else if (DerivativeDim == 2)
            {
                //Third dimension
                for (Index i1 = 0; i1 < _NI; ++i1)
                for (Index j1 = 0; j1 < _NJ; ++j1)
                for (Index k1 = 0; k1 < _NK; ++k1)
                for (Index k2 = 0; k2 < _NK; ++k2)
                LocalSolution[i1 + _NI*j1 + _NI*_NJ*k1] += _DerivativeValuesK[k1][k2] * LocalSolutionGradient[i1 + _NI*j1 + _NI*_NJ*k2];
            }
        }
    }
     
};
// -----------------------------------------------------------------------------------------//



} // OndoMathX





