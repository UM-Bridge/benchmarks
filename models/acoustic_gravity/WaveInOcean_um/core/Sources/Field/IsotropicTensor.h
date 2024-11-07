#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
    namespace Field {
        
        template <
        Field::Regularity _Regularity = Field::Regularity::Constant
        >
        class IsotropicTensor {
            
        private:
            
            // For piecewise C0 field (the second argument is a label used to defined each reagion)
            std::function<void(const RealVector&, const Index&, Real& Lambda, Real& mu)> _Func_PiecewiseC0;
            
            // Constant field
            Real _Lambda;
            Real _Mu;
            
            // function -> C0 field
            std::function<void(const RealVector&, Real& Lambda, Real& mu)> _Func_C0;
            
            // List of scalars -> piecewise constant field
            std::map<Index,std::pair<Real,Real>> _List;
  
        public:
            
            static const Field::Regularity Regularity = _Regularity;

			// ---------------------------------------------------------------------------------//
			IsotropicTensor() = delete;

            
            IsotropicTensor(std::function<void(const RealVector&, Real&, Real&)> aFunc) : _Func_C0(aFunc)
			{
                assert(Regularity==Field::Regularity::C0);
			}
            
            
            IsotropicTensor(std::map<Index,std::pair<Real,Real>> aList) : _List(aList) 
            {
                 assert(Regularity==Field::Regularity::PiecewiseConstant);
            }
            
            IsotropicTensor(Real aLambda, Real aMu) : _Lambda(aLambda), _Mu(aMu)
            {
                assert(Regularity== Field::Regularity::Constant);
            }
            
            IsotropicTensor(std::function<void(const RealVector&, const Index&, Real&, Real&)> aFunc) : _Func_PiecewiseC0(aFunc) 
            {
                 assert(Regularity==Field::Regularity::PiecewiseC0);
                 
            }

            // ---------------------------------------------------------------------------------//
            template<Index Dim> void Eval(const std::array<std::array<Real, Dim>, Dim>& U,
                                            std::array<std::array<Real, Dim>, Dim>& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                Real lambda;
                Real mu;

                if constexpr(Regularity == Field::Regularity::Constant) 
                {
                    lambda = _Lambda;
                    mu = _Mu;
                }

                if constexpr(Regularity == Field::Regularity::PiecewiseConstant) 
                {
                    auto lambdamu = _List.at(Label);
                    lambda = lambdamu.first;
                    mu = lambdamu.second;
                }

                if constexpr(Regularity == Field::Regularity::C0) 
                {
                    _Func_C0(xyz, lambda, mu);
                }

                if constexpr(Regularity == Field::Regularity::PiecewiseC0) 
                {
                    _Func_PiecewiseC0(xyz, Label, lambda, mu);
                }

                if constexpr(Dim==2)
                {
                   Real LambdaTraceStrainTens = lambda*(U[0][0]+U[1][1]);
                    
                    V[0][0] =  LambdaTraceStrainTens + 2.0 * mu * U[0][0];
                    V[1][1] =  LambdaTraceStrainTens + 2.0 * mu * U[1][1];
                    V[0][1] =  mu * (U[0][1]+U[1][0]);
                    V[1][0] =  V[0][1];
                }

                if constexpr(Dim==3)
                {
                    Real LambdaTraceStrainTens = lambda*(U[0][0]+U[1][1]+U[2][2]);
                    
                    V[0][0] =  LambdaTraceStrainTens + 2.0 * mu * U[0][0];
                    V[1][1] =  LambdaTraceStrainTens + 2.0 * mu * U[1][1];
                    V[2][2] =  LambdaTraceStrainTens + 2.0 * mu * U[2][2];
                    V[0][1] =  mu * (U[0][1]+U[1][0]);
                    V[0][2] =  mu * (U[0][2]+U[2][0]);
                    V[1][2] =  mu * (U[1][2]+U[2][1]);
                    V[1][0] =  V[0][1];
                    V[2][0] =  V[0][2];
                    V[2][1] =  V[1][2];
                }
            }
        };

    } // Field

} // OndoMathX
 
