#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
    namespace Field {
        
        template <
        Field::Regularity _Regularity = Field::Regularity::Constant
        >
        class NeoHookeanTensor {
            
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

            void _getLambdaMu(const RealVector & xyz, const Index & Label, Real & lambda, Real & mu) const
            {
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
            }




			// ---------------------------------------------------------------------------------//
			NeoHookeanTensor() = delete;

            
            NeoHookeanTensor(std::function<void(const RealVector&, Real&, Real&)> aFunc) : _Func_C0(aFunc)
			{
                assert(Regularity==Field::Regularity::C0);
			}
            
            
            NeoHookeanTensor(std::map<Index,std::pair<Real,Real>> aList) : _List(aList) 
            {
                 assert(Regularity==Field::Regularity::PiecewiseConstant);
            }
            
            NeoHookeanTensor(Real aLambda, Real aMu) : _Lambda(aLambda), _Mu(aMu)
            {
                assert(Regularity== Field::Regularity::Constant);
            }
            
            NeoHookeanTensor(std::function<void(const RealVector&, const Index&, Real&, Real&)> aFunc) : _Func_PiecewiseC0(aFunc) 
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

                _getLambdaMu(xyz,Label,lambda,mu);
              

                if constexpr(Dim==2)
                {
                    assert(false);
                }

                if constexpr(Dim==3)
                {
                  /*
                    Real LambdaTraceStrainTens = lambda*(U[0][0]+U[1][1]+U[2][2]);
                    
                    V[0][0] =  LambdaTraceStrainTens + 2.0 * mu * U[0][0];
                    V[1][1] =  LambdaTraceStrainTens + 2.0 * mu * U[1][1];
                    V[2][2] =  LambdaTraceStrainTens + 2.0 * mu * U[2][2];
                    V[0][1] =  mu * (U[0][1]+U[1][0]);
                    V[0][2] =  mu * (U[0][2]+U[2][0]);
                    V[1][2] =  mu * (U[1][2]+U[2][1]);
                    V[1][0] =  V[0][1];
                    V[2][0] =  V[0][2];
                    V[2][1] =  V[1][2];  */
 
 
                    std::array<std::array<Real,Dim>, Dim> F; 
                    std::array<std::array<Real,Dim>, Dim> invF; 
                    Real J;

                    F = U;

                    F[0][0]+=1.0; 
                    F[1][1]+=1.0; 
                    F[2][2]+=1.0; 

                    ArrayAlgebra::Inv(F,invF);
                    J = ArrayAlgebra::Det(F);
                    ArrayAlgebra::Transpose(invF,V,lambda*J*(J-1)-mu);
                    ArrayAlgebra::MatScaleAdd(mu,F,1.0,V);   
                }
            }

            template<Index Dim> void Eval(const std::array<Real, Dim*Dim>& U,
                                            std::array<std::array<Real, Dim>, Dim>& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                Real lambda;
                Real mu;

                _getLambdaMu(xyz,Label,lambda,mu);
              

                if constexpr(Dim!=3)
                {
                    assert(false);
                }

                if constexpr(Dim==3)
                {
                    
                    std::array<std::array<Real,Dim>, Dim> F; 
                    std::array<std::array<Real,Dim>, Dim> invF; 
                    Real J;

                    F[0][0] = U[0]+1.0;
                    F[0][1] = U[1];
                    F[0][2] = U[2];
                    F[1][0] = U[3];
                    F[1][1] = U[4]+1.0;
                    F[1][2] = U[5];
                    F[2][0] = U[6];
                    F[2][1] = U[7];
                    F[2][2] = U[8]+1.0;

                    ArrayAlgebra::Inv(F,invF);
                    J = ArrayAlgebra::Det(F);
                    ArrayAlgebra::Transpose(invF,V,lambda*J*(J-1)-mu);
                    ArrayAlgebra::MatScaleAdd(mu,F,1.0,V);   
                }
            }

            // ---------------------------------------------------------------------------------//
            template<Index Dim> void Eval(const std::array<std::array<Real, Dim>, Dim>& U,
                                            Real& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                Real lambda;
                Real mu;
                
                _getLambdaMu(xyz,Label,lambda,mu);
              

                if constexpr(Dim==2)
                {
                    assert(false);
                }

                if constexpr(Dim==3)
                {
                    /* 
                    std::array<std::array<Real,Dim>, Dim> F;


                    Real LambdaTraceStrainTens = lambda*(U[0][0]+U[1][1]+U[2][2]);
                    
                    F[0][0] =  LambdaTraceStrainTens + 2.0 * mu * U[0][0];
                    F[1][1] =  LambdaTraceStrainTens + 2.0 * mu * U[1][1];
                    F[2][2] =  LambdaTraceStrainTens + 2.0 * mu * U[2][2];
                    F[0][1] =  mu * (U[0][1]+U[1][0]);
                    F[0][2] =  mu * (U[0][2]+U[2][0]);
                    F[1][2] =  mu * (U[1][2]+U[2][1]);
                    F[1][0] =  F[0][1];
                    F[2][0] =  F[0][2];
                    F[2][1] =  F[1][2];  

                    V = 0.5*ArrayAlgebra::ContractionProduct(F,U);*/
                   
              
                    std::array<std::array<Real,Dim>, Dim> F; 
                  
                    Real J;

                    F = U;

                    F[0][0]+=1.0; 
                    F[1][1]+=1.0; 
                    F[2][2]+=1.0; 

                    J = ArrayAlgebra::Det(F);

                    Real I3 = J*J;
                    Real I1 = ArrayAlgebra::NormFSquared(F);
                   
                    V =  0.5*lambda*(I3 - 2*J +1) + 0.5*mu*(I1-3-log(I3)); 
                }
            }       


            // ---------------------------------------------------------------------------------//
            template<Index Dim> void Eval(const std::array<Real, Dim>& U,
                                            Real& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                Real lambda;
                Real mu;
                
                _getLambdaMu(xyz,Label,lambda,mu);
              

                if constexpr(Dim!=9)
                {
                    assert(false);
                }

                if constexpr(Dim==9)
                {
                    std::array<std::array<Real,3>, 3> F; 
                  
                    Real J;

                        
                    F[0][0] = U[0]+1.0;
                    F[0][1] = U[1];
                    F[0][2] = U[2];
                    F[1][0] = U[3];
                    F[1][1] = U[4]+1.0;
                    F[1][2] = U[5];
                    F[2][0] = U[6];
                    F[2][1] = U[7];
                    F[2][2] = U[8]+1.0;
                  
            

                    J = ArrayAlgebra::Det(F);

                    Real I3 = J*J;
                    Real I1 = ArrayAlgebra::NormFSquared(F);
                   
                    V =  0.5*lambda*(I3 - 2*J +1) + 0.5*mu*(I1-3-log(I3)); 
                }
            }                  
        };

    } // Field

} // OndoMathX
 
