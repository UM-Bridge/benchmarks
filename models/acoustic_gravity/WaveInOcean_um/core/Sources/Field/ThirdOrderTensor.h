#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
    namespace Field {
        
        template <
        Index DIM_I,
        Index DIM_J = DIM_I,
        Index DIM_K = DIM_I
        >
        class ThirdOrderTensor {
            
        private:
            

            //Function that takes as input a matrix and outputs a vector 
            std::function<void(const RealVector & xyz,
                               const std::array<std::array<Real, DIM_I>, DIM_J> & A,
                               std::array<Real, DIM_K> & b)> _func;

            //Its adjoint 
            std::function<void(const RealVector & xyz,
                               const std::array<Real, DIM_K> & b,
                               std::array<std::array<Real, DIM_I>, DIM_J> & A)> _adj_func;
        
        public:
            
            static const Field::Regularity Regularity = Field::Regularity::C0;

            // ---------------------------------------------------------------------------------//
            ThirdOrderTensor() = delete;
            

            // Paraxial equation 45°
            ThirdOrderTensor(std::function<void(const RealVector & xyz,
                                                const std::array<std::array<Real, DIM_I>, DIM_J> & A,
                                                std::array<Real, DIM_K> & b)> aFunc,
                            std::function<void(const RealVector & xyz,
                            const std::array<Real, DIM_K> & b,
                            std::array<std::array<Real, DIM_I>, DIM_J> & A)> aAdjFunc) :  _func(aFunc), _adj_func(aAdjFunc)
            {
                
            }
            

            // ---------------------------------------------------------------------------------//

            template<Index Dim_I, Index Dim_J, Index Dim_K> void Eval(const std::array<std::array<Real, Dim_I>, Dim_J>& U,
                                          std::array<Real, Dim_K>& V, 
                                          const RealVector & xyz, 
                                          const Index & Label) const
            {
                _func(xyz,U,V);
            }
            
            template<Index Dim_I, Index Dim_J, Index Dim_K> void AdjointEval(const std::array<Real, DIM_K>& U,  
                                                       std::array<std::array<Real, Dim_I>, Dim_J>& V,
                                                 const RealVector & xyz, 
                                                 const Index & Label) const
            {
                _adj_func(xyz,U,V);
            }

        };

    } // Field

} // OndoMathX
 
