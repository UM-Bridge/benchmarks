#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
    namespace Field {
        
        template <
        Index DIM,
        Field::Regularity _Regularity = Field::Regularity::C0
        >
        class FourthOrderTensor {
            
        private:
            
            //In Elasticity : function that takes as input a gradient of a displacement and outputs a stress (DU is overwritten)
            std::function<void(const RealVector & xyz,
                               std::array<std::array<Real, DIM>, DIM> & DU)> _func_C0;

            std::function<void(const Index & i,
                               std::array<std::array<Real, DIM>, DIM> & DU)> _func_PointData;
            
            
        public:
            
            static const Field::Regularity Regularity = _Regularity;

            // ---------------------------------------------------------------------------------//
            FourthOrderTensor() = delete;
            
       
            
            FourthOrderTensor(std::function<void(const RealVector & xyz,
                                                       std::array<std::array<Real, DIM>, DIM> & DU)> aFunc) :  _func_C0(aFunc)
            {
                assert(Regularity == Field::Regularity::C0);
            }


            FourthOrderTensor(std::function<void(const Index & i,
                                                std::array<std::array<Real, DIM>, DIM> & DU)> aFunc) :  _func_PointData(aFunc)
            {
                assert(Regularity == Field::Regularity::PointData);
            }

            // ---------------------------------------------------------------------------------//
            
            template<Index Dim> void Eval(const std::array<std::array<Real, Dim>, Dim>& U,
                                            std::array<std::array<Real, Dim>, Dim>& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                V = U;
             
                if constexpr(Regularity == Field::Regularity::C0)
                {
                    _func_C0(xyz,V);
                }
                else
                {
                    assert(false);
                }
            }

            template<Index Dim> void Eval(const std::array<std::array<Real, Dim>, Dim>& U,
                                            std::array<std::array<Real, Dim>, Dim>& V, 
                                            const Index & i) const
            {
                V = U;
             
                if constexpr(Regularity == Field::Regularity::PointData)
                {
                    _func_PointData(i,V);
                }
                else
                {
                    assert(false);
                }
            }
        };

    } // Field

} // OndoMathX
 
