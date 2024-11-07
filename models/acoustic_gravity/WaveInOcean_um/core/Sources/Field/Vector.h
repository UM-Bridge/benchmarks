#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
    namespace Field {
        
        template <
        Index Dim,
        Field::Regularity _Regularity = Field::Regularity::Constant
        >
        class Vector {
            
        private:
            
            // For piecewise C0 field (the second argument is a label used to defined each reagion)
            std::function<std::array<Real,Dim>(const RealVector&, const Index&)> _Func_PiecewiseC0;
            
           // Constant field
            std::array<Real,Dim> _Vector;
            
            // function -> C0 field
            std::function<std::array<Real,Dim>(const RealVector&)> _Func_C0;
            
            // List of Matrix -> piecewise constant field
            std::map<Index,std::array<Real,Dim>> _Vector_List;

            // function -> Point data
            std::function<std::array<Real,Dim>(const Index&)> _Func_Point_Data;
             
        public:

            static const Field::Regularity Regularity = _Regularity;

			// ---------------------------------------------------------------------------------//
			Vector() = delete;

            Vector(std::function<std::array<Real,Dim>(const Index&)> aFunc ) : _Func_Point_Data(aFunc)
			{
                assert(Regularity==PointData);
			}

            
            Vector(std::function<std::array<Real,Dim>(const RealVector&)> aFunc ) : _Func_C0(aFunc)
			{
                assert(Regularity==C0);
			}
            
            Vector(std::map<Index,std::array<Real,Dim>> aVectorList) : _Vector_List(aVectorList) 
            {
                 assert(Regularity==Field::Regularity::PiecewiseConstant);
            }
            
            Vector(std::array<Real,Dim> aVector) : _Vector(aVector)
            {
                assert(Regularity== Field::Regularity::Constant);
        
            }
            
            Vector(std::function<std::array<Real,Dim>(const RealVector&, const Index&)> aFunc) : _Func_PiecewiseC0 (aFunc)
            {
                 assert(Regularity==Field::Regularity::PiecewiseC0);
                 
            }
			// ---------------------------------------------------------------------------------//
            template<class ArrayIn,class ArrayOut> void Eval(const ArrayIn& U,
                                            ArrayOut& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                std::array<Real,Dim> A;

                if constexpr(Regularity == Field::Regularity::Constant) A = _Vector;
                if constexpr(Regularity == Field::Regularity::PiecewiseConstant) A = _Vector_List.at(Label);
                if constexpr(Regularity == Field::Regularity::C0) A = _Func_C0(xyz);
                if constexpr(Regularity == Field::Regularity::PiecewiseC0) A = _Func_PiecewiseC0(xyz,Label);
                if constexpr(Regularity == Field::Regularity::PointData) A = _Func_Point_Data(Label);

                ArrayAlgebra::VecMlt(A,U,V);
            }

            template<class ArrayIn,class ArrayOut> void AdjointEval(const ArrayIn& U,
                                                    ArrayOut& V, 
                                                    const RealVector & xyz, 
                                                    const Index & Label) const
            {
                std::array<Real,Dim> A;

                if constexpr(Regularity == Field::Regularity::Constant) A = _Vector;
                if constexpr(Regularity == Field::Regularity::PiecewiseConstant) A = _Vector_List(Label);
                if constexpr(Regularity == Field::Regularity::C0) A = _Func_C0(xyz);
                if constexpr(Regularity == Field::Regularity::PiecewiseC0) A = _Func_PiecewiseC0(xyz,Label);
                if constexpr(Regularity == Field::Regularity::PointData) A = _Func_Point_Data(Label);

                ArrayAlgebra::TransposeVecMlt(A,U,V);
            }
        };

    } // Field

} // OndoMathX
 
