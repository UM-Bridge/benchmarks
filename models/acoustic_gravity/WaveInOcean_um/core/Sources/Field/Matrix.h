#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
    namespace Field {
        
        template <
        Index DimI,
        Index DimJ,
        Field::Regularity _Regularity = Field::Regularity::Constant
        >
        class Matrix {
            
        private:
            
            // For piecewise C0 field (the second argument is a label used to defined each reagion)
            std::function<std::array<std::array<Real,DimJ>,DimI>(const RealVector&, const Index&)> _Func_PiecewiseC0;
            
           // Constant field
            std::array<std::array<Real,DimJ>,DimI> _Matrix;
            
            // function -> C0 field
            std::function<std::array<std::array<Real,DimJ>,DimI>(const RealVector&)> _Func_C0;
            
            // List of Matrix -> piecewise constant field
            std::map<Index,std::array<std::array<Real,DimJ>,DimI>> _Matrix_List;

            // function -> Point data
            std::function<std::array<std::array<Real,DimJ>,DimI>(const Index&)> _Func_Point_Data;
            
            
        public:

            static const Field::Regularity Regularity = _Regularity;

			// ---------------------------------------------------------------------------------//
			Matrix() = delete;

            Matrix(std::function<std::array<std::array<Real,DimJ>,DimI>(const Index&)> aFunc ) : _Func_Point_Data(aFunc)
			{
                assert(Regularity==PointData);
			}

            
            Matrix(std::function<std::array<std::array<Real,DimJ>,DimI>(const RealVector&)> aFunc ) : _Func_C0(aFunc)
			{
                assert(Regularity==C0);
			}
            
            
            Matrix(std::map<Index,std::array<std::array<Real,DimJ>,DimI>> aMatrixList) : _Matrix_List(aMatrixList) 
            {
                 assert(Regularity==Field::Regularity::PiecewiseConstant);
            }
            
            Matrix(std::array<std::array<Real,DimJ>,DimI> aMatrix) : _Matrix(aMatrix)
            {
                assert(Regularity== Field::Regularity::Constant);
        
            }
            
            Matrix(std::function<std::array<std::array<Real,DimJ>,DimI>(const RealVector&, const Index&)> aFunc) : _Func_PiecewiseC0 (aFunc)
            {
                 assert(Regularity==Field::Regularity::PiecewiseC0);
                 
            }
			// ---------------------------------------------------------------------------------//
            template<class Array> void Eval(const Array& U,
                                            Array& V, 
                                            const RealVector & xyz, 
                                            const Index & Label) const
            {
                std::array<std::array<Real,DimJ>,DimI> A;

                if constexpr(Regularity == Field::Regularity::Constant) A = _Matrix;
                if constexpr(Regularity == Field::Regularity::PiecewiseConstant) A = _Matrix_List.at(Label);
                if constexpr(Regularity == Field::Regularity::C0) A = _Func_C0(xyz);
                if constexpr(Regularity == Field::Regularity::PiecewiseC0) A = _Func_PiecewiseC0(xyz,Label);
                if constexpr(Regularity == Field::Regularity::PointData) A = _Func_Point_Data(Label);

                ArrayAlgebra::MatMlt(A,U,V);
            }

            template<class Array> void AdjointEval(const Array& U,
                                                    Array& V, 
                                                    const RealVector & xyz, 
                                                    const Index & Label) const
            {
                std::array<std::array<Real,DimJ>,DimI> A;

                if constexpr(Regularity == Field::Regularity::Constant) A = _Matrix;
                if constexpr(Regularity == Field::Regularity::PiecewiseConstant) A = _Matrix_List(Label);
                if constexpr(Regularity == Field::Regularity::C0) A = _Func_C0(xyz);
                if constexpr(Regularity == Field::Regularity::PiecewiseC0) A = _Func_PiecewiseC0(xyz,Label);
                if constexpr(Regularity == Field::Regularity::PointData) A = _Func_Point_Data(Label);

                ArrayAlgebra::TransposeMatMlt(A,U,V);
            }
            // ---------------------------------------------------------------------------------//
        };

      

    } // Field

} // OndoMathX
 
