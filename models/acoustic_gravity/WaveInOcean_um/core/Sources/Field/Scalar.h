#pragma once

 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

	namespace Field {

        template
        <
        Field::Regularity _Regularity = Field::Regularity::Constant
        >
		class Scalar{
            
        private:

            // For piecewise C0 field (the second argument is a label used to defined each reagion)
            std::function<Real(const RealVector&, const Index&)> _Func_PiecewiseC0;
            
           // Constant field
            Real _Scalar;
            
            // function -> C0 field
            std::function<Real(const RealVector&)> _Func_C0;
            
            // List of scalar -> piecewise constant field
            std::map<Index,Real> _Scalar_List;

            // function -> Point data
            std::function<Real(const Index&)> _Func_Point_Data;
 
		public:

            static const Field::Regularity Regularity = _Regularity;

			// ---------------------------------------------------------------------------------//
			Scalar() = delete;

            Scalar(std::function<Real(const Index&)> aFunc ) : _Func_Point_Data(aFunc)
			{
                assert(Regularity==PointData);
			}

            Scalar(std::function<Real(const RealVector&)> aFunc ) : _Func_C0(aFunc)
			{
                assert(Regularity==C0);
			}
            
            
            Scalar(std::map<Index,Real> aScalarList) : _Scalar_List(aScalarList) 
            {
                 assert(Regularity==Field::Regularity::PiecewiseConstant);
            }
            
            Scalar(Real aScalar) : _Scalar(aScalar)
            {
                assert(Regularity== Field::Regularity::Constant);
        
            }
            
            Scalar(std::function<Real(const RealVector&, const Index&)> aFunc) : _Func_PiecewiseC0 (aFunc)
            {
                 assert(Regularity==Field::Regularity::PiecewiseC0);
                 
            }
			// ---------------------------------------------------------------------------------//
            template<class Array> void Eval(const Array& U,
                                            Array& V, 
                                            const RealVector & xyz, 
                                            Index Label) const
            {
                Real alpha;

                if constexpr(Regularity == Field::Regularity::Constant) alpha = _Scalar;
                if constexpr(Regularity == Field::Regularity::PiecewiseConstant) alpha = _Scalar_List.at(Label);
                if constexpr(Regularity == Field::Regularity::C0) alpha = _Func_C0(xyz);
                if constexpr(Regularity == Field::Regularity::PiecewiseC0) alpha = _Func_PiecewiseC0(xyz,Label);
                if constexpr(Regularity == Field::Regularity::PointData) alpha = _Func_Point_Data(Label);


                ArrayAlgebra::Scale(alpha,U,V);
            }

            template<class Array> void AdjointEval(const Array& U,
                                                    Array& V, 
                                                    const RealVector & xyz, 
                                                    const Index & Label) const
            {
                Eval(U,V,xyz,Label);
            }

            // ---------------------------------------------------------------------------------//
            
		};
    
	} //Field

} // OndoMathX

// -----------------------------------------------------------------------------------------//
 
