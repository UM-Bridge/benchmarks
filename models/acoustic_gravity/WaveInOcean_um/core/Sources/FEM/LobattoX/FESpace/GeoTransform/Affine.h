#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
		template <
			Index DIM
		>
		class Affine {

		public:
            
            
            // ---------------------------------------------------------------------------------//
            /*! \brief Definition of the dimension as a static variable.*/
            static const Index Dim = DIM;
            // ---------------------------------------------------------------------------------//
 
            static const bool isContinuouslyDifferentiable = true;

			// ---------------------------------------------------------------------------------//
			/*! \brief Constructor of a affine transform from a template base matrix and translation vector.
			*/
			Affine(std::array<std::array<Real, Dim>, Dim> Basis, std::array<Real, Dim> Translation)
				: _Basis(Basis), _Translation(Translation)
			{
				// Storing cofactor matrix.
                ArrayAlgebra::CoMat<Dim>(_Basis, _CoMatrixBasis);
			}
			// ---------------------------------------------------------------------------------//





			// ---------------------------------------------------------------------------------//
			/*! \brief Evaluating the position of a (quadrature) point after being transformed.
				\param uvw is the position of a (quadrature) point in the template finite element space.
				\param xyz is the position uvw after transformation.
			*/
			void Eval(const RealVector& uvw, RealVector& xyz) const
			{
				// Application of affine transform.
				ArrayAlgebra::MltAdd<Dim>(_CoMatrixBasis, uvw, _Translation, xyz);
			}

			bool isAffine() const
			{
				return true;
			}
            
            void getGradDef(
                            const RealVector & uvw,
                std::array<std::array<Real, Dim>, Dim>& gradDef) const
            {
                gradDef = _Basis;
            }
            
            
			// ---------------------------------------------------------------------------------//

		private:

			/*! \brief Basis matrix associated to the affine transformation.*/
			std::array<std::array<Real, Dim>, Dim> _Basis;

			/*! \brief Cofactor basis matrix associated to the affine transformation.*/
			std::array<std::array<Real, Dim>, Dim> _CoMatrixBasis;

			/*! \brief Translation vector associated to the affine transformation.*/
			std::array<Real, Dim> _Translation;
		};
 
} // OndoMathX
// -----------------------------------------------------------------------------------------//

 
