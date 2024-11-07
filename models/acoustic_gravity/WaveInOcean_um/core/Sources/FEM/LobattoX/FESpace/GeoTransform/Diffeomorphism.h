#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
    
 
		template <
			Index DIM
		>
		class Diffeomorphism {

        public:
            
            
            // ---------------------------------------------------------------------------------//
            /*! \brief Definition of the dimension as a static variable.*/
            static const Index Dim = DIM;
            // ---------------------------------------------------------------------------------//
 
            static const bool isContinuouslyDifferentiable = true;
	 
            
        private:
            
            // C^1 function for the transformation
            std::function<RealVector(const RealVector&)> _deformation;
            
            // C^0 function for the gradient of the transformation
            std::function<std::array<std::array<Real,Dim>,Dim>(const RealVector&)> _deformation_gradient;
             
		public:
            
            
     

			// ---------------------------------------------------------------------------------//
			/*! \brief Constructor of a diffeomorphic transform
			*/
            Diffeomorphism(std::function<RealVector(const RealVector & xyz)> deformation,
                           std::function<std::array<std::array<Real,Dim>,Dim>(const RealVector&)> deformation_gradient)
				: _deformation(deformation), _deformation_gradient(deformation_gradient)
			{

			}
			// ---------------------------------------------------------------------------------//


			// ---------------------------------------------------------------------------------//
			/*! \brief Evaluating the position of a point after being transformed.
				\param uvw is the position of a point in the reference finite element space.
				\param xyz is the position uvw after transformation.
			*/
			void Eval(const RealVector & uvw, RealVector& xyz) const
			{
				// Application of the transformation
                xyz = _deformation(uvw);
			}

            bool isAffine() const
			{
				return false;
			}
            
            void getGradDef(
                const RealVector & uvw,
                std::array<std::array<Real, Dim>, Dim>& gradDef) const
            {
                gradDef = _deformation_gradient(uvw);
            }
            
            
			// ---------------------------------------------------------------------------------//

		};
 
} // OndoMathX
// -----------------------------------------------------------------------------------------//

 
