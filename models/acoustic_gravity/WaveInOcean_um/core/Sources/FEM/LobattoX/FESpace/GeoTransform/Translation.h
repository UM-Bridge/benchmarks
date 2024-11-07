#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 
		template <
			Index DIM
		>
		class Translation {

		public:




			// ---------------------------------------------------------------------------------//
			/*! \brief Definition of the dimension as a static variable.*/
			static const Index Dim = DIM;

			// ---------------------------------------------------------------------------------//

        
            static const bool isContinuouslyDifferentiable = true;
		 
            // ---------------------------------------------------------------------------------//
            /*! \brief Constructor.
                \param TVector is the translation vector.
            */
            Translation(std::array<Real, Dim> TVector) : _TVector(TVector)
            {

            }
            // ---------------------------------------------------------------------------------//



			// ---------------------------------------------------------------------------------//
			/*! \brief Evaluating the position of a (quadrature) point after being transformed.
				\param uvw is the position of a (quadrature) point in the template finite element space.
				\param xyz is the position uvw after transformation.
			*/
			void Eval(const RealVector &uvw, RealVector& xyz) const
			{
				for (Index i = 0; i < Dim; i++) xyz[i] = uvw[i] + _TVector[i];
			}



			/*! \brief Extracting jacobian of the transformation at a given node, equals to 1.
			*/
			bool isAffine() const
			{
				return true;
			}

            void getGradDef(
                            const RealVector & uvw,
                std::array<std::array<Real, Dim>, Dim>& gradDef) const
            {
                for (Index i = 0; i < Dim; i++)
                    for (Index j = 0; j < Dim; j++)
                    {
                        if(i==j) gradDef[i][i] = 1.0;
                        else gradDef[i][j] = 0.0;
                    }
            }
            
 
			// ---------------------------------------------------------------------------------//

		private:

			/*! \brief Translation vector.*/
			std::array<Real, Dim> _TVector;
		};
 
 
} //OndoMathX

// -----------------------------------------------------------------------------------------//
 
