#pragma once


// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 
        template <
        Index DIM
        >
		class Identity {

		public:

			// ---------------------------------------------------------------------------------//
			/*! \brief Definition of the dimension as a static variable.*/
			static const Index Dim = DIM;
			// ---------------------------------------------------------------------------------//
 
            static const bool isContinuouslyDifferentiable = true;
	 
			// ---------------------------------------------------------------------------------//
			/*! \brief Evaluating the position of a (quadrature) point after being transformed.
				\param uvw is the position of a (quadrature) point in the template finite element space.
				\param xyz is the position uvw after transformation.
			*/
			void Eval(const RealVector & uvw, RealVector& xyz) const
			{
                xyz = uvw;
			}

 
			// ---------------------------------------------------------------------------------//
            
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
            
		};
    
 
} // OndoMathX
// -----------------------------------------------------------------------------------------//
 
