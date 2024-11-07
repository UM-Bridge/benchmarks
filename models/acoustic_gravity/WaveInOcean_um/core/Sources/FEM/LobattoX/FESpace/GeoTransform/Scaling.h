#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 
        template <
        Index DIM
        >
        class Scaling {
            
        public:
            
            
            // ---------------------------------------------------------------------------------//
            /*! \brief Definition of the dimension as a static variable.*/
            static const Index Dim = DIM;
            //----------------------------------------------------------------------------//
        
            static const bool isContinuouslyDifferentiable = true;
            
            // ---------------------------------------------------------------------------------//
            /*! \brief Constructor from a scaling factor.
             */
            Scaling(Real ScaleX = 1.0, Real ScaleY = 1.0, Real ScaleZ = 1.0)
            {
                _Scale[0]=ScaleX;
                _Scale[1]=ScaleY;
                if (DIM == 3) _Scale[2]=ScaleZ;
            }
            // ---------------------------------------------------------------------------------//
            
            

            // ---------------------------------------------------------------------------------//
            /*! \brief Evaluating the position of a (quadrature) point after being transformed.
             \param uvw is the position of a (quadrature) point in the template finite element space.
             \param xyz is the position uvw after transformation.
             */
            void Eval(const RealVector &uvw, RealVector& xyz) const
            {
                switch (Dim)
                {
                    case 1: xyz[0]=_Scale[0] * uvw[0]; xyz[1]=0.0;              xyz[2]=0.0;              break;
                    case 2: xyz[0]=_Scale[0] * uvw[0]; xyz[1]=_Scale[1]*uvw[1]; xyz[2]=0.0;              break;
                    case 3: xyz[0]=_Scale[0] * uvw[0]; xyz[1]=_Scale[1]*uvw[1]; xyz[2]=_Scale[2]*uvw[2]; break;
                }
            }
 
            Real getScaleX() const
            {
                return _Scale[0];     
            }

            Real getScaleY() const
            {
                return _Scale[1];     
            }

            Real getScaleZ() const
            {
                return _Scale[2];     
            }

        
        	bool isAffine() const
			{
				return true;
			}
            
            void getGradDef(
                            const RealVector & uvw,
                std::array<std::array<Real, Dim>, Dim>& gradDef) const
            {
                if constexpr(Dim == 1)
                {
                    gradDef[0][0] = _Scale[0]; 
                }

                if constexpr(Dim == 2)
                {
                    gradDef[0][0] = _Scale[0]; gradDef[0][1] = 0.0;
                    gradDef[1][0] = 0.0;       gradDef[1][1] = _Scale[1];
                    
                    
                }
                
                if constexpr(Dim == 3)
                {
                    gradDef[0][0] = _Scale[0]; gradDef[0][1] = 0.0;       gradDef[0][2] = 0.0;
                    gradDef[1][0] = 0.0;       gradDef[1][1] = _Scale[1]; gradDef[1][2] = 0.0;
                    gradDef[2][0] = 0.0;       gradDef[2][1] = 0.0;       gradDef[2][2] = _Scale[2]; 
                }
            }
            // ---------------------------------------------------------------------------------//
            
        private:
            
            /*! \brief Scaling factors vector associated to the transformation.*/
            std::array<Real,Dim> _Scale;
        };
    
    
    
}

// -----------------------------------------------------------------------------------------//
 
