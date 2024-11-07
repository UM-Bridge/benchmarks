#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
		template <
			Index DIM
		>
		class Homeomorphism {

        public:
            
            // ---------------------------------------------------------------------------------//
            /*! \brief Definition of the dimension as a static variable.*/
            static const Index Dim = DIM;
            // ---------------------------------------------------------------------------------//
            
            static const bool isContinuouslyDifferentiable = false;
            // ---------------------------------------------------------------------------------//
            
            
            
        private:
            
            // ---------------------------------------------------------------------------------//
            enum TypeHomeomorphism
            {
                Homeo_Functions = 0,
                Homeo_Rasterization = 1
            };

            
            // C^0 function for the transformation
            std::function<RealVector(const RealVector&)> _deformation;
            
            // Discontinuous function for the gradient of the transformation
            // The second argument is a direction used to compute the value of the function using the limit  in that direction 
            std::function<std::array<std::array<Real,Dim>,Dim>(const RealVector&, const RealVector&)> _deformation_gradient;
  
            //The transformation can also be defined using a rasterization
            Rasterization _rasterization;
          
            //Type of description of the homeomorphism
            TypeHomeomorphism _type;
            
		public:
            
      

			// ---------------------------------------------------------------------------------//
			/*! \brief Constructor of a homeomorphic transform
			*/
            Homeomorphism(std::function<RealVector(const RealVector&)> deformation,
                          std::function<std::array<std::array<Real,Dim>,Dim>(const RealVector&, const RealVector&)> deformation_gradient)
				: _deformation(deformation), _deformation_gradient(deformation_gradient)
			{
                _type = Homeo_Functions;
			}
            
            Homeomorphism(Rasterization & aRasterization)
                : _rasterization(aRasterization)
            {
                assert(Dim == 3);
                _type = Homeo_Rasterization;
            }
            
			// ---------------------------------------------------------------------------------//


			// ---------------------------------------------------------------------------------//
			/*! \brief Evaluating the position of a (quadrature) point after being transformed.
				\param uvw is the position of a (quadrature) point in the template finite element space.
				\param xyz is the position uvw after transformation.
			*/
            void Eval(const RealVector& uvw, RealVector& xyz) const
			{
                switch (_type)
                {
                    case Homeo_Functions :  xyz = _deformation(uvw); break;
                    case Homeo_Rasterization :
                    {
                        Real Lx = (_rasterization.Nx-1) * _rasterization.dx;
                        Real Ly = (_rasterization.Ny-1) * _rasterization.dy;
                        
                        xyz[0]=_rasterization.Ox + uvw[0] * Lx;
                        xyz[1]=_rasterization.Oy + uvw[1] * Ly;
                        xyz[2]=uvw[2]*Interpolate(_rasterization, xyz[0], xyz[1]);
                    }
                    break;
                }
			}

            bool isAffine() const
			{
				return false;
			}

            void getGradDef(
                const RealVector & uvw,
                const RealVector & dir,
                std::array<std::array<Real, Dim>, Dim>& gradDef) const
            {
                
                switch (_type)
                {
                    case Homeo_Functions : gradDef = _deformation_gradient(uvw,dir); break;
                    case Homeo_Rasterization : gradDef =  _deformation_gradient_from_rasterization(uvw,dir);
                    break;
                }
            }
                
            // ---------------------------------------------------------------------------------//
            
        private :
            
            std::array<std::array<Real, Dim>, Dim> _deformation_gradient_from_rasterization(const RealVector& uvw, const RealVector & dir) const
            {
                assert(Dim == 3);
                
                Real Lx = (_rasterization.Nx-1) * _rasterization.dx;
                Real Ly = (_rasterization.Ny-1) * _rasterization.dy;
                
                RealVector xyz;
                
                xyz[0]=_rasterization.Ox + uvw[0] * Lx;
                xyz[1]=_rasterization.Oy + uvw[1] * Ly;
                
                std::array<Real,2> grad_2D = InterpolateGradient(_rasterization,xyz[0],xyz[1],dir[0]/Lx,dir[1]/Ly);

                std::array<std::array<Real, Dim>, Dim> gradDef;

                gradDef[0][0] = Lx ; gradDef[1][0] = 0.0; gradDef[2][0] = 0.0;
                gradDef[0][1] = 0.0; gradDef[1][1] = Ly ; gradDef[2][1] = 0.0;
                
                gradDef[0][2] = Lx*uvw[2]*grad_2D[0];
                gradDef[1][2] = Ly*uvw[2]*grad_2D[1];
                gradDef[2][2] = Interpolate(_rasterization, xyz[0], xyz[1]);
                    
                return gradDef;
            }
            
		};
 
} // OndoMathX
// -----------------------------------------------------------------------------------------//

 
