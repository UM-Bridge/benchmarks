#pragma once



// -----------------------------------------------------------------------------------------//
// Class definition
namespace Ondolab
{

namespace GGL
{

template<class FESpace, class ParameterField>
void MltAddStiffnessBndy(FESpace & aFESpace,
                         Index BoundaryLabel,
                         ParameterField & aParameterField,
                         Real Alpha,
                         LAL::Vector& U,
                         Real Beta,
                         LAL::Vector& V)
{
    Real * dataU = LAL::GetData(U);
    Real * dataV = LAL::GetData(V);
    
    Real Copy_Beta = Beta;
    
    // Creating temporary variable to store the result of the unassembled "matrix time vector" operation.
    if (Copy_Beta == 0.0)
    {
        // Setting vector to zero.
        LAL::Scale(0.0, V);
        
        Copy_Beta = 1.;
    }
    
    Real Alpha_Beta = Alpha / Copy_Beta;
    
    
    GaussGaussLobattoElement & GLE = aFESpace.getFEBndy(BoundaryLabel);
    
    Point xyz;
    
    
    // Temporary variables for storing local gradient.
    std::array< std::vector<Real>, FESpace::UnknownDim> TmpGrad;
    for (Index iDim = 0; iDim < FESpace::UnknownDim; iDim++)
    TmpGrad[iDim].resize(GLE.GetNPoints() * (FESpace::GeometricDim-1));
    
    // Temporary variables for storing local solution.
    std::array< std::vector<Real>, FESpace::UnknownDim> TmpSol;
    for (Index iDim = 0; iDim < FESpace::UnknownDim; iDim++)
    TmpSol[iDim].resize(GLE.GetNPoints());
    
    // Temporary variables for storing global index.
    std::vector<Index> IdxGlob(GLE.GetNPoints());
    
    
    for (Index el = 0; el < aFESpace.getNElemBndy(BoundaryLabel); ++el)
    {
        // Storing the position of the barycenter of the element
        Point xyzCenter;
        aFESpace.GetBarycenterBndy(BoundaryLabel, el, xyzCenter);
        
        //Pre-compute the global indices.
        for (Index iLoc = 0; iLoc < GLE.GetNPoints() ; iLoc++)
        IdxGlob[iLoc] = aFESpace.Loc2GlobBndy(BoundaryLabel, el, iLoc)*FESpace::UnknownDim;
        
        // Extract the local solution vector.
        for (Index iLoc = 0; iLoc < GLE.GetNPoints() ; iLoc++)
        {
            Index Idx = IdxGlob[iLoc];
            
            for (Index iDim = 0; iDim < FESpace::UnknownDim; iDim++)
            {
                // Extracting solution.
                Real data = dataU[Idx];
                
                TmpSol[iDim][iLoc] = data;
                
                // Increasing global index to obtain other dimension.
                Idx++;
            }
        }
        
        
        // Compute the interpolated gradient on the quadrature points.
        for (Index iDim = 0; iDim < FESpace::UnknownDim; iDim++)
        GLE.GradientInterpolation(TmpSol[iDim], TmpGrad[iDim]);
        
        
        // Loop on quadrature points.
        for (Index I=0;I<GLE.GetNPoints();++I)
        {
            
            // Recovering the gradient of the unkowns at a quadrature point.
            std::array< std::array<Real, FESpace::GeometricDim-1>, FESpace::UnknownDim> SolutionGradientQuad;
            
            for (Index j = 0; j < FESpace::UnknownDim; ++j)
            for (Index k = 0; k < FESpace::GeometricDim-1; ++k)
            SolutionGradientQuad[j][k] = TmpGrad[j][I * (FESpace::GeometricDim-1) + k];
            
            
            // Multiplying by the comatrix of the transform of the template space.
            if (FESpace::GeometricDim != 2) assert(false);
            
            // Applying constitutive Parameter of the material.
            switch (aParameterField.GetProperty())
            {
                case Field::Constant:
                    aParameterField.Eval(SolutionGradientQuad);
                    break;
                case Field::PiecewiseConstant:
                    aParameterField.Eval(xyzCenter,SolutionGradientQuad);
                    break;
                case Field::Smooth:
                    aFESpace.GetDoFCoordinateBndy(BoundaryLabel, el, I, xyz);
                    aParameterField.Eval(xyz,SolutionGradientQuad);
                    break;
                case Field::None:
                    aFESpace.GetDoFCoordinateBndy(BoundaryLabel, el, I, xyz);
                    aParameterField.Eval(xyz,xyzCenter,SolutionGradientQuad);
                    break;
            }
            
            // Multiplying by the transposed comatrix of the transform of the template space.
            if (FESpace::GeometricDim != 2) assert(false);
            
            // Computing scaling coefficient (quadrature weight divided by the jacobian).
            Real wi_J = GLE.GetQuadratureWeight(I) / aFESpace.GetJacobianBndy(BoundaryLabel, el, I);
            
            // Scaling result and assigns it.
            for (Index j = 0; j < FESpace::UnknownDim; ++j)
            for (Index k = 0; k < FESpace::GeometricDim-1; ++k)
            TmpGrad[j][I * (FESpace::GeometricDim-1) + k] = SolutionGradientQuad[j][k] * wi_J;
        }
        
        
        // Compute the interpolated gradient on the quadrature points.
        for (Index iDim = 0; iDim < FESpace::UnknownDim; iDim++)
        GLE.TransposeGradientInterpolation(TmpGrad[iDim], TmpSol[iDim]);
        
        
        // Summing into the global solution vector.
        for (Index iLoc = 0; iLoc < GLE.GetNPoints(); iLoc++)
        {
            // Extracting global index.
            Index Idx = IdxGlob[iLoc];
            
            for (Index iDim = 0; iDim < FESpace::UnknownDim; iDim++)
            {
                // Adding local contribution.
                dataV[Idx] += Alpha_Beta * TmpSol[iDim][iLoc];
                
                // Updating global index.
                Idx++;
            }
        }
        
        
        
        
        
        
    }
    
    if (Copy_Beta != 1.) LAL::Scale(Copy_Beta, V);
}



} //GGL

} // Ondolab
// -----------------------------------------------------------------------------------------//




