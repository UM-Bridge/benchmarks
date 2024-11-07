#pragma once

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX
{

    namespace Core
    {

        template <
            Index SysDim = 1,
            Index NQuadNL = 1,
            class FESpace,
            class ParameterField_C>
        void NonLinearStiffnessMltAdd(FESpace &aFESpace,
                                      ParameterField_C &aParameterField_C,
                                      Real Alpha,
                                      LAL::Vector &U,
                                      LAL::Vector &W,
                                      Real Beta,
                                      LAL::Vector &V,
                                      Real tn_m,
                                      Real tn,
                                      Real tn_p)
        {
            constexpr Index DimD = FESpace::Dim;

            assert(SysDim > 0);
            assert(DimD > 1);

            Real *dataU = LAL::getData(U);
            Real *dataV = LAL::getData(V);
            Real *dataW = LAL::getData(W);
            Real Copy_Beta = Beta;

            // Creating temporary variable to store the result of the unassembled "matrix time vector" operation.
            if (Copy_Beta == 0.0)
            {
                // Setting vector to zero.
                LAL::Scale(0.0, V);

                Copy_Beta = 1.;
            }

            Real Alpha_Beta = Alpha / Copy_Beta;

            const auto &FE = aFESpace.getFE();

            // Boolean used to know if the label of the element is used for the parameters computation
            constexpr bool req_label_buffer =
                (ParameterField_C::Regularity == Field::Regularity::PiecewiseC0) ||
                (ParameterField_C::Regularity == Field::Regularity::PiecewiseConstant);

            constexpr bool req_xyz_quadrature_buffer =
                (ParameterField_C::Regularity == Field::Regularity::C0) ||
                (ParameterField_C::Regularity == Field::Regularity::PiecewiseC0);

            // Number of quadrature/interpolation points
            Index NLoc = FE.getNumPoints();

            Index NumColors = aFESpace.getNumColors();

            Real w_0 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;
            Real w_1 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
            Real w_2 = 128.0 / 225.0;
            Real w_3 = (322.0 + 13.0 * sqrt(70.0)) / 900.0;
            Real w_4 = (322.0 - 13.0 * sqrt(70.0)) / 900.0;

            Real xi_0 = -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;
            Real xi_1 = -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
            Real xi_2 = 0.0;
            Real xi_3 = sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0;
            Real xi_4 = sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0;

            Real dt = tn_p - tn_m;

            Real t0 = tn_m + 0.5 * dt * (1.0 + xi_0);
            Real t1 = tn_m + 0.5 * dt * (1.0 + xi_1);
            Real t2 = tn_m + 0.5 * dt * (1.0 + xi_2);
            Real t3 = tn_m + 0.5 * dt * (1.0 + xi_3);
            Real t4 = tn_m + 0.5 * dt * (1.0 + xi_4);

            for (Index iColor = 0; iColor < NumColors; iColor++)
            {
#pragma omp parallel
                {
                    // Buffer for storing the local solution.
                    std::array<std::vector<Real>, SysDim> Buffer_Sol_U;
                    for (Index u = 0; u < SysDim; u++)
                        Buffer_Sol_U[u].resize(NLoc);

                    std::array<std::vector<Real>, SysDim> Buffer_Sol_W;
                    for (Index u = 0; u < SysDim; u++)
                        Buffer_Sol_W[u].resize(NLoc);

                    std::array<std::vector<Real>, SysDim> Buffer_Sol_V;
                    for (Index v = 0; v < SysDim; v++)
                        Buffer_Sol_V[v].resize(NLoc);

                    // Buffer used for gradients
                    std::array<std::vector<Real>, SysDim> Buffer_Grad_U;
                    for (Index u = 0; u < SysDim; u++)
                        Buffer_Grad_U[u].resize(NLoc * DimD);

                    std::array<std::vector<Real>, SysDim> Buffer_Grad_W;
                    for (Index u = 0; u < SysDim; u++)
                        Buffer_Grad_W[u].resize(NLoc * DimD);

                    std::array<std::vector<Real>, SysDim> Buffer_Grad_V;
                    for (Index v = 0; v < SysDim; v++)
                        Buffer_Grad_V[v].resize(NLoc * DimD);

                    // Temporary variables for storing global index.
                    std::vector<Index> IdxGlobU(NLoc);

                    auto Buffer_IJ = InitBuffer<SysDim, DimD>();
                    auto Buffer_IJ_ACC = InitBuffer<SysDim, DimD>();
                    auto Buffer_UD = InitBuffer<SysDim, DimD>();
                    auto Buffer_WD = InitBuffer<SysDim, DimD>();
                    auto Buffer_VD = InitBuffer<SysDim, DimD>();
                    auto Buffer_UWD = InitBuffer<SysDim, DimD>();

                    std::array<std::array<Real, DimD>, DimD> Buffer_GradDef;
                    std::array<std::array<Real, DimD>, DimD> Buffer_ComatGradDef;
                    Real Buffer_Jacobian;
                    Real Buffer_Weight_Jacobian;

                    // Buffer for position informations
                    Index buffer_label;
                    RealVector buffer_xyz_quadrature;

                    // Boolean used to check if the current element is an affine deformation of the reference element
                    bool is_not_affine = true;

                    // Number of elements with the corresponding color
                    Index NElemColor;

                    NElemColor = aFESpace.getNumElements(iColor);

                    // Loop on elements.
#pragma omp for
                    for (Index iElt = 0; iElt < NElemColor; iElt++)
                    {
                        // Extracting global numbering of the element.
                        Index iEltGlob;

                        iEltGlob = aFESpace.EltColorNum2EltGlobalNum(iColor, iElt);

                        // Pre-compute the global indices.
                        aFESpace.Loc2Glob(iEltGlob, IdxGlobU);

                        // Compute the label of the element if needed
                        if (req_label_buffer)
                            buffer_label = aFESpace.getLabel(iEltGlob);

                        // Extract the local solution vector.
                        for (Index iLoc = 0; iLoc < NLoc; iLoc++)
                        {
                            Index Idx = SysDim * IdxGlobU[iLoc];

                            for (Index u = 0; u < SysDim; u++)
                            {
                                // Extracting solution.
                                Buffer_Sol_U[u][iLoc] = dataU[Idx];
                                if constexpr (NQuadNL != 1)
                                    Buffer_Sol_W[u][iLoc] = dataW[Idx];

                                // Increasing global index to obtain other dimension.
                                Idx++;
                            }
                        }

                        for (Index v = 0; v < SysDim; v++)
                            for (Index iLoc = 0; iLoc < NLoc; ++iLoc)
                                Buffer_Sol_V[v][iLoc] = 0.0;

                        for (Index u = 0; u < SysDim; u++)
                        {
                            FE.GradientInterpolation(Buffer_Sol_U[u], Buffer_Grad_U[u]);
                            if constexpr (NQuadNL != 1)
                                FE.GradientInterpolation(Buffer_Sol_W[u], Buffer_Grad_W[u]);
                        }

                        if (aFESpace.isAffine(iEltGlob))
                        {
                            is_not_affine = false;
                            aFESpace.getGradDef(iEltGlob, 0, Buffer_GradDef);
                            Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                            ArrayAlgebra::CoMat(Buffer_GradDef, Buffer_ComatGradDef);
                        }

                        // Loop on quadrature points.
                        for (Index iLoc = 0; iLoc < NLoc; ++iLoc)
                        {
                            Real Weight = FE.getQuadratureWeight(iLoc);

                            if (is_not_affine)
                            {
                                aFESpace.getGradDef(iEltGlob, iLoc, Buffer_GradDef);
                                Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                ArrayAlgebra::CoMat(Buffer_GradDef, Buffer_ComatGradDef);
                            }

                            if constexpr (req_xyz_quadrature_buffer)
                                aFESpace.getDoFCoordinate(iEltGlob, iLoc, buffer_xyz_quadrature);

                            if constexpr (SysDim == 1)
                            {
                                // Recovering the gradient of the unkowns at a quadrature point
                                for (Index d = 0; d < DimD; d++)
                                    Buffer_UD[d] = Buffer_Grad_U[0][iLoc * DimD + d];

                                if constexpr (NQuadNL != 1)
                                    for (Index d = 0; d < DimD; d++)
                                        Buffer_WD[d] = Buffer_Grad_W[0][iLoc * DimD + d];
                            }
                            else
                            {
                                // Recovering the gradient of the unkowns at a quadrature point
                                for (Index u = 0; u < SysDim; ++u)
                                    for (Index d = 0; d < DimD; d++)
                                        Buffer_UD[u][d] = Buffer_Grad_U[u][iLoc * DimD + d];

                                if constexpr (NQuadNL != 1)
                                    for (Index u = 0; u < SysDim; ++u)
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_WD[u][d] = Buffer_Grad_W[u][iLoc * DimD + d];
                            }

                            // Multiply each gradient by the comatrix of the deformation from the reference element to the current element
                            if constexpr (SysDim == 1)
                            {
                                ArrayAlgebra::MatMlt(Buffer_ComatGradDef, Buffer_UD);
                                ArrayAlgebra::Scale(1.0 / Buffer_Jacobian, Buffer_UD);

                                if constexpr (NQuadNL != 1)
                                {
                                    ArrayAlgebra::MatMlt(Buffer_ComatGradDef, Buffer_WD);
                                    ArrayAlgebra::Scale(1.0 / Buffer_Jacobian, Buffer_WD);
                                }
                            }
                            else
                            {
                                for (Index u = 0; u < SysDim; ++u)
                                {
                                    ArrayAlgebra::MatMlt(Buffer_ComatGradDef, Buffer_UD[u]);
                                    ArrayAlgebra::Scale(1.0 / Buffer_Jacobian, Buffer_UD[u]);
                                }

                                if constexpr (NQuadNL != 1)
                                {
                                    for (Index u = 0; u < SysDim; ++u)
                                    {
                                        ArrayAlgebra::MatMlt(Buffer_ComatGradDef, Buffer_WD[u]);
                                        ArrayAlgebra::Scale(1.0 / Buffer_Jacobian, Buffer_WD[u]);
                                    }
                                }
                            }

                            if constexpr (NQuadNL == 1)
                            {
                                aParameterField_C.Eval(Buffer_UD, Buffer_VD, buffer_xyz_quadrature, buffer_label);
                                ArrayAlgebra::Scale(tn_p - tn_m, Buffer_VD);
                            }

                            // if constexpr (NQuadNL == 2)
                            // {
                            //     Real dt = tn_p - tn_m;
                            //     Real xi = 1.0 / sqrt(3);
                            //     // Second order gauss rule
                            //     ArrayAlgebra::Scale(0.0, Buffer_VD);

                            //     Buffer_UWD = Buffer_UD;
                            //     ArrayAlgebra::MatScaleAdd((0.5 - xi * 0.5) * dt, Buffer_WD, 1.0, Buffer_UWD);
                            //     aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                            //     ArrayAlgebra::MatScaleAdd(0.5 * dt, Buffer_IJ, 1.0, Buffer_VD);

                            //     Buffer_UWD = Buffer_UD;
                            //     ArrayAlgebra::MatScaleAdd(+0.5 * dt / sqrt(3), Buffer_WD, 1.0, Buffer_UWD);
                            //     aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                            //     ArrayAlgebra::MatScaleAdd(0.5 * dt, Buffer_IJ, 1.0, Buffer_VD);
                            // }

                            if constexpr (NQuadNL == 5)
                            {

                                //  gauss rule with 5 points
                                ArrayAlgebra::Scale(0.0, Buffer_VD);

                                Buffer_UWD = Buffer_UD;
                                ArrayAlgebra::MatScaleAdd(t0 - tn, Buffer_WD, 1.0, Buffer_UWD);
                                aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                                ArrayAlgebra::MatScaleAdd(0.5 * dt * w_0, Buffer_IJ, 1.0, Buffer_VD);

                                Buffer_UWD = Buffer_UD;
                                ArrayAlgebra::MatScaleAdd(t1 - tn, Buffer_WD, 1.0, Buffer_UWD);
                                aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                                ArrayAlgebra::MatScaleAdd(0.5 * dt * w_1, Buffer_IJ, 1.0, Buffer_VD);

                                Buffer_UWD = Buffer_UD;
                                ArrayAlgebra::MatScaleAdd(t2 - tn, Buffer_WD, 1.0, Buffer_UWD);
                                aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                                ArrayAlgebra::MatScaleAdd(0.5 * dt * w_2, Buffer_IJ, 1.0, Buffer_VD);

                                Buffer_UWD = Buffer_UD;
                                ArrayAlgebra::MatScaleAdd(t3 - tn, Buffer_WD, 1.0, Buffer_UWD);
                                aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                                ArrayAlgebra::MatScaleAdd(0.5 * dt * w_3, Buffer_IJ, 1.0, Buffer_VD);

                                Buffer_UWD = Buffer_UD;
                                ArrayAlgebra::MatScaleAdd(t4 - tn, Buffer_WD, 1.0, Buffer_UWD);
                                aParameterField_C.Eval(Buffer_UWD, Buffer_IJ, buffer_xyz_quadrature, buffer_label);
                                ArrayAlgebra::MatScaleAdd(0.5 * dt * w_4, Buffer_IJ, 1.0, Buffer_VD);
                            }

                            if constexpr (SysDim == 1)
                            {
                                ArrayAlgebra::TransposeMatMlt(Buffer_ComatGradDef, Buffer_VD);
                            }
                            else
                            {
                                for (Index v = 0; v < SysDim; ++v)
                                    ArrayAlgebra::TransposeMatMlt(Buffer_ComatGradDef, Buffer_VD[v]);
                            }

                            if constexpr (SysDim == 1)
                            {
                                for (Index d = 0; d < DimD; d++)
                                    Buffer_Grad_V[0][iLoc * DimD + d] = Weight * Buffer_VD[d];
                            }
                            else
                            {
                                for (Index v = 0; v < SysDim; ++v)
                                    for (Index d = 0; d < DimD; d++)
                                        Buffer_Grad_V[v][iLoc * DimD + d] = Weight * Buffer_VD[v][d];
                            }
                        }

                        for (Index v = 0; v < SysDim; v++)
                            FE.template TransposeGradientInterpolation<false>(Buffer_Grad_V[v], Buffer_Sol_V[v]);

                        // Summing into the global solution vector.
                        for (Index iLoc = 0; iLoc < NLoc; iLoc++)
                        {
                            Index Idx;

                            Idx = SysDim * IdxGlobU[iLoc];

                            for (Index v = 0; v < SysDim; v++)
                            {
                                // Adding local contribution.
                                dataV[Idx] += Alpha_Beta * Buffer_Sol_V[v][iLoc];

                                // Updating global index.
                                Idx++;
                            }
                        }
                    }
                }
            }
        }

    } // Core

} // OndoMathX
