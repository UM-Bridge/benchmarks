#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
namespace Core {


template<Field::Regularity R, class ParameterField, class DataIn, class DataOut> void
    Eval(ParameterField & field,DataIn &in, DataOut &out, RealVector &xyz, Index &label, Index &iglob )
    {
        if constexpr(R == Field::Regularity::PointData)
        {
            field.Eval(in,out,iglob);
        }
        else
        {
            field.Eval(in,out,xyz,label);
        }
    }

template<
    class FESpaceU, 
    class FESpaceV,  
    bool Same_FESpace,
    bool ContinuousU, 
    bool ContinuousV,
    QuadratureOpt QuadOpt,
    Index DimU,
    Index DimV,
    Index DimI,
    Index DimJ,
    class ParameterField_C,
    DiffOp E1,
    DiffOp F1,
    class ParameterField_L1 = OndoMathX::Field::Identity,
    class ParameterField_M1 = OndoMathX::Field::Identity,
    DiffOp E2 = Zero,
    DiffOp F2 = Zero,
    class ParameterField_L2 = OndoMathX::Field::Identity,
    class ParameterField_M2 = OndoMathX::Field::Identity,
    DiffOp E3 = Zero,
    DiffOp F3 = Zero,
    class ParameterField_L3 = OndoMathX::Field::Identity,
    class ParameterField_M3 = OndoMathX::Field::Identity,
    bool Check_zero = false,
    class FEQuad = DefaultFE
    >
void MltAdd(FESpaceU & aFESpaceU,
            FESpaceV & aFESpaceV,
            ParameterField_C & aParameterField_C,
            Real Alpha,
            const LAL::Vector& U,
            Real Beta,
            LAL::Vector& V,
            ParameterField_L1 & aParameterField_L1 = Field::IdentityField,
            ParameterField_M1 & aParameterField_M1 = Field::IdentityField,
            ParameterField_L2 & aParameterField_L2 = Field::IdentityField,
            ParameterField_M2 & aParameterField_M2 = Field::IdentityField,
            ParameterField_L3 & aParameterField_L3 = Field::IdentityField,
            ParameterField_M3 & aParameterField_M3 = Field::IdentityField,
            FEQuad & FE_Q = _DefaultFE
            )
            {
                constexpr Index DimD = FESpaceU::Dim;

                static_assert(DimD == FESpaceV::Dim);

                static_assert(DimU>0);
                static_assert(DimV>0);
                static_assert(DimI>0);
                static_assert(DimJ>0);

                if constexpr(Same_FESpace == true)
                {
                    assert(&aFESpaceU == &aFESpaceV);
                }

                // Recovering pointer to data and compute scaling coefficient//////////////////////////////////////////
                const Real * dataU = LAL::getData(U);
                Real * dataV = LAL::getData(V);
    
                if (Beta == 0.0)
                {
                    // Setting vector to zero.
                    LAL::Scale(0.0, V);
        
                    Beta = 1.;
                }
    
                Real Alpha_Beta = Alpha / Beta;
                ///////////////////////////////////////////////////////////////////////////////////////////////////////

 
    
                //Boolean used to know if the label of the element is used for the parameters computation//////////////
                constexpr bool req_index_buffer =  
                            (   ParameterField_C::Regularity == Field::Regularity::PointData
                            ||  ParameterField_L1::Regularity == Field::Regularity::PointData
                            ||  ParameterField_L2::Regularity == Field::Regularity::PointData
                            ||  ParameterField_L3::Regularity == Field::Regularity::PointData
                            ||  ParameterField_M1::Regularity == Field::Regularity::PointData
                            ||  ParameterField_M2::Regularity == Field::Regularity::PointData
                            ||  ParameterField_M3::Regularity == Field::Regularity::PointData);

                constexpr bool req_label_buffer = 
                             (   ParameterField_C::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_L1::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_L2::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_L3::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_M1::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_M2::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_M3::Regularity == Field::Regularity::PiecewiseC0) || 
                            (    ParameterField_C::Regularity == Field::Regularity::PiecewiseConstant
                            ||  ParameterField_L1::Regularity == Field::Regularity::PiecewiseConstant
                            ||  ParameterField_L2::Regularity == Field::Regularity::PiecewiseConstant
                            ||  ParameterField_L3::Regularity == Field::Regularity::PiecewiseConstant
                            ||  ParameterField_M1::Regularity == Field::Regularity::PiecewiseConstant
                            ||  ParameterField_M2::Regularity == Field::Regularity::PiecewiseConstant
                            ||  ParameterField_M3::Regularity == Field::Regularity::PiecewiseConstant);

                constexpr bool req_xyz_quadrature_buffer = 
                             (   ParameterField_C::Regularity == Field::Regularity::C0
                            ||  ParameterField_L1::Regularity == Field::Regularity::C0
                            ||  ParameterField_L2::Regularity == Field::Regularity::C0
                            ||  ParameterField_L3::Regularity == Field::Regularity::C0
                            ||  ParameterField_M1::Regularity == Field::Regularity::C0
                            ||  ParameterField_M2::Regularity == Field::Regularity::C0
                            ||  ParameterField_M3::Regularity == Field::Regularity::C0) || 
                             (   ParameterField_C::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_L1::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_L2::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_L3::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_M1::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_M2::Regularity == Field::Regularity::PiecewiseC0
                            ||  ParameterField_M3::Regularity == Field::Regularity::PiecewiseC0); 

                constexpr bool compute_gradient_U = 
                    (   E1 == DiffOp::Gradient 
                    ||  E2 == DiffOp::Gradient 
                    ||  E3 == DiffOp::Gradient);

                constexpr bool compute_divergence_U = 
                    (   E1 == DiffOp::Divergence 
                    ||  E2 == DiffOp::Divergence 
                    ||  E3 == DiffOp::Divergence);

                constexpr bool req_gradient_buffer_U = compute_gradient_U || compute_divergence_U;
    
                constexpr bool compute_gradient_V = 
                    (   F1 == DiffOp::Gradient 
                    ||  F2 == DiffOp::Gradient 
                    ||  F3 == DiffOp::Gradient);

                constexpr bool compute_divergence_V = 
                    (   F1 == DiffOp::Divergence 
                    ||  F2 == DiffOp::Divergence 
                    ||  F3 == DiffOp::Divergence);

                constexpr bool req_gradient_buffer_V = compute_gradient_V || compute_divergence_V; 

                constexpr bool trace_operator_U = ( E1 == DiffOp::Trace       || E2 == DiffOp::Trace       || E3 == DiffOp::Trace
                                                 || E1 == DiffOp::NormalTrace || E2 == DiffOp::NormalTrace || E3 == DiffOp::NormalTrace
                                                 || E1 == DiffOp::TraceNormal || E2 == DiffOp::TraceNormal || E3 == DiffOp::TraceNormal); 

                constexpr bool trace_operator_V = ( F1 == DiffOp::Trace       || F2 == DiffOp::Trace       || F3 == DiffOp::Trace
                                                 || F1 == DiffOp::NormalTrace || F2 == DiffOp::NormalTrace || F3 == DiffOp::NormalTrace
                                                 || F1 == DiffOp::TraceNormal || F2 == DiffOp::TraceNormal || F3 == DiffOp::TraceNormal);  

                constexpr bool req_sol_buffer_U = (E1 == DiffOp::Identity || E2 == DiffOp::Identity || E3 == DiffOp::Identity) || trace_operator_U;   

                constexpr bool req_sol_buffer_V = (F1 == DiffOp::Identity || F2 == DiffOp::Identity || F3 == DiffOp::Identity) || trace_operator_V;    

                static_assert(trace_operator_U == trace_operator_V);

                constexpr bool trace_operator = trace_operator_U && trace_operator_V;                             
                ///////////////////////////////////////////////////////////////////////////////////////////////////////

                //Init quadrature and interpolation points ////////////////////////////////////////////////////////////
                const auto & FE_U = aFESpaceU.getFE();
                const auto & FE_V = aFESpaceV.getFE();

                typename FESpaceU::Interpolator FEI_U;
                typename FESpaceV::Interpolator FEI_V;

                Index NLocU = FE_U.getNumPoints();
                Index NLocV = FE_V.getNumPoints();
                Index NLocQuad;

          
                constexpr bool use_FE_Q = (QuadOpt == Other);
         
                if constexpr (QuadOpt == Default)
                {
                    NLocQuad = NLocU;

                    if constexpr(!Same_FESpace)
                    {
                        FEI_V = Interpolator(FE_V,FE_U);
                    }
                }

                if constexpr (QuadOpt == Trial) 
                {
                    NLocQuad = NLocV;

                    if constexpr(!Same_FESpace)
                    {
                        FEI_U = Interpolator(FE_U,FE_V);
                    }
                }

                if constexpr (use_FE_Q)
                {
                    assert(&FE_Q != &_DefaultFE);    

                    NLocQuad = FE_Q.getNumPoints();

                    FEI_U = Interpolator(FE_U,FE_Q);  
                    FEI_V = Interpolator(FE_V,FE_Q); 
                } 

                constexpr bool use_quad_for_U = ((QuadOpt == Trial && !Same_FESpace) || use_FE_Q);
                constexpr bool use_quad_for_V = ((QuadOpt == Default && !Same_FESpace) || use_FE_Q);

                ///////////////////////////////////////////////////////////////////////////////////////////////////////



                //Number of colors of elements ////////////////////////////////////////////////////////////////////////
                Index NumColors;

                if constexpr(ContinuousV)
                    NumColors  = aFESpaceV.getNumColors();
                else 
                    NumColors = 1;
                ///////////////////////////////////////////////////////////////////////////////////////////////////////


                 //Loop over the colors for parallel processing of each element of the same color//////////////////////
                for (Index iColor = 0; iColor < NumColors; iColor++)
                {
#pragma omp parallel
                    {
                        //Initialisation and allocation////////////////////////////////////////////////////////////////

                        // Buffer for storing the local solution and their gradients///////////////
                        std::array<std::vector<Real>,DimU> Buffer_Sol_U;
                        for (Index u = 0; u < DimU; u++) 
                            Buffer_Sol_U[u].resize(NLocU);

                        std::array<std::vector<Real>,DimU> Buffer_Sol_Quad_U;
                        if constexpr(req_sol_buffer_U && use_quad_for_U)
                        {
                            for (Index u = 0; u < DimU; u++) 
                                Buffer_Sol_Quad_U[u].resize(NLocQuad);
                        }

                        std::array<std::vector<Real>,DimU> Buffer_Grad_Quad_U;
                        if constexpr(req_gradient_buffer_U)
                        {
                            for (Index u = 0; u < DimU; u++)  
                                Buffer_Grad_Quad_U[u].resize(NLocQuad * DimD);
                        }

                        std::vector<Real> Buffer_Grad_U;
                        if constexpr(req_gradient_buffer_U && use_quad_for_U)
                        {
                            Buffer_Grad_U.resize(NLocU * DimD);
                        }
                        ///////////////////////////////////////////////////////////////////////////

                        // Buffer for storing the local solution and their gradients///////////////
                        std::array<std::vector<Real>,DimV> Buffer_Sol_V;
                        for (Index v = 0; v < DimV; v++) 
                            Buffer_Sol_V[v].resize(NLocV);

                        std::array<std::vector<Real>,DimV> Buffer_Sol_Quad_V;
                        if constexpr(req_sol_buffer_V && use_quad_for_V)
                        {
                            for (Index v = 0; v < DimV; v++) 
                                Buffer_Sol_Quad_V[v].resize(NLocQuad);
                        }

                        std::array<std::vector<Real>,DimV> Buffer_Grad_Quad_V;
                        if constexpr(req_gradient_buffer_V)
                        {
                            for (Index v = 0; v < DimV; v++)  
                                Buffer_Grad_Quad_V[v].resize(NLocQuad * DimD);
                        }

                        std::vector<Real> Buffer_Grad_V;
                        if constexpr(req_gradient_buffer_V && use_quad_for_V)
                        {
                            Buffer_Grad_V.resize(NLocV * DimD);
                        }
                        ///////////////////////////////////////////////////////////////////////////

                        //Temporary variables for storing global index.
                        std::vector<Index> IdxGlobU(NLocU);
                        std::vector<Index> IdxGlobV;

                        if constexpr( (ContinuousV != ContinuousU) || !Same_FESpace) 
                            IdxGlobV.resize(NLocV);

                        auto Buffer_IJ = InitBuffer<DimI,DimJ>();

                        auto Buffer_UD  = InitBuffer<DimU,DimD>();
                        auto Buffer_FUD = InitBuffer<DimU,DimD>();
                        auto Buffer_U   = InitBuffer<DimU,1>();
                        Real Buffer_UN;

                        auto Buffer_VD = InitBuffer<DimV,DimD>();
                        auto Buffer_V  = InitBuffer<DimV,1>();
                        

                        auto Buffer_IJ_1 = InitBuffer<DimI,DimJ>();
                        auto Buffer_IJ_2 = InitBuffer<DimI,DimJ>();
                        auto Buffer_IJ_3 = InitBuffer<DimI,DimJ>();

                        auto Buffer_VD_1 = InitBuffer<DimV,DimD>();
                        auto Buffer_VD_2 = InitBuffer<DimV,DimD>();
                        auto Buffer_VD_3 = InitBuffer<DimV,DimD>();

                        auto Buffer_V_1  = InitBuffer<DimV,1>();
                        auto Buffer_V_2  = InitBuffer<DimV,1>();
                        auto Buffer_V_3  = InitBuffer<DimV,1>();

                        //Compute the Identity matrix -> to optimize : if not needed should be skipped
                        std::array<std::array<Real,DimD>, DimD> Id; 
                        for (Index d=0;d<DimD;++d)
                        {
                            Id[d].fill(0.0);
                            Id[d][d] = 1.0;        
                        }
                   
                        //Buffer for geometric information at the quadrature point
                        std::array<std::array<Real,DimD>, DimD> Buffer_GradDef;
                        std::array<std::array<Real,DimD>, DimD> Buffer_ComatGradDef;
                        std::array<Real,DimD+1> Buffer_Normal;
                        Real Buffer_Jacobian;
                        Real Buffer_Weight_Jacobian;
                        Real Weight;

                        // Buffer for position informations
                        Index buffer_label;
                        Index buffer_index;
                        RealVector buffer_xyz_quadrature;
                        
                        //Boolean used to check if the input vector is zero
                        bool is_zero = false;

                        //Boolean used to check if the current element is an affine deformation of the reference element 
                        bool is_not_affine = true;

                        //Reecover the number of the elements in each color
                        Index NElemColor; 

                        if constexpr(ContinuousV)
                            NElemColor  = aFESpaceV.getNumElements(iColor);
                        else 
                            NElemColor = aFESpaceV.getNumElements();
                        ///////////////////////////////////////////////////////////////////////////////////////////////
   


                        // Loop on elements////////////////////////////////////////////////////////////////////////////
#pragma omp for
                        for (Index iElt = 0; iElt <NElemColor; iElt++)
                        {

                            // Extracting global numbering of the element 
                            // (here it is assumed that the numbering of elements between spaces is coherent).
                            Index iEltGlob;
                            
                            if constexpr(ContinuousV)
                                iEltGlob = aFESpaceV.EltColorNum2EltGlobalNum(iColor, iElt);
                            else 
                                iEltGlob = iElt;

                            //Pre-compute the global indices.
                            if constexpr(ContinuousU)
                            {
                                aFESpaceU.Loc2Glob(iEltGlob,IdxGlobU);
                            }
                            else
                            {
                                aFESpaceU.Loc2GlobDisc(iEltGlob,IdxGlobU);
                            }

                            if constexpr((ContinuousV != ContinuousU) || !Same_FESpace)
                            {
                                if constexpr(ContinuousV)
                                {
                                    aFESpaceV.Loc2Glob(iEltGlob,IdxGlobV);
                                }
                                else
                                {
                                    aFESpaceV.Loc2GlobDisc(iEltGlob,IdxGlobV);
                                }
                            }
                            
                            //Compute the label of the element if needed
                            if (req_label_buffer)
                                buffer_label = aFESpaceV.getLabel(iEltGlob);
                
                            //boolean indicating if the input vector is zero
                            if constexpr(Check_zero) is_zero = true;
                
                            // Extract the local solution vector 
                            for (Index iLoc = 0; iLoc < NLocU ; iLoc++)
                            {
                                Index Idx = DimU*IdxGlobU[iLoc];
                                
                                for (Index u = 0; u < DimU; u++)
                                {
                                    // Extracting solution.
                                    Real data = dataU[Idx];
                        
                                    Buffer_Sol_U[u][iLoc] = data;
                    
                                    // Checking if the input vector is zero
                                    if constexpr(Check_zero)
                                    {
                                        if (is_zero && fabs(data)>Zero_Machine) {is_zero=false;}
                                    }
                        
                                    // Increasing global index to obtain other dimension.
                                    Idx++;
                                }
                            }


                            //We skip the product if the input vector is zero
                            if constexpr(Check_zero)
                                if (is_zero) continue;

                            if constexpr(req_sol_buffer_U && use_quad_for_U)
                            {
                                for (Index u = 0; u < DimU; u++)
                                    FEI_U.Interpolation(Buffer_Sol_U[u], Buffer_Sol_Quad_U[u]);
                            }

                            if constexpr(req_gradient_buffer_U)
                            {
                                if constexpr(use_quad_for_U)
                                {
                                    for (Index u = 0; u < DimU; u++)
                                    {
                                        //If some hardcoded routines can be used it should be put here instead of the two lines below
                                        FE_U.GradientInterpolation(Buffer_Sol_U[u], Buffer_Grad_U);
                                        FEI_U.Interpolation<DimD>(Buffer_Grad_U, Buffer_Grad_Quad_U[u]);
                                    }
                                }
                                else if constexpr(!use_quad_for_U)
                                {
                                    for (Index u = 0; u < DimU; u++)
                                        FE_U.GradientInterpolation(Buffer_Sol_U[u], Buffer_Grad_Quad_U[u]);
                                }
                            }

                            //Recover geometrical informations when the element is affine
                            if (aFESpaceV.isAffine(iEltGlob))
                            {
                                is_not_affine = false;

                                if constexpr(trace_operator)
                                {
                                    aFESpaceV.getNormal(iEltGlob,0,Buffer_Normal);
                                    Buffer_Jacobian = ArrayAlgebra::Norm(Buffer_Normal);
                                }
                                else
                                {
                                    aFESpaceV.getGradDef(iEltGlob, 0, Buffer_GradDef);
                                    Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                    ArrayAlgebra::CoMat(Buffer_GradDef,Buffer_ComatGradDef);
                                }
                            } 

                            // Loop on quadrature points //////////////////////////////////////////////////////////////
                            ///////////////////////////////////////////////////////////////////////////////////////////
                            ///////////////////////////////////////////////////////////////////////////////////////////

                
                            for (Index iLoc = 0; iLoc < NLocQuad; ++iLoc)
                            {
                                if constexpr(QuadOpt == Default)
                                {
                                    Weight = FE_U.getQuadratureWeight(iLoc);
                                }

                                if constexpr(QuadOpt == Trial)
                                {
                                    Weight = FE_V.getQuadratureWeight(iLoc);
                                }

                                if constexpr(QuadOpt == Other)
                                {
                                    Weight = FE_Q.getQuadratureWeight(iLoc);
                                }
                                
                                //Recover geometric informations
                                if (is_not_affine)
                                {
                                    if constexpr(trace_operator)
                                    {
                                        if constexpr(QuadOpt == Default)
                                        {
                                            aFESpaceU.getNormal(iEltGlob,iLoc,Buffer_Normal);
                                        }

                                        if constexpr(QuadOpt == Trial)
                                        {
                                            aFESpaceV.getNormal(iEltGlob,iLoc,Buffer_Normal);
                                        }

                                        if constexpr(QuadOpt == Other)
                                        {
                                            assert(false);
                                            //???.getNormal(iEltGlob,uvw,Buffer_Normal);
                                        }

                                        Buffer_Jacobian = ArrayAlgebra::Norm(Buffer_Normal);
                                    }
                                    else
                                    {
                                        if constexpr(QuadOpt == Default)
                                        {
                                            aFESpaceU.getGradDef(iEltGlob, iLoc, Buffer_GradDef); 
                                        }

                                        if constexpr(QuadOpt == Trial)
                                        {
                                            aFESpaceV.getGradDef(iEltGlob,iLoc,Buffer_GradDef);
                                        }

                                        if constexpr(QuadOpt == Other)
                                        {
                                            assert(false);
                                            //???.getGradDef(iEltGlob,uvw,Buffer_GradDef);
                                        }

                                        Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                        ArrayAlgebra::CoMat(Buffer_GradDef,Buffer_ComatGradDef);
                                    }
                                }

                                if constexpr(req_xyz_quadrature_buffer)
                                    aFESpaceV.getDoFCoordinate(iEltGlob, iLoc, buffer_xyz_quadrature);

                                if constexpr(req_index_buffer)
                                {
                                    if constexpr(QuadOpt == Default)
                                    {
                                        buffer_index = aFESpaceU.Loc2GlobDisc(iEltGlob,iLoc);
                                    }

                                    if constexpr(QuadOpt == Trial)
                                    {
                                        buffer_index = aFESpaceV.Loc2GlobDisc(iEltGlob,iLoc);
                                    }
                                     
                                    if constexpr(QuadOpt == Other)
                                    {
                                        assert(false); 
                                    }
                                }
                                   
                                if constexpr(req_gradient_buffer_U)
                                {
                                    if constexpr(DimU==1) 
                                    {
                                        // Recovering the gradient of the unkowns at a quadrature point 
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_UD[d] = Buffer_Grad_Quad_U[0][iLoc * DimD + d]; 
                                    }  
                                    else  
                                    {
                                        // Recovering the gradient of the unkowns at a quadrature point 
                                        for (Index u = 0; u < DimU; ++u)
                                            for (Index d = 0; d < DimD; d++)
                                                Buffer_UD[u][d] = Buffer_Grad_Quad_U[u][iLoc * DimD + d];
 
                                    } 
                                }

                                if constexpr(compute_gradient_U)
                                {
                                    // Multiply each gradient by the comatrix of the deformation from the reference element to the current element
                                    if constexpr(DimU==1) 
                                    {                                    
                                        ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_UD,Buffer_FUD);
                                    }  
                                    else  
                                    {
                                        for (Index u = 0; u < DimU; ++u)
                                            ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_UD[u],Buffer_FUD[u]);
                                    } 

                                    ArrayAlgebra::Scale(1.0/Buffer_Jacobian,Buffer_FUD);

                                    if constexpr(E1 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>)
                                        {
                                            Buffer_IJ_1 = Buffer_FUD;
                                        }
                                        else
                                        {
                                            aParameterField_L1.Eval(Buffer_FUD,Buffer_IJ_1,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L2,Field::Identity>)
                                        {
                                            Buffer_IJ_2 = Buffer_FUD;
                                        }
                                        else
                                        {
                                            aParameterField_L2.Eval(Buffer_FUD,Buffer_IJ_2,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E3 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L3,Field::Identity>)
                                        {
                                            Buffer_IJ_3 = Buffer_FUD;
                                        }
                                        else
                                        {
                                            aParameterField_L3.Eval(Buffer_FUD,Buffer_IJ_3,buffer_xyz_quadrature,buffer_label);                          
                                        }
                                    }
                                }       
                                
                                if constexpr(compute_divergence_U)
                                {
                                    Real Buffer = (ArrayAlgebra::ContractionProduct(Buffer_ComatGradDef,Buffer_UD))/Buffer_Jacobian; 

                                    if constexpr(E1 == DiffOp::Divergence)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>)
                                        {
                                            Buffer_IJ_1 = Buffer;
                                        }
                                        else
                                        {
                                            aParameterField_L1.Eval(Buffer,Buffer_IJ_1,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::Divergence)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L2,Field::Identity>)
                                        {
                                            Buffer_IJ_2 = Buffer;
                                        }
                                        else
                                        {
                                            aParameterField_L2.Eval(Buffer,Buffer_IJ_2,buffer_xyz_quadrature,buffer_label);
                                        
                                        }
                                    }

                                    if constexpr(E3 == DiffOp::Divergence)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L3,Field::Identity>)
                                        {
                                            Buffer_IJ_3 = Buffer;
                                        }
                                        else
                                        {
                                            aParameterField_L3.Eval(Buffer,Buffer_IJ_3,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }
                                }          
                                

                  

                                if constexpr(req_sol_buffer_U)
                                {

                                     
                                    // Recovering the unkowns at a quadrature point 
                                    if constexpr(DimU==1) 
                                    {
                                        if constexpr(use_quad_for_U) 
                                        {
                                            Buffer_U = Buffer_Sol_Quad_U[0][iLoc];
                                        }
                                        else
                                        {
                                            Buffer_U = Buffer_Sol_U[0][iLoc];
                                        }
                                    }
                                    else
                                    {
                                        if constexpr(use_quad_for_U) 
                                        {
                                            for (Index u = 0; u < DimU; ++u) 
                                                Buffer_U[u] = Buffer_Sol_Quad_U[u][iLoc];
                                        }
                                        else
                                        {
                                            for (Index u = 0; u < DimU; ++u) 
                                                Buffer_U[u] = Buffer_Sol_U[u][iLoc];
                                        }
                                    }

                                    if constexpr(E1 == DiffOp::Identity || E1 == DiffOp::Trace) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>) 
                                        {
                                            Buffer_IJ_1 = Buffer_U;
                                        }
                                        else
                                        {
                                            Eval<ParameterField_L1::Regularity>(aParameterField_L1,
                                            Buffer_U,Buffer_IJ_1,
                                            buffer_xyz_quadrature,buffer_label,buffer_index);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::Identity || E2 == DiffOp::Trace) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L2,Field::Identity>)
                                        {
                                            Buffer_IJ_2 = Buffer_U;
                                        }
                                        else
                                        {
                                            aParameterField_L2.Eval(Buffer_U,Buffer_IJ_2,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E3 == DiffOp::Identity || E3 == DiffOp::Trace) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L3,Field::Identity>)
                                        {
                                            Buffer_IJ_3 = Buffer_U;
                                        }
                                        else
                                        {
                                            aParameterField_L3.Eval(Buffer_U,Buffer_IJ_3,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E1 == DiffOp::NormalTrace 
                                              || E2 == DiffOp::NormalTrace
                                              || E3 == DiffOp::NormalTrace) 
                                     {
                                        Buffer_UN = 0.0;

                                        for (Index u = 0; u < DimU; ++u) 
                                            Buffer_UN += Buffer_U[u]*Buffer_Normal[u];

                                        Buffer_UN /= Buffer_Jacobian; // Careful -> Is it compatible with non linear operator C or L ? 
                                     }


                                    if constexpr(E1 == DiffOp::NormalTrace) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>) 
                                        {
                                            Buffer_IJ_1 = Buffer_UN; 
                                        }
                                        else
                                        {
                                            aParameterField_L1.Eval(Buffer_UN,Buffer_IJ_1,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::NormalTrace) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L2,Field::Identity>)
                                        {
                                            Buffer_IJ_2 = Buffer_UN;
                                        }
                                        else
                                        {
                                            aParameterField_L2.Eval(Buffer_UN,Buffer_IJ_2,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E3 == DiffOp::NormalTrace) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L3,Field::Identity>)
                                        {
                                            Buffer_IJ_3 = Buffer_UN;
                                        }
                                        else
                                        {
                                            aParameterField_L3.Eval(Buffer_UN,Buffer_IJ_3,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E1 == DiffOp::TraceNormal) 
                                    {                                                
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>) 
                                        {
                                            Buffer_IJ_1 = Buffer_U; 
                                        }
                                        else
                                        {
                                            aParameterField_L1.Eval(Buffer_U,Buffer_IJ_1,Buffer_Normal,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::TraceNormal) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L2,Field::Identity>)
                                        {
                                            Buffer_IJ_2 = Buffer_U;
                                        }
                                        else
                                        {
                                            aParameterField_L2.Eval(Buffer_U,Buffer_IJ_2,Buffer_Normal,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E3 == DiffOp::TraceNormal) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L3,Field::Identity>)
                                        {
                                            Buffer_IJ_3 = Buffer_U;
                                        }
                                        else
                                        {
                                            aParameterField_L3.Eval(Buffer_U,Buffer_IJ_3,Buffer_Normal,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }
                                }

                                if constexpr(E2 != DiffOp::Zero)
                                    ArrayAlgebra::MatAdd(Buffer_IJ_2,Buffer_IJ_1);

                                if constexpr(E3 != DiffOp::Zero)
                                    ArrayAlgebra::MatAdd(Buffer_IJ_3,Buffer_IJ_1);

  

                                if constexpr(std::is_same_v<ParameterField_C,Field::Identity>)
                                {
                                    Buffer_IJ = Buffer_IJ_1;
                                }
                                else
                                {
                                    Eval<ParameterField_C::Regularity>(
                                            aParameterField_C,
                                            Buffer_IJ_1,Buffer_IJ,
                                            buffer_xyz_quadrature,buffer_label,buffer_index);
                                }

                                //Multiply by the quadrature weight
                                ArrayAlgebra::Scale(Weight,Buffer_IJ);

                                if constexpr(req_gradient_buffer_V)
                                {                        
                                    if constexpr(DimV==1) 
                                    {
                                        for (Index d = 0; d < DimD; d++) Buffer_VD[d] = 0.0;
                                    }  
                                    else  
                                    {
                                        for (Index v = 0; v < DimV; ++v)
                                            for (Index d = 0; d < DimD; d++)
                                                Buffer_VD[v][d] = 0.0;
                                    } 

                                    if constexpr(F1 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M1,Field::Identity>)
                                        {
                                            Buffer_VD_1 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M1.AdjointEval(Buffer_IJ,Buffer_VD_1,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::MatAdd(Buffer_VD_1,Buffer_VD);
                                    }


                                    if constexpr(F1 == DiffOp::Divergence)
                                    {
                                        Real Buffer; 
                                        Buffer_VD_1 = Id;

                                        if constexpr(std::is_same_v<ParameterField_M1,Field::Identity>)
                                        {
                                            Buffer = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M1.AdjointEval(Buffer_IJ,Buffer,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::Scale(Buffer,Buffer_VD_1);
                                        ArrayAlgebra::MatAdd(Buffer_VD_1,Buffer_VD);
                                    }


                                    if constexpr(F2 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M2,Field::Identity>)
                                        {
                                            Buffer_VD_2 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M2.AdjointEval(Buffer_IJ,Buffer_VD_2,buffer_xyz_quadrature,buffer_label);                                          
                                        }    

                                        ArrayAlgebra::MatAdd(Buffer_VD_2,Buffer_VD);  
                                    }


                                    if constexpr(F2 == DiffOp::Divergence)
                                    {
                                        Real Buffer; 
                                        Buffer_VD_2 = Id;

                                        if constexpr(std::is_same_v<ParameterField_M2,Field::Identity>)
                                        {
                                            Buffer = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M2.AdjointEval(Buffer_IJ,Buffer,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::Scale(Buffer,Buffer_VD_2);
                                        ArrayAlgebra::MatAdd(Buffer_VD_2,Buffer_VD);
                                    }


                                    if constexpr(F3 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M3,Field::Identity>)
                                        {
                                            Buffer_VD_3 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M3.AdjointEval(Buffer_IJ,Buffer_VD_3,buffer_xyz_quadrature,buffer_label);                                   
                                        }       

                                        ArrayAlgebra::MatAdd(Buffer_VD_3,Buffer_VD);
                                    }

                                    if constexpr(F3 == DiffOp::Divergence)
                                    {
                                        Real Buffer; 
                                        Buffer_VD_3 = Id;

                                        if constexpr(std::is_same_v<ParameterField_M3,Field::Identity>)
                                        {
                                            Buffer = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M3.AdjointEval(Buffer_IJ,Buffer,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::Scale(Buffer,Buffer_VD_3);
                                        ArrayAlgebra::MatAdd(Buffer_VD_3,Buffer_VD);
                                    }

                                    if constexpr(DimV==1) 
                                    {
                                        ArrayAlgebra::TransposeMatMlt(Buffer_ComatGradDef,Buffer_VD);  
                                    }
                                    else  
                                    {
                                       for (Index v = 0; v < DimV; ++v) 
                                             ArrayAlgebra::TransposeMatMlt(Buffer_ComatGradDef,Buffer_VD[v]);  
                                    }

                                    if constexpr(DimV==1) 
                                    {
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_Grad_Quad_V[0][iLoc * DimD + d] = Buffer_VD[d];
                                    }
                                    else
                                    {
                                        for (Index v = 0; v < DimV; ++v)
                                            for (Index d = 0; d < DimD; d++)
                                                Buffer_Grad_Quad_V[v][iLoc * DimD + d] = Buffer_VD[v][d];
                                    }
                                }

                                if constexpr(req_sol_buffer_V)
                                {                        
                                    if constexpr(DimV==1) 
                                    {
                                        Buffer_V = 0.0;
                                    }  
                                    else  
                                    {
                                        for (Index v = 0; v < DimV; ++v) Buffer_V[v] = 0.0;
                                    } 

                                    if constexpr(F1 == DiffOp::Identity || F1 == DiffOp::Trace || F1 == DiffOp::NormalTrace)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M1,Field::Identity>)
                                        {
                                            Buffer_V_1 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M1.AdjointEval(Buffer_IJ,Buffer_V_1,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        if constexpr(F1 == DiffOp::Identity || F1 == DiffOp::Trace) 
                                        {
                                            ArrayAlgebra::MatAdd(Buffer_V_1,Buffer_V);  
                                        }
                                        else if constexpr(F1 == DiffOp::NormalTrace)
                                        {
                                            for (Index v = 0; v < DimV; ++v) Buffer_V[v] += Buffer_V_1*Buffer_Normal[v];
                                        }
                                    }

                                    if constexpr(F1 == DiffOp::TraceNormal)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M1,Field::Identity>)
                                        {
                                            Buffer_V_1 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M1.AdjointEval(Buffer_IJ,Buffer_V_1,Buffer_Normal,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::MatAdd(Buffer_V_1,Buffer_V);  
                                    }

                                    if constexpr(F2 == DiffOp::Identity || F2 == DiffOp::Trace || F2 == DiffOp::NormalTrace)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M2,Field::Identity>)
                                        {
                                            Buffer_V_2 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M2.AdjointEval(Buffer_IJ,Buffer_V_2,buffer_xyz_quadrature,buffer_label);                                          
                                        }    

                                        if constexpr(F2 == DiffOp::Identity || F2 == DiffOp::Trace) 
                                        {
                                            ArrayAlgebra::MatAdd(Buffer_V_2,Buffer_V);  
                                        }
                                        else if constexpr(F2 == DiffOp::NormalTrace)
                                        {
                                            for (Index v = 0; v < DimV; ++v) Buffer_V[v] += Buffer_V_2*Buffer_Normal[v];
                                        }
                                    }

                                    if constexpr(F2 == DiffOp::TraceNormal)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M2,Field::Identity>)
                                        {
                                            Buffer_V_2 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M2.AdjointEval(Buffer_IJ,Buffer_V_2,Buffer_Normal,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::MatAdd(Buffer_V_2,Buffer_V);  
                                    }

                                    if constexpr(F3 == DiffOp::Identity || F3 == DiffOp::Trace || F3 == DiffOp::NormalTrace)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M3,Field::Identity>)
                                        {
                                            Buffer_V_3 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M3.AdjointEval(Buffer_IJ,Buffer_V_3,buffer_xyz_quadrature,buffer_label);                                   
                                        }       

                                        if constexpr(F3 == DiffOp::Identity || F3 == DiffOp::Trace) 
                                        {
                                            ArrayAlgebra::MatAdd(Buffer_V_3,Buffer_V);  
                                        }
                                        else if constexpr(F3 == DiffOp::NormalTrace)
                                        {
                                            for (Index v = 0; v < DimV; ++v) Buffer_V[v] += Buffer_V_3*Buffer_Normal[v];
                                        }
                                    }

                                    if constexpr(F3 == DiffOp::TraceNormal)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M3,Field::Identity>)
                                        {
                                            Buffer_V_3 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M3.AdjointEval(Buffer_IJ,Buffer_V_3,Buffer_Normal,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::MatAdd(Buffer_V_3,Buffer_V);  
                                    }

                                    if constexpr(use_quad_for_V)
                                    {
                                        if constexpr(DimV==1) 
                                        {
                                            Buffer_Sol_Quad_V[0][iLoc] = Buffer_Jacobian*Buffer_V;
                                        }
                                        else
                                        {
                                            for (Index v = 0; v < DimV; ++v)
                                                Buffer_Sol_Quad_V[v][iLoc] = Buffer_Jacobian*Buffer_V[v];
                                        }
                                    }
                                    else
                                    {
                                        if constexpr(DimV==1) 
                                        {
                                            Buffer_Sol_V[0][iLoc] = Buffer_Jacobian*Buffer_V;
                                        }
                                        else
                                        {
                                            for (Index v = 0; v < DimV; ++v)
                                                Buffer_Sol_V[v][iLoc] = Buffer_Jacobian*Buffer_V[v];
                                        }
                                    }
                                }
                            }

                            ///////////////////////////////////////////////////////////////////////////////////////////
                            ///////////////////////////////////////////////////////////////////////////////////////////
                            ///////////////////////////////////////////////////////////////////////////////////////////

                            if constexpr(use_quad_for_V && req_sol_buffer_V)
                            {
                                for (Index v = 0; v < DimV; v++)
                                    FEI_V.TransposeInterpolation(Buffer_Sol_Quad_V[v],Buffer_Sol_V[v]);
                            }

                            if constexpr(req_gradient_buffer_V)
                            {

                                if constexpr(use_quad_for_V)
                                {
                                    for (Index v = 0; v < DimV; v++)
                                    {
                                        //Again hardcoded versions could be used here for efficiency
                                        FEI_V.template TransposeInterpolation<DimD>(Buffer_Grad_Quad_V[v], Buffer_Grad_V);

                                        FE_V.template TransposeGradientInterpolation<!req_sol_buffer_V>(Buffer_Grad_V, Buffer_Sol_V[v]);
                                    }
                                }
                                else
                                {
                                    for (Index v = 0; v < DimV; v++)
                                    {
                                        FE_V.template TransposeGradientInterpolation<!req_sol_buffer_V>(Buffer_Grad_Quad_V[v], Buffer_Sol_V[v]);
                                    }
                                }
                            }

                            // Summing into the global solution vector.
                            for (Index iLoc = 0; iLoc < NLocV; iLoc++)
                            {
                                Index Idx;

                                if constexpr((ContinuousV != ContinuousU) || !Same_FESpace)
                                    Idx =  DimV*IdxGlobV[iLoc];
                                else
                                    Idx =  DimV*IdxGlobU[iLoc];
                                
                                for (Index v = 0; v < DimV; v++)
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

} //Core
 
}//OndoMathX
