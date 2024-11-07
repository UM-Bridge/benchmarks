#pragma once




// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
namespace Core {

template<
    class FESpace,  
    bool ContinuousU, 
    bool ContinuousV,
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
    bool Check_zero = true 
    >
void MltAdd(FESpace & aFESpace,
            ParameterField_C & aParameterField_C,
            Real Alpha,
            LAL::Vector& U,
            Real Beta,
            LAL::Vector& V,
            ParameterField_L1 & aParameterField_L1 = Field::IdentityField,
            ParameterField_M1 & aParameterField_M1 = Field::IdentityField,
            ParameterField_L2 & aParameterField_L2 = Field::IdentityField,
            ParameterField_M2 & aParameterField_M2 = Field::IdentityField,
            ParameterField_L3 & aParameterField_L3 = Field::IdentityField,
            ParameterField_M3 & aParameterField_M3 = Field::IdentityField)
            {
                constexpr Index DimD = FESpace::Dim;

                assert(DimU>0);
                assert(DimV>0);
                assert(DimI>0);
                assert(DimJ>0);
                assert(DimD>1);

                Real * dataU = LAL::getData(U);
                Real * dataV = LAL::getData(V);
                Real Copy_Beta = Beta; 
    
                // Creating temporary variable to store the result of the unassembled "matrix time vector" operation.
                if (Copy_Beta == 0.0)
                {
                    // Setting vector to zero.
                    LAL::Scale(0.0, V);
        
                    Copy_Beta = 1.;
                }
    
                Real Alpha_Beta = Alpha / Copy_Beta;
    
                const GaussLobattoElement & GLE = aFESpace.getFE();
    
                //Boolean used to know if the label of the element is used for the parameters computation
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


                constexpr bool req_sol_buffer_U = (E1 == DiffOp::Identity || E2 == DiffOp::Identity || E3 == DiffOp::Identity);   

                constexpr bool req_sol_buffer_V = (F1 == DiffOp::Identity || F2 == DiffOp::Identity || F3 == DiffOp::Identity);       
              
                //Number of quadrature/interpolation points
                Index NLoc = GLE.getNumPoints();

                Index NumColors;

                if constexpr(ContinuousV)
                    NumColors  = aFESpace.getNumColors();
                else 
                    NumColors = 1;

                for (Index iColor = 0; iColor < NumColors; iColor++)
                {
#pragma omp parallel
                    {
                        // Buffer for storing the local solution.
                        std::array<std::vector<Real>,DimU> Buffer_Sol_U;
                        for (Index u = 0; u < DimU; u++) 
                            Buffer_Sol_U[u].resize(NLoc);

                        
                        // Buffer for storing the local solution.
                        std::array<std::vector<Real>,DimV> Buffer_Sol_V;
                            for (Index v = 0; v < DimV; v++) 
                                Buffer_Sol_V[v].resize(NLoc);

                        // Buffer used for gradients
                        std::array<std::vector<Real>,DimU> Buffer_Grad_U;
                        if constexpr(req_gradient_buffer_U)
                            for (Index u = 0; u < DimU; u++)  
                                Buffer_Grad_U[u].resize(NLoc * DimD);

                        // Buffer used for gradients
                        std::array<std::vector<Real>,DimV> Buffer_Grad_V;
                        if constexpr(req_gradient_buffer_V)
                            for (Index v = 0; v < DimV; v++)  
                                Buffer_Grad_V[v].resize(NLoc * DimD);

                        // Temporary variables for storing global index.
                        std::vector<Index> IdxGlobU(NLoc);
                        std::vector<Index> IdxGlobV;

                        if constexpr(ContinuousV != ContinuousU) 
                            IdxGlobV.resize(NLoc);

                        auto Buffer_IJ = InitBuffer<DimI,DimJ>();
                        auto Buffer_UD = InitBuffer<DimU,DimD>();
                        auto Buffer_VD = InitBuffer<DimV,DimD>();
                        auto Buffer_U  = InitBuffer<DimU,1>();
                        auto Buffer_V  = InitBuffer<DimV,1>();
                        
                        auto Buffer_IJ_1 = InitBuffer<DimI,DimJ>();
                        auto Buffer_IJ_2 = InitBuffer<DimI,DimJ>();
                        auto Buffer_IJ_3 = InitBuffer<DimI,DimJ>();

                        auto Buffer_UD_1 = InitBuffer<DimU,DimD>();

                        auto Buffer_VD_1 = InitBuffer<DimV,DimD>();
                        auto Buffer_VD_2 = InitBuffer<DimV,DimD>();
                        auto Buffer_VD_3 = InitBuffer<DimV,DimD>();

                        auto Buffer_V_1  = InitBuffer<DimV,1>();
                        auto Buffer_V_2  = InitBuffer<DimV,1>();
                        auto Buffer_V_3  = InitBuffer<DimV,1>();

                        std::array<std::array<Real,DimD>, DimD> Id; 
                        for (Index d=0;d<DimD;++d)
                        {
                            Id[d].fill(0.0);
                            Id[d][d] = 1.0;        
                        }
                   
                        std::array<std::array<Real,DimD>, DimD> Buffer_GradDef;
                        std::array<std::array<Real,DimD>, DimD> Buffer_ComatGradDef;
                        Real Buffer_Jacobian;
                        Real Buffer_Weight_Jacobian;

                        // Buffer for position informations
                        Index buffer_label;
                      
                        RealVector buffer_xyz_quadrature;
                    
                        //Boolean used to check if the input vector is zero
                        bool is_zero;

                        //Boolean used to check if the current element is an affine deformation of the reference element 
                        bool is_not_affine = true;

                        //Number of elements with the corresponding color
                        Index NElemColor; 

                        if constexpr(ContinuousV)
                            NElemColor  = aFESpace.getNumElements(iColor);
                        else 
                            NElemColor = aFESpace.getNumElements();

            // Loop on elements.
#pragma omp for
                        for (Index iElt = 0; iElt <NElemColor; iElt++)
                        {
                            // Extracting global numbering of the element.
                            Index iEltGlob;
                            
                            if constexpr(ContinuousV)
                                iEltGlob = aFESpace.EltColorNum2EltGlobalNum(iColor, iElt);
                            else 
                                iEltGlob = iElt;

                            //Pre-compute the global indices.
                            if constexpr(ContinuousU)
                            {
                                aFESpace.Loc2Glob(iEltGlob,IdxGlobU);
                            }
                            else
                            {
                                aFESpace.Loc2GlobDisc(iEltGlob,IdxGlobU);
                            }

                            if constexpr(ContinuousV != ContinuousU)
                            {
                                if constexpr(ContinuousV)
                                {
                                    aFESpace.Loc2Glob(iEltGlob,IdxGlobV);
                                }
                                else
                                {
                                    aFESpace.Loc2GlobDisc(iEltGlob,IdxGlobV);
                                }
                            }
                            
                            //Compute the label of the element if needed
                            if (req_label_buffer)
                                buffer_label = aFESpace.getLabel(iEltGlob);
                
                            //boolean indicating if the input vector is zero
                            if constexpr(Check_zero) is_zero = true;
                
                            // Extract the local solution vector.
                            for (Index iLoc = 0; iLoc < NLoc ; iLoc++)
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

                            for (Index v = 0; v < DimV; v++)
                                for (Index iLoc = 0; iLoc < NLoc; ++iLoc)
                                    Buffer_Sol_V[v][iLoc] = 0.0;
                                
                            //We skip the product if the input vector is zero
                            if constexpr(Check_zero)
                                if (is_zero) continue;

                            if constexpr(req_gradient_buffer_U)
                            {
                                for (Index u = 0; u < DimU; u++)
                                    GLE.GradientInterpolation(Buffer_Sol_U[u], Buffer_Grad_U[u]);
                            }
                            
                            if (aFESpace.isAffine(iEltGlob))
                            {
                                is_not_affine = false;
                                aFESpace.getGradDef(iEltGlob, 0, Buffer_GradDef);
                                Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                ArrayAlgebra::CoMat(Buffer_GradDef,Buffer_ComatGradDef);
                            } 

                            // Loop on quadrature points.
                            for (Index iLoc = 0; iLoc < NLoc; ++iLoc)
                            {
                                Real Weight = GLE.getQuadratureWeight(iLoc);

                                if (is_not_affine)
                                {
                                    aFESpace.getGradDef(iEltGlob, iLoc, Buffer_GradDef); 
                                    Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                    ArrayAlgebra::CoMat(Buffer_GradDef,Buffer_ComatGradDef);
                                }

                                if constexpr(req_xyz_quadrature_buffer)
                                    aFESpace.getDoFCoordinate(iEltGlob, iLoc, buffer_xyz_quadrature);

                                if constexpr(req_gradient_buffer_U)
                                {
                                    if constexpr(DimU==1) 
                                    {
                                        // Recovering the gradient of the unkowns at a quadrature point 
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_UD[d] = Weight*Buffer_Grad_U[0][iLoc * DimD + d]; 
                                    }  
                                    else  
                                    {
                                        // Recovering the gradient of the unkowns at a quadrature point 
                                        for (Index u = 0; u < DimU; ++u)
                                            for (Index d = 0; d < DimD; d++)
                                                Buffer_UD[u][d] = Weight*Buffer_Grad_U[u][iLoc * DimD + d];
 
                                    } 
                                }


                                if constexpr(compute_gradient_U)
                                {
                                    // Multiply each gradient by the comatrix of the deformation from the reference element to the current element
                                    if constexpr(DimU==1) 
                                    {                                    
                                        ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_UD,Buffer_UD_1);
                                    }  
                                    else  
                                    {
                                        for (Index u = 0; u < DimU; ++u)
                                            ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_UD[u],Buffer_UD_1[u]);
                                    } 

                                    if constexpr(E1 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>)
                                        {
                                            Buffer_IJ_1 = Buffer_UD_1;
                                        }
                                        else
                                        {
                                            aParameterField_L1.Eval(Buffer_UD_1,Buffer_IJ_1,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L2,Field::Identity>)
                                        {
                                            Buffer_IJ_2 = Buffer_UD_1;
                                        }
                                        else
                                        {
                                            aParameterField_L2.Eval(Buffer_UD_1,Buffer_IJ_2,buffer_xyz_quadrature,buffer_label);
                                        
                                        }
                                    }

                                    if constexpr(E3 == DiffOp::Gradient)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L3,Field::Identity>)
                                        {
                                            Buffer_IJ_3 = Buffer_UD_1;
                                        }
                                        else
                                        {
                                            aParameterField_L3.Eval(Buffer_UD_1,Buffer_IJ_3,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }
                                }       
                                
                                if constexpr(compute_divergence_U)
                                {
                                    Real Buffer = ArrayAlgebra::ContractionProduct(Buffer_ComatGradDef,Buffer_UD); 

                                    
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
                                    Buffer_Weight_Jacobian = Weight*Buffer_Jacobian;
                                    
                                    // Recovering the unkowns at a quadrature point 
                                    if constexpr(DimU==1) 
                                    {
                                        Buffer_U = Buffer_Weight_Jacobian*Buffer_Sol_U[0][iLoc];
                                    }
                                    else
                                    {
                                        for (Index u = 0; u < DimU; ++u) 
                                            Buffer_U[u] = Buffer_Weight_Jacobian*Buffer_Sol_U[u][iLoc];
                                    }

                                    if constexpr(E1 == DiffOp::Identity) 
                                    {
                                        if constexpr(std::is_same_v<ParameterField_L1,Field::Identity>) 
                                        {
                                            Buffer_IJ_1 = Buffer_U;
                                        }
                                        else
                                        {
                                            aParameterField_L1.Eval(Buffer_U,Buffer_IJ_1,buffer_xyz_quadrature,buffer_label);
                                        }
                                    }

                                    if constexpr(E2 == DiffOp::Identity) 
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

                                    if constexpr(E3 == DiffOp::Identity) 
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
                                    aParameterField_C.Eval(Buffer_IJ_1,Buffer_IJ,buffer_xyz_quadrature,buffer_label);
                                }



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

                                    ArrayAlgebra::Scale(1.0/Buffer_Jacobian,Buffer_VD);

                                    if constexpr(DimV==1) 
                                    {
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_Grad_V[0][iLoc * DimD + d] = Buffer_VD[d];
                                    }
                                    else
                                    {
                                        for (Index v = 0; v < DimV; ++v)
                                            for (Index d = 0; d < DimD; d++)
                                                Buffer_Grad_V[v][iLoc * DimD + d] = Buffer_VD[v][d];
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

                                    if constexpr(F1 == DiffOp::Identity)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M1,Field::Identity>)
                                        {
                                            Buffer_V_1 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M1.AdjointEval(Buffer_IJ,Buffer_V_1,buffer_xyz_quadrature,buffer_label);                                              
                                        }

                                        ArrayAlgebra::MatAdd(Buffer_V_1,Buffer_V);
                                    }

                                    if constexpr(F2 == DiffOp::Identity)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M2,Field::Identity>)
                                        {
                                            Buffer_V_2 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M2.AdjointEval(Buffer_IJ,Buffer_V_2,buffer_xyz_quadrature,buffer_label);                                          
                                        }    

                                        ArrayAlgebra::MatAdd(Buffer_V_2,Buffer_V);  
                                    }

                                    if constexpr(F3 == DiffOp::Identity)
                                    {
                                        if constexpr(std::is_same_v<ParameterField_M3,Field::Identity>)
                                        {
                                            Buffer_V_3 = Buffer_IJ;
                                        }
                                        else
                                        {    
                                            aParameterField_M3.AdjointEval(Buffer_IJ,Buffer_V_3,buffer_xyz_quadrature,buffer_label);                                   
                                        }       

                                        ArrayAlgebra::MatAdd(Buffer_V_3,Buffer_V);
                                    }

                                    if constexpr(DimV==1) 
                                    {
                                        Buffer_Sol_V[0][iLoc] = Buffer_V;
                                    }
                                    else
                                    {
                                        for (Index v = 0; v < DimV; ++v)
                                                Buffer_Sol_V[v][iLoc] = Buffer_V[v];
                                    }
                                }
                            }

                            if constexpr(req_gradient_buffer_V)
                                for (Index v = 0; v < DimV; v++)
                                    GLE.TransposeGradientInterpolation<false>(Buffer_Grad_V[v], Buffer_Sol_V[v]);


                            // Summing into the global solution vector.
                            for (Index iLoc = 0; iLoc < NLoc; iLoc++)
                            {
                                Index Idx;

                                if constexpr(ContinuousV != ContinuousU)
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




template<
    class FESpace,  
    bool ContinuousU, 
    bool ContinuousV,
    Index DimU,
    Index DimV,
    Index DimI,
    Index DimJ,
    class ParameterField_C,
    class SparseMatrix,
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
    class ParameterField_M3 = OndoMathX::Field::Identity
    >
    void Assemble(FESpace & aFESpace,
            ParameterField_C & aParameterField_C,
            SparseMatrix & aMatrix,
            Real Alpha,
            ParameterField_L1 & aParameterField_L1 = Field::IdentityField,
            ParameterField_M1 & aParameterField_M1 = Field::IdentityField,
            ParameterField_L2 & aParameterField_L2 = Field::IdentityField,
            ParameterField_M2 & aParameterField_M2 = Field::IdentityField,
            ParameterField_L3 & aParameterField_L3 = Field::IdentityField,
            ParameterField_M3 & aParameterField_M3 = Field::IdentityField)
            {
                constexpr bool lumped_matrix = 
                             (   E1 == Identity
                            &&   F1 == Identity
                            &&   (E2 == Zero || E2 == Identity) 
                            &&   (F2 == Zero || F2 == Identity) 
                            &&   (E3 == Zero || E3 == Identity) 
                            &&   (F3 == Zero || F3 == Identity));

                //Todo : use the knowledge that the matrix is lumped to simplify some loops

                const GaussLobattoElement & GLE = aFESpace.getFE();

                Index NLoc = GLE.getNumPoints();
                Index NElem = aFESpace.getNumElements();

                LAL::Vector U;
                LAL::Vector V;

                Index nU = aFESpace.getNumDoFsDisc()*DimU;
                Index nV = aFESpace.getNumDoFsDisc()*DimV; 

                LAL::Allocate(U,nU);
                LAL::Allocate(V,nV);

                Real * pU = LAL::getData(U);
                Real * pV = LAL::getData(V);

                // Loop on the dimension of the unknowns
                for (Index u=0; u<DimU; ++u)
                {
                    // Loop on local number of dof
                    for (Index iLoc=0; iLoc<NLoc; ++iLoc)   
                    {
                        {
                            //Loop on the elements to assign 1 to the corresponding dof
                            for (Index iElt=0;iElt<NElem;++iElt)
                            {
                                Index IdxU = aFESpace.Loc2GlobDisc(iElt,iLoc)*DimU+u;
                                pU[IdxU] = Alpha;
                            }
                            
                            //MltAdd in order to recover one row of the elementary matrix of each elements
                            MltAdd<FESpace,false,false,
                            DimU,
                            DimV,
                            DimI,
                            DimJ,
                            ParameterField_C,
                            E1,F1,ParameterField_L1,ParameterField_M1,
                            E2,F2,ParameterField_L2,ParameterField_M2,
                            E3,F3,ParameterField_L3,ParameterField_M3> (aFESpace, aParameterField_C, 1.0, U, 0.0, V,
                            aParameterField_L1, aParameterField_M1,
                            aParameterField_L2, aParameterField_M2,
                            aParameterField_L3, aParameterField_M3);
                            
                

                            //Loop on the elements and assembling procedure
                            for (Index iElt=0;iElt<NElem;++iElt)
                            {
                                Index IdxU;

                                //Pre-compute the global indices.
                                if constexpr(ContinuousU)
                                {
                                    IdxU = aFESpace.Loc2Glob(iElt,iLoc)*DimU+u;
                                }
                                else
                                {
                                    IdxU = aFESpace.Loc2GlobDisc(iElt,iLoc)*DimU+u;
                                }

                                for (Index jLoc=0; jLoc<NLoc; ++jLoc)
                                {
                                    for (Index v=0; v<DimV; ++v)
                                    {
                                        Index IdxV;
                                        Index IdxV_D;
                                        
                                        IdxV_D = aFESpace.Loc2GlobDisc(iElt,jLoc)*DimV+v;
                                        
                                        if constexpr(ContinuousV == true)
                                            IdxV = aFESpace.Loc2Glob(iElt,jLoc)*DimV+v;
                                        else
                                            IdxV = IdxV_D;
                                        
                                        Real value = pV[IdxV_D];
                                        
                                        if (value != 0.0)
                                        {
                                            LAL::AddInteraction(aMatrix,IdxV,IdxU,value);
                                        }
                                    }
                                }
                            }
                            
                            //Loop on the elements to assign 1 to the corresponding dof
                            for (Index iElt=0;iElt<NElem;++iElt)
                            {
                                Index IdxU = aFESpace.Loc2GlobDisc(iElt,iLoc)*DimU+u;
                                pU[IdxU] = 0.0;
                            }
                        }
                    }
                }

            }

template<   
    Index SysDim = 1,
    Index NQuadNL = 1,
    class FESpace,  
    class ParameterField_C
    >
void NonLinearStiffnessMltAdd(FESpace & aFESpace,
            ParameterField_C & aParameterField_C,
            Real Alpha,
            LAL::Vector& U,
            LAL::Vector& W,
            Real Beta,
            LAL::Vector& V,
            Real dt)
            {
                constexpr Index DimD = FESpace::Dim;

                assert(SysDim>0);
                assert(DimD>1);

                Real * dataU = LAL::getData(U);
                Real * dataV = LAL::getData(V);
                Real * dataW = LAL::getData(W);
                Real Copy_Beta = Beta; 
    
                // Creating temporary variable to store the result of the unassembled "matrix time vector" operation.
                if (Copy_Beta == 0.0)
                {
                    // Setting vector to zero.
                    LAL::Scale(0.0, V);
        
                    Copy_Beta = 1.;
                }
    
                Real Alpha_Beta = Alpha / Copy_Beta;
    
                const GaussLobattoElement & GLE = aFESpace.getFE();
    
                //Boolean used to know if the label of the element is used for the parameters computation
                constexpr bool req_label_buffer = 
                             (   ParameterField_C::Regularity == Field::Regularity::PiecewiseC0) || 
                             (   ParameterField_C::Regularity == Field::Regularity::PiecewiseConstant);

                constexpr bool req_xyz_quadrature_buffer = 
                             (   ParameterField_C::Regularity == Field::Regularity::C0) || 
                             (   ParameterField_C::Regularity == Field::Regularity::PiecewiseC0); 

                //Number of quadrature/interpolation points
                Index NLoc = GLE.getNumPoints();

                Index NumColors  = aFESpace.getNumColors();
  
                for (Index iColor = 0; iColor < NumColors; iColor++)
                {
#pragma omp parallel
                    {
                        // Buffer for storing the local solution.
                        std::array<std::vector<Real>,SysDim> Buffer_Sol_U;
                        for (Index u = 0; u < SysDim; u++) 
                            Buffer_Sol_U[u].resize(NLoc);

                        std::array<std::vector<Real>,SysDim> Buffer_Sol_W;
                        for (Index u = 0; u < SysDim; u++) 
                            Buffer_Sol_W[u].resize(NLoc);

                        std::array<std::vector<Real>,SysDim> Buffer_Sol_V;
                            for (Index v = 0; v < SysDim; v++) 
                                Buffer_Sol_V[v].resize(NLoc);

                        // Buffer used for gradients
                        std::array<std::vector<Real>,SysDim> Buffer_Grad_U;
                        for (Index u = 0; u < SysDim; u++)  
                            Buffer_Grad_U[u].resize(NLoc * DimD);

                        std::array<std::vector<Real>,SysDim> Buffer_Grad_W;
                        for (Index u = 0; u < SysDim; u++)  
                            Buffer_Grad_W[u].resize(NLoc * DimD);

                        std::array<std::vector<Real>,SysDim> Buffer_Grad_V;
                        for (Index v = 0; v < SysDim; v++)  
                            Buffer_Grad_V[v].resize(NLoc * DimD);

                        // Temporary variables for storing global index.
                        std::vector<Index> IdxGlobU(NLoc);
                        

                        auto Buffer_IJ     = InitBuffer<SysDim,DimD>();
                        auto Buffer_IJ_ACC = InitBuffer<SysDim,DimD>();
                        auto Buffer_UD     = InitBuffer<SysDim,DimD>();
                        auto Buffer_WD     = InitBuffer<SysDim,DimD>();
                        auto Buffer_VD     = InitBuffer<SysDim,DimD>();
                        auto Buffer_UWD    = InitBuffer<SysDim,DimD>();
                   
                        std::array<std::array<Real,DimD>, DimD> Buffer_GradDef;
                        std::array<std::array<Real,DimD>, DimD> Buffer_ComatGradDef;
                        Real Buffer_Jacobian;
                        Real Buffer_Weight_Jacobian;

                        // Buffer for position informations
                        Index buffer_label;                      
                        RealVector buffer_xyz_quadrature;
                
                        //Boolean used to check if the current element is an affine deformation of the reference element 
                        bool is_not_affine = true;

                        //Number of elements with the corresponding color
                        Index NElemColor; 

                        NElemColor  = aFESpace.getNumElements(iColor);
                 

            // Loop on elements.
#pragma omp for
                        for (Index iElt = 0; iElt <NElemColor; iElt++)
                        {
                            // Extracting global numbering of the element.
                            Index iEltGlob;
                            
                            iEltGlob = aFESpace.EltColorNum2EltGlobalNum(iColor, iElt);
                           
                            //Pre-compute the global indices.
                            aFESpace.Loc2Glob(iEltGlob,IdxGlobU);
                           
                            //Compute the label of the element if needed
                            if (req_label_buffer)
                                buffer_label = aFESpace.getLabel(iEltGlob);
                
                          
                            // Extract the local solution vector.
                            for (Index iLoc = 0; iLoc < NLoc ; iLoc++)
                            {
                                Index Idx = SysDim*IdxGlobU[iLoc];
                                
                                for (Index u = 0; u < SysDim; u++)
                                {
                                    // Extracting solution. 
                                    Buffer_Sol_U[u][iLoc] = dataU[Idx];
                                    if constexpr(NQuadNL != 1) Buffer_Sol_W[u][iLoc] = dataW[Idx];
                    
                                    // Increasing global index to obtain other dimension.
                                    Idx++;
                                }
                            }

                            for (Index v = 0; v < SysDim; v++)
                                for (Index iLoc = 0; iLoc < NLoc; ++iLoc)
                                    Buffer_Sol_V[v][iLoc] = 0.0;
                                
                            for (Index u = 0; u < SysDim; u++)
                            {
                                GLE.GradientInterpolation(Buffer_Sol_U[u], Buffer_Grad_U[u]);
                                if constexpr(NQuadNL != 1) GLE.GradientInterpolation(Buffer_Sol_W[u], Buffer_Grad_W[u]);
                            }
                          
                            
                            if (aFESpace.isAffine(iEltGlob))
                            {
                                is_not_affine = false;
                                aFESpace.getGradDef(iEltGlob, 0, Buffer_GradDef);
                                Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                ArrayAlgebra::CoMat(Buffer_GradDef,Buffer_ComatGradDef);
                            } 

                            // Loop on quadrature points.
                            for (Index iLoc = 0; iLoc < NLoc; ++iLoc)
                            {
                                Real Weight = GLE.getQuadratureWeight(iLoc);

                                if (is_not_affine)
                                {
                                    aFESpace.getGradDef(iEltGlob, iLoc, Buffer_GradDef); 
                                    Buffer_Jacobian = ArrayAlgebra::Det(Buffer_GradDef);
                                    ArrayAlgebra::CoMat(Buffer_GradDef,Buffer_ComatGradDef);
                                }

                                if constexpr(req_xyz_quadrature_buffer)
                                    aFESpace.getDoFCoordinate(iEltGlob, iLoc, buffer_xyz_quadrature);


                                if constexpr(SysDim==1) 
                                {
                                    // Recovering the gradient of the unkowns at a quadrature point 
                                    for (Index d = 0; d < DimD; d++)
                                        Buffer_UD[d] = Buffer_Grad_U[0][iLoc * DimD + d]; 

                                    if constexpr(NQuadNL != 1)
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_WD[d] = Buffer_Grad_W[0][iLoc * DimD + d]; 
                                }  
                                else  
                                {
                                    // Recovering the gradient of the unkowns at a quadrature point 
                                    for (Index u = 0; u < SysDim; ++u)
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_UD[u][d] = Buffer_Grad_U[u][iLoc * DimD + d];

                                    if constexpr(NQuadNL != 1)
                                        for (Index u = 0; u < SysDim; ++u)
                                            for (Index d = 0; d < DimD; d++)
                                                Buffer_WD[u][d] = Buffer_Grad_W[u][iLoc * DimD + d];

                                } 
                                
                                // Multiply each gradient by the comatrix of the deformation from the reference element to the current element
                                if constexpr(SysDim==1) 
                                {                                    
                                    ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_UD);
                                    ArrayAlgebra::Scale(1.0/Buffer_Jacobian,Buffer_UD);

                                    if constexpr(NQuadNL != 1)
                                    {
                                        ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_WD);
                                        ArrayAlgebra::Scale(1.0/Buffer_Jacobian,Buffer_WD);
                                    }
                                }  
                                else  
                                {
                                    for (Index u = 0; u < SysDim; ++u)
                                    {
                                        ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_UD[u]);
                                        ArrayAlgebra::Scale(1.0/Buffer_Jacobian,Buffer_UD[u]);
                                    }
                                    
                                    if constexpr(NQuadNL != 1)
                                    {
                                        for (Index u = 0; u < SysDim; ++u)
                                        {
                                            ArrayAlgebra::MatMlt(Buffer_ComatGradDef,Buffer_WD[u]);
                                            ArrayAlgebra::Scale(1.0/Buffer_Jacobian,Buffer_WD[u]);
                                        }
                                    }
                                } 
                      

                                if constexpr(NQuadNL==1) 
                                {
                                    aParameterField_C.Eval(Buffer_UD,Buffer_VD,buffer_xyz_quadrature,buffer_label);
                                    ArrayAlgebra::Scale(dt,Buffer_VD);
                                }

                                if constexpr(NQuadNL==2) 
                                {
                                    //Second order gauss rule
                                    ArrayAlgebra::Scale(0.0,Buffer_VD); 

                                    Buffer_UWD = Buffer_UD;
                                    ArrayAlgebra::MatScaleAdd(-0.5*dt/sqrt(3),Buffer_WD,1.0,Buffer_UWD);
                                    aParameterField_C.Eval(Buffer_UWD,Buffer_IJ,buffer_xyz_quadrature,buffer_label);
                                    ArrayAlgebra::MatScaleAdd(0.5*dt,Buffer_IJ,1.0,Buffer_VD); 

                                    Buffer_UWD = Buffer_UD;
                                    ArrayAlgebra::MatScaleAdd(+0.5*dt/sqrt(3),Buffer_WD,1.0,Buffer_UWD);
                                    aParameterField_C.Eval(Buffer_UWD,Buffer_IJ,buffer_xyz_quadrature,buffer_label);
                                    ArrayAlgebra::MatScaleAdd(0.5*dt,Buffer_IJ,1.0,Buffer_VD);    
                                }

                                if constexpr(SysDim==1) 
                                {
                                    ArrayAlgebra::TransposeMatMlt(Buffer_ComatGradDef,Buffer_VD);  
                                }
                                else  
                                {
                                    for (Index v = 0; v < SysDim; ++v) 
                                            ArrayAlgebra::TransposeMatMlt(Buffer_ComatGradDef,Buffer_VD[v]);  
                                }

                                ArrayAlgebra::Scale(Weight,Buffer_VD);

                                if constexpr(SysDim==1) 
                                {
                                    for (Index d = 0; d < DimD; d++)
                                        Buffer_Grad_V[0][iLoc * DimD + d] = Weight*Buffer_VD[d];
                                }
                                else
                                {
                                    for (Index v = 0; v < SysDim; ++v)
                                        for (Index d = 0; d < DimD; d++)
                                            Buffer_Grad_V[v][iLoc * DimD + d] = Weight*Buffer_VD[v][d];
                                }
                              
                            }

                            for (Index v = 0; v < SysDim; v++)
                                GLE.TransposeGradientInterpolation<false>(Buffer_Grad_V[v], Buffer_Sol_V[v]);


                            // Summing into the global solution vector.
                            for (Index iLoc = 0; iLoc < NLoc; iLoc++)
                            {
                                Index Idx;

                                Idx =  SysDim*IdxGlobU[iLoc];
                                
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




} //Core
 
}//OndoMathX
