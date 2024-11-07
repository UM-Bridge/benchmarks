
namespace OndoMathX {
 
TetriXP2Space::TetriXP2Space(const Mesh & aMesh) : _mesh(aMesh)
{
   //const auto& quad_rule_list = QuadratureNS::GetTetrahedronQuadratureRules();
    //const auto& quad_rule = *(quad_rule_list[QuadratureNS::Type::Tet15]);
    
    assert(_mesh.getDim() == 3);
    
    _elements = _mesh.getElementList(3);
    _NElem = _elements.size();
 
    std::vector<RealVector> DofSet;
    std::vector<std::vector<Index>> all_loc2global;
    
    DofSet.resize(_NLocDoF);
    for (Index i=0; i<DofSet.size();i++)
    {
       // DofSet[i] = quad_rule.GetPoint(i).GetCoordinates();
       DofSet[i] = TetriX::Tet_15::Point[i];
    }

    _NDoF = Numbering<RTetrahedron>(_mesh,DofSet,all_loc2global);
    
    //Construct the array barycenter, the jacobian at each element points and the deformation gradient as well
    _Barycenter.resize(_NElem);
    _loc2glob.resize(_NElem);
    _GradDef_Origin.resize(_NElem);
    _ComatGradDef_Origin.resize(_NElem);
    _Jacobian_Origin.resize(_NElem);
    _SymGradDef_Origin.resize(_NElem);
    
    _Affine.resize(_NElem,true);

    //Construct the array storing the coordinates of the dof
    _DofCoordinates.resize(_NDoF);


    std::vector<Index> glob2glob_P2(_NDoF,0);
    _P2_NDoF = 0;
       
    for (Index el=0;el<_NElem;++el)
    {
        Index elem = _elements[el];

        assert(_mesh.getElement(elem).getType() == GeoElement::GeoTetrahedron4);

        _loc2glob[el] = all_loc2global[elem];
   
        //Get the barycenter
        _Barycenter[el]  = _mesh.getElement(elem).getBarycenter();

        RealVector uvw = {0.0,0.0,0.0};
        
        _mesh.getElement(elem).getGradDef(uvw,_GradDef_Origin[el]);

        ArrayAlgebra::CoMat(_GradDef_Origin[el],_ComatGradDef_Origin[el]);
        _Jacobian_Origin[el] = ArrayAlgebra::Det(_GradDef_Origin[el]);
        _SymGradDef_Origin[el]= _ComatGradDef_Origin[el]; 
        ArrayAlgebra::TransposeMatMlt(_ComatGradDef_Origin[el],_SymGradDef_Origin[el]);
        ArrayAlgebra::Scale(1.0/_Jacobian_Origin[el],_SymGradDef_Origin[el]);


        for (Index i=0;i<_P2_NLocDoF;i++) // P2 dof
        {
            Index globalDof = _loc2glob[el][i];

            if (glob2glob_P2[globalDof] == 0)
            {
                _P2_NDoF++;
                glob2glob_P2[globalDof] = _P2_NDoF;
            }
        }

        for (Index i=0;i<_NLocDoF;i++)
        {
            Index globalDof = _loc2glob[el][i];
            
            //Recover the position of the Dof
            //RealVector uvw = quad_rule.GetPoint(i).GetCoordinates();
            RealVector uvw = TetriX::Tet_15::Point[i];
            
            //Recover the coordinates of the dof (same as quadrature because we are doing specrtal finite elements)
            _mesh.getElement(elem).uvw2xyz(uvw,_DofCoordinates[globalDof]);
        }
    }


    //Correspondance with P2 finite elements
    _P2_glob2glob.resize(_P2_NDoF);

    for (Index i=0;i<_NDoF;++i)
    {
        if (glob2glob_P2[i] != 0)
            _P2_glob2glob[glob2glob_P2[i]-1] = i;
    }

    _P2_loc2glob.resize(_NElem);
    
    for (Index el=0;el<_NElem;++el)
    {
        _P2_loc2glob[el].resize(_P2_NLocDoF);

        for (Index i=0;i<_P2_NLocDoF;i++)
        {
            _P2_loc2glob[el][i] = glob2glob_P2[_loc2glob[el][i]]-1;
        }
    }


    
    // //Deal with the boundary
    // std::vector<Index> labels = _mesh->getLabels(DIM-1);
    
    // for (Index l:labels)
    // {
    //     std::vector<Index> elem_list = _mesh->getElementList(DIM-1,l);
        
    //     _mesh_bndy[l] =  std::make_shared<Mesh>(_mesh->extract(elem_list));

    //     _Jacobians_bndy[l].resize(elem_list.size());
    //     _Barycenter_bndy[l].resize(elem_list.size());
    //     _loc2glob_bndy[l].resize(elem_list.size());
        
    //     for (int el=0;el<elem_list.size();++el)
    //     {
    //         Index elem = elem_list[el];
    //         Index NLocDoF = all_loc2global[elem].size();
    //         _loc2glob_bndy[l][el] = all_loc2global[elem];
    //         _Jacobians_bndy[l][el].resize(NLocDoF);
    //         _Barycenter_bndy[l][el] = _mesh->getElement(elem).getBarycenter();
            
    //         for (Index i=0;i<NLocDoF;i++)
    //         {
    //             //Recover the position of the Dof
    //             Point uvw = _GLE_bndy.GetPointCoordinate(i);
            
    //             //Recover the jacobian
    //             std::array<std::array<Real, DIM>, DIM > F;
    //             std::array<Real, DIM> S;
                
    //             _mesh->getElement(elem).getJacobian(uvw,F);
                
    //             if (DIM == 2) _Jacobians_bndy[l][el][i]
    //                 = sqrt(F[0][0]*F[0][0]+F[1][0]*F[1][0]);
    //             else
    //             {
    //                 //To check !!!!!!!!!!!!!!!!!!!!! TODO
    //                 //assert(false);
    //                 //Should be wrong, use the transpose of F instead
    //                 ArrayAlgebra::CrossProduct(F[0],F[1],S);
    //                 _Jacobians_bndy[l][el][i] = fabs(ArrayAlgebra::Norm(S));
    //             }
    //         }
    //     }
        
    // }
    
    //Construct the coloring of the mesh and construct the array
    std::vector<Index> colors;
    getColoring(_mesh,colors,_NColor);

    _elemByColor.resize(_NColor);
    
    for (Index el=0;el<_NElem;el++)
    {
        Index elem = _elements[el];
        Index color = colors[elem];
        
        _elemByColor[color].push_back(el);
    }
    
    //Init the array glob2elem
    // _glob2elem.resize(_NDoF);
    
    // for (int el=0;el<_NElem;++el)
    // {
    //     for (int i=0;i<_loc2glob[el].size();++i)
    //     {
    //         _glob2elem[_loc2glob[el][i]].push_back(el);
    //     }
    // }
    
    // //Sort and remove duplicate
    // for (int el=0;el<_NElem;++el)
    // {
    //     std::sort(_glob2elem[el].begin(),_glob2elem[el].end());
    //     auto last = std::unique(_glob2elem[el].begin(),_glob2elem[el].end());
    //     _glob2elem[el].erase(last, _glob2elem[el].end());
    // }
    
}

void TetriXP2Space::AssembleMass(LAL::DiagonalMatrix &mass, Real alpha) const
{
    //Index SysDim = LAL::getDimension(mass)/_NDoF;
    Index SysDim = 1;

    //const auto& quad_rule_list = QuadratureNS::GetTetrahedronQuadratureRules();
    //const auto& quad_rule = *(quad_rule_list[QuadratureNS::Type::Tet15]);

    for (Index iElt = 0; iElt < _NElem; iElt++)
    {
        for (Index iLoc = 0; iLoc < _NLocDoF ; iLoc++)
        {
            Index iGlob = _loc2glob[iElt][iLoc];
        
            Real J = ArrayAlgebra::Det(_GradDef_Origin[iElt]);

            //Real weight = quad_rule.GetPoint(iLoc).GetWeight();
            Real weight = TetriX::Tet_15::Weight[iLoc];

            LAL::AddInteraction(mass,
                                iGlob,
                                iGlob,
                                alpha*J*weight);
        }
    }  
}



void TetriXP2Space::MltStiffness(const LAL::Vector &U, LAL::Vector &V, Real alpha) const
{
    //const auto& quad_rule_list = QuadratureNS::GetTetrahedronQuadratureRules();
    //const auto& quad_rule = *(quad_rule_list[QuadratureNS::Type::Tet15]);

    LAL::Scale(0.0, V);

    const Real * dataU = LAL::getData(U);
    Real * dataV = LAL::getData(V);

    constexpr Index SysDim = 1;
    constexpr Index NQuadPoints = 14;

    for (Index iColor = 0; iColor < _NColor; iColor++)
    {
#pragma omp parallel
    {
        // // Buffer for storing the local solution.
        // std::array<std::vector<Real>,SysDim> Buffer_Sol_U;
        // for (Index u = 0; u < SysDim; u++) 
        //     Buffer_Sol_U[u].resize(_NLocDoF);
 
        // std::array<std::vector<Real>,SysDim> Buffer_Sol_V;
        // for (Index u = 0; u < SysDim; u++) 
        //     Buffer_Sol_V[u].resize(_NLocDoF);
 
        // // Buffer used for gradients
        // std::array<std::vector<Real>,SysDim> Buffer_Grad_U;
        // for (Index u = 0; u < SysDim; u++)  
        //     Buffer_Grad_U[u].resize(NQuadPoint * 3);


        // Buffer for storing the local solution.
        Real Buffer_Sol_U[_NLocDoF];
        Real Buffer_Sol_V[_NLocDoF];
   
        // Buffer used for gradients
        Real Buffer_Grad_U[NQuadPoints * 3];
        Real Buffer_Grad_V[NQuadPoints * 3];

        RealVector Buffer_UD;
   
        Index NElemColor  = _elemByColor[iColor].size();

        // Loop on elements.
#pragma omp for
        for (Index iElt = 0; iElt <NElemColor; iElt++)
        {
            Index iEltGlob = _elemByColor[iColor][iElt];
 
            const std::vector<Index> &Loc2Glob = _loc2glob[iEltGlob];

            // Extract the local solution vector.
            for (Index iLoc = 0; iLoc < _NLocDoF ; iLoc++)
            {
                //Index Idx = SysDim*Loc2Glob[iLoc];
                
                Index Idx = Loc2Glob[iLoc];

                Buffer_Sol_U[iLoc] = dataU[Idx];

             
                // for (Index u = 0; u < SysDim; u++)
                // {
                //     // Extracting solution. 
                //     Buffer_Sol_U[u][iLoc] = dataU[Idx];
        
                //     // Increasing global index to obtain other dimension.
                //     Idx++;
                // }
            }

            // for (Index u = 0; u < SysDim; u++)
            // {
            //     GradientInterpolation(Buffer_Sol_U[u], Buffer_Grad_U[u]);
            // }

            TetriX::Tet_15_14::GradientInterpolation_Opt(Buffer_Sol_U, Buffer_Grad_U);

            const RealMatrix3x3 & SymGrad = _SymGradDef_Origin[iEltGlob];


            for (Index iQuad = 0; iQuad < NQuadPoints; iQuad++)
            {
                //Real weight = quad_rule.GetPoint(iQuad).GetWeight();

                for (Index d = 0; d < 3; d++)
                    Buffer_UD[d] = Buffer_Grad_U[iQuad * 3 + d]; 

                ArrayAlgebra::MatMlt(SymGrad,Buffer_UD);

                Real weight = TetriX::Tet_14::Weight[iQuad];

                for (Index d = 0; d < 3; d++)
                    Buffer_Grad_V[iQuad * 3 + d] = weight*Buffer_UD[d]; 
            }
            
            TetriX::Tet_15_14::TransposeGradientInterpolation(Buffer_Grad_V, Buffer_Sol_V);

/*
            Real Buffer_Jacobian = ArrayAlgebra::Det(_GradDef_Origin[iEltGlob]);
            for (Index iQuad = 0; iQuad < 15; iQuad++)
            {
                Real weight = TetriX::Tet_15::Weight[iQuad];
                Buffer_Sol_V[iQuad] =  weight*Buffer_Jacobian*Buffer_Sol_U[iQuad];
            }*/

            // Summing into the global solution vector.
            for (Index iLoc = 0; iLoc < _NLocDoF; iLoc++)
            {
                Index Idx = Loc2Glob[iLoc];

                dataV[Idx] += alpha*Buffer_Sol_V[iLoc];

               

                // Index Idx;

                // Idx =  SysDim*IdxGlobU[iLoc];
                
                // for (Index v = 0; v < SysDim; v++)
                // {
                //     // Adding local contribution.
                //     dataV[Idx] += Alpha_Beta * Buffer_Sol_V[v][iLoc];
        
                //     // Updating global index.
                //     Idx++;
                // }
            }


        }




    }




    }


}


} 




