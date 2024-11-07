
namespace OndoMathX {
 
template<Index Dim> FEMesh<Dim>::FEMesh(const Mesh & aMesh, Index order) : _mesh(aMesh)
{
 
    assert(Dim == _mesh.getDim());
    _Order = order;
    
    _elements = _mesh.getElementList(Dim);
    _NElem = _elements.size();

    std::vector<RealVector> DofSet;
    std::vector<std::vector<Index>> all_loc2global;
    
    if (Dim == 2)
    {
        _GLE = GaussLobattoElement(_Order+1,_Order+1,0);
        //_GLE_bndy = GaussGaussLobattoElement(_Order+1,0,0);
    }
    else if (Dim == 3)
    {
        _GLE = GaussLobattoElement(_Order+1,_Order+1,_Order+1);
        //_GLE_bndy = GaussGaussLobattoElement(_Order+1,_Order+1,0);
    }
    else {assert(false);}

    _NLocDoF = _GLE.getNumPoints();
    
    DofSet.resize(_GLE.getNumPoints());
    for (Index i=0; i<DofSet.size();i++)
        DofSet[i] =  _GLE.getPointCoordinate(i);


    if (Dim == 2) _NDoF = Numbering<RQuadrangle>(_mesh,DofSet,all_loc2global);
    else if (Dim == 3) _NDoF = Numbering<RHexahedron>(_mesh,DofSet,all_loc2global);
  
    //Construct the array barycenter, the jacobian at each quadrature points, the coordinates of the dof and the deformation gradient as well
    _Barycenter.resize(_NElem);
    FESpaceT<FEMesh<Dim>,Dim>::_loc2glob.resize(_NElem);
    _DofCoordinates.resize(_NDoF);
   
    
    for (Index el=0;el<_NElem;++el)
    {
        Index elem = _elements[el];
        FESpaceT<FEMesh<Dim>,Dim>::_loc2glob[el] = all_loc2global[elem];
   
        //Get the barycenter
        _Barycenter[el]  = _mesh.getElement(elem).getBarycenter();
       
        for (Index i=0;i<_NLocDoF;i++)
        {
            Index globalDof = FESpaceT<FEMesh<Dim>,Dim>::_loc2glob[el][i];
            
            //Recover the position of the Dof
            RealVector uvw = _GLE.getPointCoordinate(i);    
            
            //Recover the coordinates of the dof (same as quadrature because we are doing specrtal finite elements)
            _mesh.getElement(elem).uvw2xyz(uvw,_DofCoordinates[globalDof]);
        }
    }

    FESpaceT<FEMesh<Dim>,Dim>::_Loc2GlobPreComputed = true;
    
    /*
    //Deal with the boundary
    std::vector<Index> labels = _mesh->getLabels(DIM-1);
    
    for (Index l:labels)
    {
        std::vector<Index> elem_list = _mesh->getElementList(DIM-1,l);
        
        _mesh_bndy[l] =  std::make_shared<Mesh>(_mesh->extract(elem_list));

        _Jacobians_bndy[l].resize(elem_list.size());
        _Barycenter_bndy[l].resize(elem_list.size());
        _loc2glob_bndy[l].resize(elem_list.size());
        
        for (int el=0;el<elem_list.size();++el)
        {
            Index elem = elem_list[el];
            Index NLocDoF = all_loc2global[elem].size();
            _loc2glob_bndy[l][el] = all_loc2global[elem];
            _Jacobians_bndy[l][el].resize(NLocDoF);
            _Barycenter_bndy[l][el] = _mesh->getElement(elem).getBarycenter();
            
            for (Index i=0;i<NLocDoF;i++)
            {
                //Recover the position of the Dof
                Point uvw = _GLE_bndy.GetPointCoordinate(i);
            
                //Recover the jacobian
                std::array<std::array<Real, DIM>, DIM > F;
                std::array<Real, DIM> S;
                
                _mesh->getElement(elem).getJacobian(uvw,F);
                
                if (DIM == 2) _Jacobians_bndy[l][el][i]
                    = sqrt(F[0][0]*F[0][0]+F[1][0]*F[1][0]);
                else
                {
                    //To check !!!!!!!!!!!!!!!!!!!!! TODO
                    //assert(false);
                    //Should be wrong, use the transpose of F instead
                    ArrayAlgebra::CrossProduct(F[0],F[1],S);
                    _Jacobians_bndy[l][el][i] = fabs(ArrayAlgebra::Norm(S));
                }
            }
        }
        
    }*/
    
    //Construct the coloring of the mesh and construct the array
    std::vector<Index> colors;
    getColoring(_mesh,colors, _NColor);

    _elemByColor.resize(_NColor);
    
    for (int el=0;el<_NElem;el++)
    {
        Index elem = _elements[el];
        Index color = colors[elem];
        if (color>=0)
            _elemByColor[color].push_back(el);
    }
    
    //Init the array glob2elem
    //_initGlob2Elem();
    
}

/*
template<Index Dim> void FEMesh<Dim>::_initGlob2Elem()
{
    _glob2elem.resize(_NDoF);
    
    for (int el=0;el<_NElem;++el)
    {
        for (int i=0;i<_loc2glob[el].size();++i)
        {
            _glob2elem[_loc2glob[el][i]].push_back(el);
        }
    }
    
    //Sort and remove duplicate
    for (int el=0;el<_NElem;++el)
    {
        std::sort(_glob2elem[el].begin(),_glob2elem[el].end());
        auto last = std::unique(_glob2elem[el].begin(),_glob2elem[el].end());
        _glob2elem[el].erase(last, _glob2elem[el].end());
    }
}
*/
} 




