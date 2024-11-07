#pragma once


#include <set>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include <cassert>


namespace OndoMathX {

template<
Index DIM
>
class FEMesh : public FESpaceT<FEMesh<DIM>,DIM>
{
    
protected:
    
    //Basic informations for the generic numbering routine/////////////////////////////////////////////
    //! Reference on the mesh
    const Mesh & _mesh;

    Index _Order; 

    GaussLobattoElement _GLE;
   
    
    //! List of geometric element of dimension DIM
    std::vector<Index> _elements;

    //!List of element per global dof (somewhat the reverse of the loc2glob aray)
    //std::vector< std::vector<Index> > _glob2elem;
    
////////////////////////////////////////////////////////////////////////////////////////////////////
    
    //Data for the finite element space interface///////////////////////////////////////////////////////
    /*! \brief Number of elements by color. */
    Index _NColor;
    std::vector< std::vector<Index> > _elemByColor;
    
    /*! \brief Total number of elements. */
    Index _NElem;
    
    /*! \brief Number of DoF. */
    Index _NDoF;
    Index _NLocDoF;
    
    // Quadrature weight, barycenter, coordinates, Jacobian for each element
    std::vector<RealVector> _Barycenter;
    std::vector<RealVector> _DofCoordinates;
    
    //! Everything related to boundaries (with map by label)
    //std::map<Index, std::shared_ptr<Mesh>> _mesh_bndy;
    //std::map<Index, std::vector<std::vector<Index>>> _loc2glob_bndy;
    //std::map<Index, std::vector<Point>> _Barycenter_bndy;
    //std::map<Index, std::vector<std::vector<Real>>> _Jacobians_bndy;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    FEMesh() {}
    
public:
    
    static const Index Dim = DIM;

    //Constructor !
    FEMesh(const Mesh &aMesh, Index order);
    
    
    const GaussLobattoElement & getFE() const
    {
        return _GLE;
    }
    
    /*
    GaussGaussLobattoElement & getFEBndy(Index Label el)
    {
        return _GLE_bndy;
    }*/
    
    //Print the dof (only for debug mode)
    //void print(std::ostream& out) const;
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Extracting number of color in the template geometry, used for parallelization.
     \return The number of colors.
     */
    Index getNumColors() const
    {
        return _NColor;
    }
    
    /*! \brief Extracting number of element associated to a given color.
     \return The number of element of a color group.
     */
    Index getNumElements(Index iColor) const
    {
        return _elemByColor[iColor].size();
    }
    
    
    /*! \brief Extracting number of global number of element.
     \return The number of element in the template.
     */
    Index getNumElements() const
    {
        return _NElem;
    }
    
    
    //Mesh & getMeshBndy(Index label)
    //{
    //    return *_mesh_bndy.at(label);
    //}
    
    
    // ---------------------------------------------------------------------------------//
    
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Extracting the total number of DoF.
     \return The number of  DoF.
     */
    Index getNumDoFs() const
    {
        return _NDoF;
    }

 	Index getNumLocDoFs() const
	{
		return _NLocDoF;
	}   
    
    bool isAffine(const Index iEltGlob) const
    {
        return false; //Could be optimised
    }
    
    /*! \brief fonction that maps the numbering of the element by color to a global numbering.
     \param c is the index of the color.
     \param e is the numbering of the element in the color e.
     \return the global numbering of the element.
     */
    Index EltColorNum2EltGlobalNum(Index c, Index e) const
    {
        return _elemByColor[c][e];
    }
    
    /*! \brief Extracting global DoF index from local DoF index.
     \param iEltGlob is the index of the element in the color group.
     \param iLoc is the local DoF index.
     \return The global index.
     */
    Index _Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        return FESpaceT<FEMesh<Dim>,Dim>::_loc2glob[iEltGlob][iLoc];
    }
    
    //Index Loc2GlobBndy(Index label,Index iEltBndy, Index iLoc) const
    //{
    //    return _loc2glob_bndy.at(label)[iEltBndy][iLoc];
    //}
    
    void _Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        loc2glob = FESpaceT<FEMesh<Dim>,Dim>::_loc2glob[iEltGlob];
    }
    
    // Get the label of the element (here 0)
    Index getLabel(const Index iEltGlob) const
    {
        return _mesh.getElement(_elements[iEltGlob]).getLabel();
    }

    
    //std::vector<Index> Loc2GlobBndy(Index label, Index iEltBndy, std::vector<Index> & loc2glob) const
    //{
    //    loc2glob = _loc2glob_bndy.at(label)[iEltBndy];
    //}
    
    
    
    
    /*! Return the elements that contains the input dof
     */
    //std::vector<Index> Glob2Elem(Index iEltGlob) const
    //{
    //    return _glob2elem[iEltGlob];
    //}
    
  

    
    // ---------------------------------------------------------------------------------//
    /*! \brief Gives a degree of freedom coordinate associated to a global DoF.
     \param iGlob is a global index of an DoF.
     \param xyz are the coordinates of the degree of freedo massociated to iLoc which will be filled.
     */
    void getDoFCoordinate(Index iGlob, RealVector & xyz) const
    {
        xyz = _DofCoordinates[iGlob];
    }
    
    void getDoFCoordinate(Index iElem, Index iloc, RealVector & xyz) const
    {
        xyz = _DofCoordinates[FESpaceT<FEMesh<Dim>,Dim>::_loc2glob[iElem][iloc]];
    }
    
 
    
    
    /*! \brief Gives the barycenter of the element
     \param iEltGlob is a global index of an element.
     \param xyzCenter are the coordinates of the barycenter of the element
     */
    void getBarycenter(const Index iEltGlob, RealVector & xyzCenter) const
    {
        xyzCenter = _Barycenter[iEltGlob];
    }
    
    
 
   // Real GetJacobian(Index iEltGlob, Index iLoc) const
   // {
   //     return _Jacobians[iEltGlob][iLoc];
   // }
    
    //Real GetJacobianBndy(Index label, Index iEltGlob, Index iLoc) const
    //{
    //    return _Jacobians_bndy.at(label)[iEltGlob][iLoc];
    //}
    
    void _getGradDef(Index iEltGlob, Index iLoc, std::array<std::array<Real, Dim>, Dim >& GradDef) const
    {            
        //Recover the position of the Dof on the reference element
        RealVector uvw = _GLE.getPointCoordinate(iLoc);

        _mesh.getElement(iEltGlob).getGradDef(uvw,GradDef);  

    }     
    // ---------------------------------------------------------------------------------//
    
    



    

    
    // ---------------------------------------------------------------------------------//
    /*! \brief Get the number of element of the interpolated mesh used for Q1 representation
     of high order finite element space.
     */

    Index GetNInterpolatedElement() const
    {
        if (Dim==2)
            return _NElem * _Order * _Order;
        else
            return _NElem * _Order * _Order * _Order;
    }


    
    
    /*! \brief Get the index of the DoFs of an element of the interpolated mesh used for P1 representation
     of high order finite element space.
     \param iInterpElt is the number of an element of the interpolated mesh.
     \param ElementInd are the vector gathering the index  the DoFs of the element.
     */
    void GetInterpolatedElementConnection(Index iInterpElt, std::vector<Index>& ElementInd) const
    {
        Index PX = _Order;

        std::vector<std::vector<Index>> & l2g = FESpaceT<FEMesh<Dim>,Dim>::_loc2glob;
        
        if (Dim==2)
        {
            //Recover the element index that is being used for interpolation
            Index elem = iInterpElt / (PX * PX);
            Index elem_loc = iInterpElt - elem * (PX * PX);
            
            Index offsetY = elem_loc/PX;
            Index offsetX = elem_loc - offsetY*PX;
            
            ElementInd[0] = l2g[elem][offsetX+offsetY*(PX+1)];
            ElementInd[1] = l2g[elem][(offsetX+1)+offsetY*(PX+1)];
            ElementInd[2] = l2g[elem][(offsetX+1)+(offsetY+1)*(PX+1)];
            ElementInd[3] = l2g[elem][offsetX+(offsetY+1)*(PX+1)];
            
            
        }
        else
        {
            //Recover the element index that is being used for interpolation
            Index elem = iInterpElt / (PX * PX * PX);
            Index elem_loc = iInterpElt - elem * (PX * PX * PX);
            
            Index offsetZ = elem_loc/(PX * PX);
            Index offsetY = (elem_loc - offsetZ*(PX * PX))/PX;
            Index offsetX = elem_loc - offsetZ*(PX * PX) - offsetY*PX;
            
            ElementInd[0] = l2g[elem][offsetX+offsetY*(PX+1)+offsetZ*(PX+1)*(PX+1)];
            ElementInd[1] = l2g[elem][(offsetX+1)+offsetY*(PX+1)+offsetZ*(PX+1)*(PX+1)];
            ElementInd[2] = l2g[elem][(offsetX+1)+(offsetY+1)*(PX+1)+offsetZ*(PX+1)*(PX+1)];
            ElementInd[3] = l2g[elem][offsetX+(offsetY+1)*(PX+1)+offsetZ*(PX+1)*(PX+1)];
            
            ElementInd[4] = l2g[elem][offsetX+offsetY*(PX+1)+(offsetZ+1)*(PX+1)*(PX+1)];
            ElementInd[5] = l2g[elem][(offsetX+1)+offsetY*(PX+1)+(offsetZ+1)*(PX+1)*(PX+1)];
            ElementInd[6] = l2g[elem][(offsetX+1)+(offsetY+1)*(PX+1)+(offsetZ+1)*(PX+1)*(PX+1)];
            ElementInd[7] = l2g[elem][offsetX+(offsetY+1)*(PX+1)+(offsetZ+1)*(PX+1)*(PX+1)];
            
        }
    }
    // ---------------------------------------------------------------------------------//
    
//private:
    
   // void _initGlob2Elem() ;
    
};

} 

#include "FEMesh.hxx"

