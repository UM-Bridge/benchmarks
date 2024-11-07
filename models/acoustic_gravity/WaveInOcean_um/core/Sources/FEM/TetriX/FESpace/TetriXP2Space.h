#pragma once


#include <set>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include <cassert>


namespace OndoMathX {

 
class TetriXP2Space  
{
public:
    static constexpr Index Dim = 3;
    static constexpr Index Order = 2;
    
private:

    static constexpr Index _NLocDoF = 15;   // Local dimension of the enriched P2 space
    static constexpr Index _P2_NLocDoF = 10; // Local dimension of the P2 space

    //Basic informations for the generic numbering routine/////////////////////////////////////////////
    //! Reference on the mesh
    const Mesh & _mesh;

 
    
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
   
    
    // Quadrature weight, barycenter, coordinates, Jacobian for each element
    std::vector<RealVector> _Barycenter;
    std::vector<RealVector> _DofCoordinates;
    std::vector<RealMatrix3x3> _GradDef_Origin;

    std::vector<RealMatrix3x3> _ComatGradDef_Origin;
    std::vector<Real> _Jacobian_Origin;

    std::vector<RealMatrix3x3> _SymGradDef_Origin;

    std::vector<bool> _Affine;

    // Loc to glob
    std::vector<std::vector<Index>> _loc2glob;

    //Correspondance with P2 finite elements
    Index _P2_NDoF;
    std::vector<std::vector<Index>> _P2_loc2glob;
    std::vector<Index> _P2_glob2glob;


    //! Everything related to boundaries (with map by label)
    //std::map<Index, std::shared_ptr<Mesh>> _mesh_bndy;
    //std::map<Index, std::vector<std::vector<Index>>> _loc2glob_bndy;
    //std::map<Index, std::vector<Point>> _Barycenter_bndy;
    //std::map<Index, std::vector<std::vector<Real>>> _Jacobians_bndy;
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
public:

    TetriXP2Space(const Mesh &aMesh);
    

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
    
    // ---------------------------------------------------------------------------------//
    
    
    
    // ---------------------------------------------------------------------------------//
    /*! \brief Extracting the total number of DoF.
     \return The number of  DoF.
     */
    Index getNumDoFs() const
    {
        return _NDoF;
    }

    Index getP2NumDoFs() const
    {
        return _P2_NDoF;
    }

 	Index getNumLocDoFs() const
	{
		return _NLocDoF;
	}   

    Index getNumLocP2Dofs() const
    {
        return _P2_NLocDoF;
    }

    const RealVector& getLocDofCoordinates(Index iloc) const
    {
        return TetriX::Tet_15::Point[iloc];
    }
    
    bool isAffine(const Index iEltGlob) const
    {
        return _Affine[iEltGlob];  
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
    Index Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        return _loc2glob[iEltGlob][iLoc];
    }

    Index P2Loc2Glob(Index iEltGlob, Index iLoc) const
    {
        return _P2_loc2glob[iEltGlob][iLoc];
    }
    
    //Index Loc2GlobBndy(Index label,Index iEltBndy, Index iLoc) const
    //{
    //    return _loc2glob_bndy.at(label)[iEltBndy][iLoc];
    //}
    
    void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        loc2glob = _loc2glob[iEltGlob];
    }

    void P2Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
    {
        loc2glob = _P2_loc2glob[iEltGlob];
    }

    Index P2Glob2Glob(Index iEltGlob) const
    {
        return _P2_glob2glob[iEltGlob];
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
     \param xyz are the coordinates of the degree of freedom associated to iLoc which will be filled.
     */
    void getDoFCoordinate(Index iGlob, RealVector & xyz) const
    {
        xyz = _DofCoordinates[iGlob];
    }

    void getDoFCoordinate(Index iElem, Index iloc, RealVector & xyz) const
    {
        xyz = _DofCoordinates[_loc2glob[iElem][iloc]];
    }
    
    void getP2DoFCoordinate(Index iGlob, RealVector & xyz) const
    {
         xyz = _DofCoordinates[_P2_glob2glob[iGlob]];
    }
 
    /*! \brief Gives the barycenter of the element
     \param iEltGlob is a global index of an element.
     \param xyzCenter are the coordinates of the barycenter of the element
     */
    void getBarycenter(const Index iEltGlob, RealVector & xyzCenter) const
    {
        xyzCenter = _Barycenter[iEltGlob];
    }
    
    /*! \brief Gives the barycentric coordinates of a local coords.
     \param uvw_local are the local coords.
     \param barycentric_coords is the updated barycentric coords array. 
     */
    void getBarycentricCoordinates(const RealVector& uvw_local, std::array<Real, 4>& barycentric_coords) const
    {
        const auto u = uvw_local[0];
        const auto v = uvw_local[1];
        const auto w = uvw_local[2];

        barycentric_coords[0] = 1 - u - v - w;
        barycentric_coords[1] = u;
        barycentric_coords[2] = v;
        barycentric_coords[3] = w;
    }
    
   // Real GetJacobian(Index iEltGlob, Index iLoc) const
   // {
   //     return _Jacobians[iEltGlob][iLoc];
   // }
    
    //Real GetJacobianBndy(Index label, Index iEltGlob, Index iLoc) const
    //{
    //    return _Jacobians_bndy.at(label)[iEltGlob][iLoc];
    //}
    
    void getGradDef(Index iEltGlob, Index iLoc, RealMatrix3x3 & GradDef) const
    {            
       GradDef = _GradDef_Origin[iEltGlob];   
    }     
    // ---------------------------------------------------------------------------------//
    

    void AssembleMass(LAL::DiagonalMatrix &mass, Real alpha) const;
    void MltStiffness(const LAL::Vector &U, LAL::Vector &V, Real alpha) const;



};

} 

#include "TetriXP2Space.hxx"

