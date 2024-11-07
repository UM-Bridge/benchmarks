#pragma once


// -----------------------------------------------------------------------------------------//
// Misc. includes
#include <memory>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <thread>
// -----------------------------------------------------------------------------------------//



// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {
 
namespace _VIZIR {

 

// ---------------------------------------------------------------------------------//
inline void WriteSolution(std::string FileName, const TetriXP2Space & aFESpace, const Mesh& aMesh, const LAL::Vector & solutionCopy, Index N)
{
    //Computation of the system dimension
    Index DimSol = LAL::getDimension(solutionCopy)/aFESpace.getNumDoFs();
    
    char buffer[256];
    std::snprintf(buffer, 256, "%.5u",(unsigned int)N);

    std::string filePath = FileName + "." + buffer + ".sol";
    std::ofstream ofs(filePath);
    
    Index dimension = 3;
    
    // Writing header
    ofs << "MeshVersionFormatted 2" << std::endl;
    ofs << "Dimension" << std::endl << dimension << std::endl;
    ofs << "HOSolAt";

    // Get geometric element type. Assumes that we only do I/O for the highest dimension.
    std::string geo_element_tag = "TetrahedraP1";//{""};

/*
    const std::vector<Index>& element_id_list = aMesh.getElementList(dimension);
    const auto& first_element = aMesh.getElement(element_id_list[0]);
    const auto element_type = first_element.getType();
    
    switch (element_type)
    {
        case GeoSegment2:
            geo_element_tag = "EdgesP1";
            break;
        case GeoSegment3:
            geo_element_tag = "EdgesP2";
            break;
        case GeoTriangle3:
            geo_element_tag = "TrianglesP1";
            break;
        case GeoTriangle6:
            geo_element_tag = "TrianglesP2";
        case GeoQuadrangle4:
            geo_element_tag = "QuadrilateralsQ1";
            break;
        case GeoQuadrangle9:
            geo_element_tag = "QuadrilateralsQ2";
            break;
        case GeoTetrahedron4:
            geo_element_tag = "TetrahedraP1";
            break;
        case GeoHexahedron8:
            geo_element_tag = "HexahedraQ1";
            break;
        case GeoHexahedron27:
            geo_element_tag = "HexahedraQ2";
            break;
        case GeoPoint:
        case GeoNone:         
            {assert(false);} 
            break;
    }
    */

    // Writing local finite element
    ofs << geo_element_tag << "NodesPositions\n";
    ofs << aFESpace.getNumLocP2Dofs() << std::endl;

    std::array<Real, 4> barycentric_coords {0., 0., 0., 0.};
    for (size_t iloc  = 0; iloc <  aFESpace.getNumLocP2Dofs(); iloc++)
    {
        const RealVector & uvw = aFESpace.getLocDofCoordinates(iloc);
        aFESpace.getBarycentricCoordinates(uvw, barycentric_coords);
        for (auto barycentric_coord : barycentric_coords)
            ofs << barycentric_coord << " ";
        ofs << '\n';
    }
    
    // Writing finite element space.
    ofs << "HOSolAt" << geo_element_tag << '\n';
    ofs << aFESpace.getNumElements() << std::endl;
    ofs << "1 " << DimSol << std::endl;
    ofs << TetriXP2Space::Order << " " << aFESpace.getNumLocP2Dofs() << std::endl;
                
    const Real* solutionData = LAL::getData(solutionCopy);
    
    for (size_t iElt = 0; iElt < aFESpace.getNumElements(); iElt++)
    {
        for (Index iloc  = 0; iloc < aFESpace.getNumLocP2Dofs(); iloc++)
        {
            for (Index k  = 0; k < DimSol; k++)
            {
                Index iGlob_P2 = aFESpace.P2Loc2Glob(iElt,iloc);
                Index iGlob = aFESpace.P2Glob2Glob(iGlob_P2);
                ofs << solutionData[DimSol*iGlob+k] << " ";
            }

            
        }

        ofs << '\n';

    }
    
    ofs << "End";
    ofs.close();
    }
    
// ---------------------------------------------------------------------------------//




} // _VIZIR


inline void WriteSolutionVIZIR(const std::string& FileName, const TetriXP2Space & aFESpace, const Mesh& aMesh, const LAL::Vector & aSolutionVect, const Index N = 0)
{
    _VIZIR::WriteSolution(FileName, aFESpace, aMesh, aSolutionVect, N);
}


 

} // OndoMathX

// -----------------------------------------------------------------------------------------//

