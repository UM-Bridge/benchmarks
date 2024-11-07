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
namespace OndoMathX
{

    namespace _VIZIR
    {

        inline std::thread _WriterThread;

        // ---------------------------------------------------------------------------------//
        template <class FESpace>
        inline void WriteSolution(std::string FileName, const FESpace &aFESpace, const Mesh &aMesh, const LAL::Vector &solutionCopy, Index N)
        {
            // Computation of the system dimension
            Index DimSol = LAL::getDimension(solutionCopy) / aFESpace.getNumDoFs();

            char buffer[256];
            std::snprintf(buffer, 256, "%.5u", (unsigned int)N);

            std::string filePath = FileName + "." + buffer + ".sol";
            std::ofstream ofs(filePath);

            const auto dimension = FESpace::Dim;

            // Writing header
            ofs << "MeshVersionFormatted 2" << std::endl;
            ofs << "Dimension" << std::endl
                << dimension << std::endl;
            ofs << "HOSolAt";

            // Get geometric element type. Assumes that we only do I/O for the highest dimension.
            std::string geo_element_tag{""};
            const std::vector<Index> &element_id_list = aMesh.getElementList(dimension);
            const auto &first_element = aMesh.getElement(element_id_list[0]);
            const auto element_type = first_element.getType(); // Assumption: Only one type of elements in the mesh

            switch (element_type)
            {
            case GeoQuadrangle4:
                geo_element_tag = "QuadrilateralsQ1";
                break;
            case GeoQuadrangle9:
                geo_element_tag = "QuadrilateralsQ2";
                break;
            case GeoHexahedron8:
                geo_element_tag = "HexahedraQ1";
                break;
            case GeoHexahedron27:
                geo_element_tag = "HexahedraQ2";
                break;
            case GeoPoint:
            case GeoSegment2:
            case GeoTriangle3:
            case GeoTriangle6:
            case GeoTetrahedron4:
            case GeoNone:
            {
                assert(false);
            }
            break;
            }

            // Writing local finite element
            ofs << geo_element_tag << "NodesPositions\n";

            // Get the reference element definition
            const GaussLobattoElement &gle = aFESpace.getFE();
            assert(gle.getNI() == gle.getNJ());
            if (dimension == 3)
                assert(gle.getNI() == gle.getNK());

            ofs << gle.getNumPoints() << std::endl;
            for (size_t iloc = 0; iloc < gle.getNumPoints(); iloc++)
            {
                const RealVector &uvw = gle.getPointCoordinate(iloc);
                for (Index idim = 0; idim < dimension; ++idim)
                    ofs << uvw[idim] << " ";
                ofs << '\n';
            }

            // Writing finite element space.
            ofs << "HOSolAt" << geo_element_tag << '\n';
            ofs << aFESpace.getNumElements() << std::endl;
            ofs << "1 " << DimSol << std::endl;
            ofs << gle.getNI() - 1 << " " << gle.getNumPoints() << std::endl;

            const Real *solutionData = LAL::getData(solutionCopy);

            for (size_t iElt = 0; iElt < aFESpace.getNumElements(); iElt++)
            {
                for (Index iloc = 0; iloc < gle.getNumPoints(); iloc++)
                {
                    for (Index k = 0; k < DimSol; k++)
                    {
                        ofs << solutionData[DimSol * aFESpace.Loc2Glob(iElt, iloc) + k] << " ";
                    }
                }

                ofs << '\n';
            }

            ofs << "End";
            ofs.close();
        }

        // ---------------------------------------------------------------------------------//

    } // _VIZIR

    inline void ParallelWriteVIZIRTerminate()
    {
        if (_VIZIR::_WriterThread.joinable())
            _VIZIR::_WriterThread.join();
    }

    template <class FESpace>
    inline void ParallelWriteSolutionVIZIR(const std::string &FileName, const FESpace &aFESpace, const Mesh &aMesh, const LAL::Vector &aSolutionVect, const Index N = 0)
    {

        if (_VIZIR::_WriterThread.joinable())
            _VIZIR::_WriterThread.join();

        // Writing solution.
        _VIZIR::_WriterThread = std::thread([FileName, &aFESpace, &aMesh, &aSolutionVect, N]
                                            {
        LAL::Vector copy = aSolutionVect;
        _VIZIR::WriteSolution(FileName, aFESpace, aMesh, copy, N); });
    }

    template <class FESpace>
    inline void WriteSolutionVIZIR(const std::string &FileName, const FESpace &aFESpace, const Mesh &aMesh, const LAL::Vector &aSolutionVect, const Index N = 0)
    {
        _VIZIR::WriteSolution(FileName, aFESpace, aMesh, aSolutionVect, N);
    }

} // OndoMathX

// -----------------------------------------------------------------------------------------//
