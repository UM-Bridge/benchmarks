#pragma once

//#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <thread>


// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

namespace _VTK
{


// ---------------------------------------------------------------------------------//
// Writing solution in .vtk file format.
// ofs is a stream created from the file.
inline void Write(std::ostream& ofs, TetriXP2Space & aFESpace, const LAL::Vector & aSolutionVect)
{
    //Computation of the system dimension
    Index nDof = aFESpace.getP2NumDoFs();
    Index nDofGlob = aFESpace.getNumDoFs();
    Index SysDim = LAL::getDimension(aSolutionVect)/nDofGlob;

    // Interpolated mesh points coordinates.
    RealVector xyz;

    //
    std::vector<Index> loc2glob;
    
    // 
    Index NumTet = aFESpace.getNumElements();


    // Writing header.
    ofs << "# vtk DataFile Version 2.0\n";
    ofs << "OndoMathX TetriX VTK IO\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
    
    ofs << "POINTS " << nDof << " double \n";

    ofs.precision(8);
    ofs << std::scientific;

    for (Index iDoF = 0; iDoF < nDof; iDoF++)
    {
        // Extracting coordinates of the DoF of the template finite element space.
        aFESpace.getP2DoFCoordinate(iDoF, xyz);
    
        ofs << double(xyz[0]) << " " << double(xyz[1]) << " " << double(xyz[2]) << " \n";
    }
    
    // Writing elements.
    ofs << "CELLS " << NumTet << " " << 11 * NumTet << "\n";
    
    for (Index elem = 0; elem < NumTet; elem++)
    {

        aFESpace.P2Loc2Glob(elem,loc2glob);
        
        // Storing to stream.
        ofs << 10 << " "
        << loc2glob[0] << " " << loc2glob[1] << " " << loc2glob[2] << " " << loc2glob[3] << " "
        << loc2glob[4] << " " << loc2glob[5] << " " << loc2glob[6] << " " << loc2glob[7] << " " << loc2glob[8] << " " << loc2glob[9] <<  "\n";
    }

    // Writing element types.
    ofs << "CELL_TYPES " << NumTet << "\n";

    for (Index i = 0; i < NumTet; i++)    ofs << 24 << "\n";
    
    const Real * aSolutionData = LAL::getData(aSolutionVect);
    
    // if (SysDim == 1)
    // {
        
        
        ofs << "POINT_DATA " << nDof << std::endl;
        ofs << "SCALARS SOLUTION double 1" << std::endl;
        ofs << "LOOKUP_TABLE default" << std::endl;
        
        for (Index iDoF = 0; iDoF < nDof; iDoF++)
        {
            Index iGlob = aFESpace.P2Glob2Glob(iDoF);

            double v =  double(aSolutionData[iGlob]) ;

            if (fabs(v)<1e-40) v = 0.0;

            ofs << v << std::endl; //Add the conversion to float to avoid a bug in paraview
        }

    // }
    // else
    // {
        
    //     Index N = LAL::getDimension(aSolutionVect) / SysDim;
        
    //     ofs << "POINT_DATA " << N << std::endl;
    //     ofs << "VECTORS SOLUTION double" << std::endl;
        
    //     for (Index iDoF = 0; iDoF < N; iDoF++)
    //     {
    //         switch (SysDim){
                    
    //             case 2: {
                    
    //                 double v1 = double(aSolutionData[SysDim*iDoF]);
    //                 double v2 = double(aSolutionData[SysDim*iDoF+1]);

    //                 if (fabs(v1)<1e-40) v1 = 0.0;
    //                 if (fabs(v2)<1e-40) v2 = 0.0;

    //                 ofs << v1 << " " << v2 << " " << 0.0 << std::endl;
                    
    //                 break; }
    //             case 3: {

    //                 double v1 = double(aSolutionData[SysDim*iDoF]);
    //                 double v2 = double(aSolutionData[SysDim*iDoF+1]);
    //                 double v3 = double(aSolutionData[SysDim*iDoF+2]);

    //                 if (fabs(v1)<1e-40) v1 = 0.0;
    //                 if (fabs(v2)<1e-40) v2 = 0.0;
    //                 if (fabs(v3)<1e-40) v3 = 0.0;
                    
    //                 ofs << v1 << " " << v2 << " " << v3 << std::endl;
                    
    //                 break; }
                    
    //             default : {assert(false);};
    //         }
    //     }
    // }
    
}
// ---------------------------------------------------------------------------------//

} // Namespace _VTK

inline void WriteVTK(const std::string& FilePath, TetriXP2Space & aFESpace, LAL::Vector & aSolution)
{
    std::ofstream ofs(FilePath);
    _VTK::Write(ofs, aFESpace, aSolution);
    ofs.close();
}
// ---------------------------------------------------------------------------------//

 
} 
 