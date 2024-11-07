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
inline std::thread _WriterThread;


// ---------------------------------------------------------------------------------//
// Writing finite element space in .vtk file format.
//ofs is a stream created from the file.
template<class FESpace>
inline void WriteFESpace(std::ostream& ofs, FESpace & aFESpace, bool continuous = true)
{
    assert(FESpace::Dim==2 || FESpace::Dim==3);
    
    // Get the reference element definition
    const GaussLobattoElement & gle = aFESpace.getFE();
    
    // Interpolated mesh points coordinates.
    RealVector xyz;
    
    // Counter if elements in the interpolated mesh.
    Index NumEdges = 0;
    Index NumQuadrangles = 0;
    Index NumHexahedra = 0;
    
    // Flux for connections in the interpolated mesh.
    std::ostringstream FluxEdges(std::ostringstream::out);
    std::ostringstream FluxQuadrangles(std::ostringstream::out);
    std::ostringstream FluxHexahedra(std::ostringstream::out);
    
    // Writing header.
    ofs << "# vtk DataFile Version 2.0\n";
    ofs << "OndoMathX LobattoX VTK IO\n";
    ofs << "ASCII\n";
    ofs << "DATASET UNSTRUCTURED_GRID\n";
    
    // Writing points.
    Index nDof;

    if (continuous) 
    {
        nDof = aFESpace.getNumDoFs();
    }
    else 
    {
        nDof = aFESpace.getNumDoFsDisc();
    }

    ofs << "POINTS " << nDof << " double \n";

    for (Index iDoF = 0; iDoF < nDof; iDoF++)
    {
        // Extracting coordinates of the DoF of the template finite element space.
        if (continuous) 
        {
            aFESpace.getDoFCoordinate(iDoF, xyz);
        }
        else
        {
            aFESpace.getDoFCoordinateDisc(iDoF, xyz);
        }
        
        // Writing DoF coordinates.
        switch (FESpace::Dim)
        {
            case 2:{ ofs << double(xyz[0]) << " " << double(xyz[1]) << " " << 0. << " \n"; break; }
            case 3:{ ofs << double(xyz[0]) << " " << double(xyz[1]) << " " << double(xyz[2]) << " \n"; break; }
        }
    }
    
    // Subdivide each High-order element for Q1 vizualisation
    Index PX = gle.getNI()-1;;
    Index PY = gle.getNJ()-1;;
    Index PZ = gle.getNK()-1;
    
    if (FESpace::Dim==2) PZ=1;
    
    Index NumSubdivisions=PX*PY*PZ;

    std::vector<Index> loc2glob(gle.getNumPoints());
    
    if (FESpace::Dim == 2)
    {
        for (Index elem = 0; elem < aFESpace.getNumElements(); elem++)
        {
            if (continuous) 
            {
                aFESpace.Loc2Glob(elem,loc2glob);
            }
            else
            {
                aFESpace.Loc2GlobDisc(elem,loc2glob);
            }
            
            for (Index sub = 0; sub < NumSubdivisions; sub++)
            {
                Index offsetY = sub/PX;
                Index offsetX = sub - offsetY*PX;
                
                //Extracting element numbering
                std::vector<Index> ElementNum(4);
            
                ElementNum[0] = loc2glob[offsetX+offsetY*(PX+1)];
                ElementNum[1] = loc2glob[(offsetX+1)+offsetY*(PX+1)];
                ElementNum[2] = loc2glob[(offsetX+1)+(offsetY+1)*(PX+1)];
                ElementNum[3] = loc2glob[offsetX+(offsetY+1)*(PX+1)];
                
                // Updating number of quadrangles.
                NumQuadrangles++;
                
                // Storing to stream.
                FluxQuadrangles << 4 << " "
                << ElementNum[0] << " " << ElementNum[1] << " "
                << ElementNum[2] << " " << ElementNum[3] << "\n";
            }
        }
    }
    else
    {
        for (Index elem = 0; elem < aFESpace.getNumElements(); elem++)
        {
            if (continuous) 
            {
                aFESpace.Loc2Glob(elem,loc2glob);
            }
            else
            {
                aFESpace.Loc2GlobDisc(elem,loc2glob);
            }
            
            for (Index sub = 0; sub < NumSubdivisions; sub++)
            {
            
                Index offsetZ = sub/(PX * PY);
                Index offsetY = (sub - offsetZ*(PX * PY))/PX;
                Index offsetX = sub - offsetZ*(PX * PY) - offsetY*PX;
            
                //Extracting element numbering
                std::vector<Index> ElementNum(8);
            
                ElementNum[0] = loc2glob[offsetX+offsetY*(PX+1)+offsetZ*(PX+1)*(PY+1)];
                ElementNum[1] = loc2glob[(offsetX+1)+offsetY*(PX+1)+offsetZ*(PX+1)*(PY+1)];
                ElementNum[2] = loc2glob[(offsetX+1)+(offsetY+1)*(PX+1)+offsetZ*(PX+1)*(PY+1)];
                ElementNum[3] = loc2glob[offsetX+(offsetY+1)*(PX+1)+offsetZ*(PX+1)*(PY+1)];
                
                ElementNum[4] = loc2glob[offsetX+offsetY*(PX+1)+(offsetZ+1)*(PX+1)*(PY+1)];
                ElementNum[5] = loc2glob[(offsetX+1)+offsetY*(PX+1)+(offsetZ+1)*(PX+1)*(PY+1)];
                ElementNum[6] = loc2glob[(offsetX+1)+(offsetY+1)*(PX+1)+(offsetZ+1)*(PX+1)*(PY+1)];
                ElementNum[7] = loc2glob[offsetX+(offsetY+1)*(PX+1)+(offsetZ+1)*(PX+1)*(PY+1)];
                
                // Updating number of Hexahedra
                NumHexahedra++;
            
                // Storing to stream.
                FluxHexahedra << 8 << " "
                << ElementNum[0] << " " << ElementNum[1] << " " << ElementNum[2] << " " << ElementNum[3] << " "
                << ElementNum[4] << " " << ElementNum[5] << " " << ElementNum[6] << " " << ElementNum[7] << "\n";
            
            }
        }
    }
    
    
    // Writing elements.
    ofs << "CELLS " << NumEdges  + NumQuadrangles  + NumHexahedra
    << " " << 3 * NumEdges  + 5 * NumQuadrangles  + 9 * NumHexahedra << "\n";
    
    ofs << FluxEdges.str();
    ofs << FluxQuadrangles.str();
    ofs << FluxHexahedra.str();
    
    // Writing element types.
    ofs << "CELL_TYPES " << NumEdges  + NumQuadrangles + NumHexahedra << "\n";
    
    for (Index i = 0; i < NumEdges; i++)    ofs << 3 << "\n";
    for (Index i = 0; i < NumQuadrangles; i++)    ofs << 9 << "\n";
    for (Index i = 0; i < NumHexahedra; i++)    ofs << 12 << "\n";
}


// ---------------------------------------------------------------------------------//
// Writing solution in .vtk file format.
// ofs is a stream created from the file.
template<class FESpace>
inline void WriteSolution(std::ostream& ofs, FESpace & aFESpace, const LAL::Vector & aSolutionVect, bool continuous = true)
{
    Index nDof;

    if (continuous) 
    {
        nDof = aFESpace.getNumDoFs();
    }
    else 
    {
        nDof = aFESpace.getNumDoFsDisc();
    }

    //Computation of the system dimension
    Index SysDim = LAL::getDimension(aSolutionVect)/nDof;
    
    // Writing finite element space.
    WriteFESpace(ofs, aFESpace, continuous);
    
    ofs.precision(8);
    ofs << std::scientific;
    
    const Real* aSolutionData = LAL::getData(aSolutionVect);
    
    if (SysDim == 1){
        
        Index N =LAL::getDimension(aSolutionVect);
        
        ofs << "POINT_DATA " << N << std::endl;
        ofs << "SCALARS SOLUTION double 1" << std::endl;
        ofs << "LOOKUP_TABLE default" << std::endl;
        
        for (Index iDoF = 0; iDoF < N; iDoF++)
        {
            double v =  double(aSolutionData[iDoF]) ;

            if (fabs(v)<1e-40) v = 0.0;

            ofs << v << std::endl; //Add the conversion to float to avoid a bug in paraview
        }
    }
    else{
        
        Index N = LAL::getDimension(aSolutionVect) / SysDim;
        
        ofs << "POINT_DATA " << N << std::endl;
        ofs << "VECTORS SOLUTION double" << std::endl;
        
        for (Index iDoF = 0; iDoF < N; iDoF++)
        {
            switch (SysDim){
                    
                case 2: {
                    
                    double v1 = double(aSolutionData[SysDim*iDoF]);
                    double v2 = double(aSolutionData[SysDim*iDoF+1]);

                    if (fabs(v1)<1e-40) v1 = 0.0;
                    if (fabs(v2)<1e-40) v2 = 0.0;

                    ofs << v1 << " " << v2 << " " << 0.0 << std::endl;
                    
                    break; }
                case 3: {

                    double v1 = double(aSolutionData[SysDim*iDoF]);
                    double v2 = double(aSolutionData[SysDim*iDoF+1]);
                    double v3 = double(aSolutionData[SysDim*iDoF+2]);

                    if (fabs(v1)<1e-40) v1 = 0.0;
                    if (fabs(v2)<1e-40) v2 = 0.0;
                    if (fabs(v3)<1e-40) v3 = 0.0;
                    
                    ofs << v1 << " " << v2 << " " << v3 << std::endl;
                    
                    break; }
                    
                default : {assert(false);};
            }
        }
    }
    
}
// ---------------------------------------------------------------------------------//



// ---------------------------------------------------------------------------------//
// Writing solution in .vtk file format.
// ofs is a stream created from the file.
inline void WriteSolution(std::ostream& ofs, Mesh & aMesh, const LAL::Vector & aSolutionVect)
{
    Index nDof;

    nDof = aMesh.getNumElements();

    //Computation of the system dimension
    Index SysDim = LAL::getDimension(aSolutionVect)/nDof;
    
    print_VTK(aMesh,ofs);

    ofs.precision(8);
    ofs << std::scientific;
    
    const Real* aSolutionData = LAL::getData(aSolutionVect);
    
    if (SysDim == 1)
    {
       
        Index N =LAL::getDimension(aSolutionVect);
        
        ofs << "CELL_DATA  " << N << std::endl;
        ofs << "SCALARS SOLUTION double 1" << std::endl;
        ofs << "LOOKUP_TABLE default" << std::endl;
        
        for (Index iDoF = 0; iDoF < N; iDoF++)
        {
            double v =  double(aSolutionData[iDoF]) ;

            if (fabs(v)<1e-40) v = 0.0;

            ofs << v << std::endl; //Add the conversion to float to avoid a bug in paraview
        } 
    }
    else{
        
        Index N = LAL::getDimension(aSolutionVect) / SysDim;
        
        ofs << "CELL_DATA  " << N << std::endl;
        ofs << "VECTORS SOLUTION double" << std::endl;
        
        for (Index iDoF = 0; iDoF < N; iDoF++)
        {
            switch (SysDim){
                    
                case 2: {
                    
                    double v1 = double(aSolutionData[SysDim*iDoF]);
                    double v2 = double(aSolutionData[SysDim*iDoF+1]);

                    if (fabs(v1)<1e-40) v1 = 0.0;
                    if (fabs(v2)<1e-40) v2 = 0.0;

                    ofs << v1 << " " << v2 << " " << 0.0 << std::endl;
                    
                    break; }
                case 3: {

                    double v1 = double(aSolutionData[SysDim*iDoF]);
                    double v2 = double(aSolutionData[SysDim*iDoF+1]);
                    double v3 = double(aSolutionData[SysDim*iDoF+2]);

                    if (fabs(v1)<1e-40) v1 = 0.0;
                    if (fabs(v2)<1e-40) v2 = 0.0;
                    if (fabs(v3)<1e-40) v3 = 0.0;
                    
                    ofs << v1 << " " << v2 << " " << v3 << std::endl;
                    
                    break; }
                    
                default : {assert(false);};
            }
        }
    }
    
}
// ---------------------------------------------------------------------------------//





} //_VTK
// ---------------------------------------------------------------------------------//

// ---------------------------------------------------------------------------------//
inline void ParallelWriteVTKTerminate()
{
    if (_VTK::_WriterThread.joinable())
    _VTK::_WriterThread.join();
    
}

template<class FESpace>
inline void ParallelWriteVTK(const std::string& FilePath, FESpace & aFESpace, LAL::Vector & aSolutionVect, bool continuous = true)
{
    
    
    if (_VTK::_WriterThread.joinable())
        _VTK::_WriterThread.join();

    // Writing solution.
    _VTK::_WriterThread = std::thread([FilePath,&aFESpace,&aSolutionVect,continuous]{
        std::ofstream ofs(FilePath);
        LAL::Vector copy = aSolutionVect;
        _VTK::WriteSolution(ofs, aFESpace, copy, continuous);
        ofs.close();
    });
    
}

template<class FESpace>
inline void WriteVTK(const std::string& FilePath, FESpace & aFESpace,LAL::Vector & aSolutionVect, bool continuous = true)
{
    std::ofstream ofs(FilePath);
    _VTK::WriteSolution(ofs, aFESpace, aSolutionVect, continuous);
    ofs.close();
}


inline void WriteVTK(const std::string& FilePath, Mesh & aMesh,LAL::Vector & aSolution)
{
    std::ofstream ofs(FilePath);
    _VTK::WriteSolution(ofs, aMesh, aSolution);
    ofs.close();
}
// ---------------------------------------------------------------------------------//

 
}

// -----------------------------------------------------------------------------------------//
