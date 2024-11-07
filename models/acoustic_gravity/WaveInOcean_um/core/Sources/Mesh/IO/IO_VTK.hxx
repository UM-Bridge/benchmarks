
namespace OndoMathX {
    
        inline void print_VTK(const Mesh & mesh, std::ostream & flux, std::string name = "OndoMathX VTK IO")
        {
            
            //Write the supported elements
            Index NumSegments2=0;     //
            Index NumTriangles3=0;    // Various counter
            Index NumQuadrangles4=0;  //
            Index NumTetrahedra4=0;   //
            Index NumHexahedra8=0;    //
            
            std::ostringstream fluxSegments2(std::ostringstream::out);    //
            std::ostringstream fluxTriangles3(std::ostringstream::out);   // Temporary fluxes
            std::ostringstream fluxQuadrangles4(std::ostringstream::out); //
            std::ostringstream fluxTetrahedra4(std::ostringstream::out);  //
            std::ostringstream fluxHexahedra8(std::ostringstream::out);   //
            
            std::vector<Index> labelSegments2;                                //
            std::vector<Index> labelTriangles3;                               //
            std::vector<Index> labelQuadrangles4;                             // Store the labels
            std::vector<Index> labelTetrahedra4;                              //
            std::vector<Index> labelHexahedra8;                               //
            
            
            //intermediate numbering 
            std::map<Index, Index> num;
            for (Index i=0;i<mesh.getNumNodes();++i)
                num[mesh.getNode(i).getIndex()] = i;
            
            //Construct the temporary fluxes
            for (Index el=0;el<mesh.getNumElements();++el)
            {
                const Element & elem = mesh.getElement(el);
                
                switch (elem.getType())
                {
                    case GeoSegment3 :
                    case GeoSegment2 :
                        ++NumSegments2;
                        labelSegments2.push_back(elem.getLabel());
                        fluxSegments2 <<2<<" "<< num[elem.getNode(0).getIndex()] << " " << num[elem.getNode(1).getIndex()] << " \n";
                        break;
                    case GeoTriangle6:
                    case GeoTriangle3:
                        ++NumTriangles3;
                        labelTriangles3.push_back(elem.getLabel());
                        fluxTriangles3 <<3<<" "<< num[elem.getNode(0).getIndex()] << " " << num[elem.getNode(1).getIndex()] << " "
                        << num[elem.getNode(2).getIndex()] << "\n";
                        break;
                    case GeoQuadrangle9:
                    case GeoQuadrangle4:
                        ++NumQuadrangles4;
                        labelQuadrangles4.push_back(elem.getLabel());
                        fluxQuadrangles4 <<4<<" "<< num[elem.getNode(0).getIndex()] << " " << num[elem.getNode(1).getIndex()] << " "
                        << num[elem.getNode(2).getIndex()] << " " << num[elem.getNode(3).getIndex()] << "\n" ;
                        break;
                    case GeoTetrahedron4:
                        ++NumTetrahedra4;
                        labelTetrahedra4.push_back(elem.getLabel());
                        fluxTetrahedra4 <<4<<" "<< num[elem.getNode(0).getIndex()] << " " << num[elem.getNode(1).getIndex()] << " "
                        << num[elem.getNode(2).getIndex()] << " " << num[elem.getNode(3).getIndex()] << "\n";
                        break;
                    case GeoHexahedron27:
                    case GeoHexahedron8:
                        ++NumHexahedra8;
                        labelHexahedra8.push_back(elem.getLabel());
                        fluxHexahedra8 <<8<<" "<< num[elem.getNode(0).getIndex()] << " " << num[elem.getNode(1).getIndex()] << " "
                        << num[elem.getNode(2).getIndex()] << " " << num[elem.getNode(3).getIndex()] << " "
                        << num[elem.getNode(4).getIndex()] << " " << num[elem.getNode(5).getIndex()] << " "
                        << num[elem.getNode(6).getIndex()] << " " << num[elem.getNode(7).getIndex()] << " \n";
                        break;
                    case GeoPoint:
                    case GeoNone:         {assert(false);} break;
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            //Write the mesh in the flux
            
            // HEADER
            flux<<"# vtk DataFile Version 2.0\n";
            flux<< name <<"\n";
            flux<<"ASCII\n";
            flux<<"DATASET UNSTRUCTURED_GRID \n";
            
            //POINTS
            flux <<"POINTS "<< mesh.getNumNodes() <<" double \n";
            for (Index i=0;i<mesh.getNumNodes();i++)
                flux<< mesh.getNode(i).x() << " " << mesh.getNode(i).y() << " " << mesh.getNode(i).z() <<" \n";
            
            //CELLS
            flux<<"CELLS "<<NumSegments2+NumTriangles3+NumQuadrangles4+NumTetrahedra4+NumHexahedra8
            <<" "<<3*NumSegments2+4*NumTriangles3+5*NumQuadrangles4+5*NumTetrahedra4+9*NumHexahedra8<<"\n";
            
            flux<<fluxSegments2.str();
            flux<<fluxTriangles3.str();
            flux<<fluxQuadrangles4.str();
            flux<<fluxTetrahedra4.str();
            flux<<fluxHexahedra8.str();
            
            //CELLS TYPES
            flux<<"CELL_TYPES "<<NumSegments2+NumTriangles3+NumQuadrangles4+NumTetrahedra4+NumHexahedra8<<"\n";
            
            for (Index i=0;i<NumSegments2;i++)      flux<<3<<"\n";
            for (Index i=0;i<NumTriangles3;i++)  	flux<<5<<"\n";
            for (Index i=0;i<NumQuadrangles4;i++)	flux<<9<<"\n";
            for (Index i=0;i<NumTetrahedra4;i++) 	flux<<10<<"\n";
            for (Index i=0;i<NumHexahedra8;i++)  	flux<<12<<"\n";   
        }
} // namespace OndoMathX




