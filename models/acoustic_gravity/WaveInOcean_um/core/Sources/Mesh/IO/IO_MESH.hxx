

namespace OndoMathX {
    
        
        inline void read_MESH(const char *fileName, std::vector<std::shared_ptr<Point> > &nodes,
                            std::vector<std::shared_ptr<Element> > &elements)
        {
            nodes.clear();
            elements.clear();

            std::string s_tmp;
            int i_tmp,j_tmp, dim;
            Index Indexmp;
            Real d_tmp1,d_tmp2,d_tmp3;
            
            //Open the file
            std::ifstream file(fileName);
            if (!file) assert(false);
            
            //Read the format
            file >> s_tmp; if (s_tmp != "MeshVersionFormatted") assert(false);
            file >> i_tmp; if (i_tmp != 0 && i_tmp != 1 && i_tmp != 2) assert(false);
            
            //Read the dimension
            file >> s_tmp; if (s_tmp != "Dimension") assert(false);
            file >> dim;
            
            //Read the vertices
            file >> s_tmp; if (s_tmp != "Vertices") assert(false);
            file >> Indexmp; nodes.resize(Indexmp);
            
            
            for (Index i=0;i<Indexmp;++i)
            {
                file >> d_tmp1;
                file >> d_tmp2;
                if (dim == 3) file >> d_tmp3; else d_tmp3 = 0.;
                file >> j_tmp;
                nodes[i] = std::make_shared<Point>(d_tmp1,d_tmp2,d_tmp3,i);
            }
            
            //Read the elements
            int elements_left = true;  //indicates if there is elements left to read
            GeoElement type = GeoNone;
            Index numV=0;                  //gives the number of Coord to read
            Index n;
            
            //temporary vector of pointers on vertices (used for construction)
            std::vector<std::shared_ptr<Point>> nodes_tmp;
            
            while(elements_left)
            {
                //Read the kind of element
                file >> s_tmp;
                
                if      (s_tmp == "Edges")           { type = GeoSegment2;    numV=2; }
                else if (s_tmp == "Triangles")       { type = GeoTriangle3;   numV=3;}
                else if (s_tmp == "Quadrangles" ||
                         s_tmp == "Quadrilaterals")  { type = GeoQuadrangle4; numV=4;}
                else if (s_tmp == "Tetrahedra")      { type = GeoTetrahedron4;numV=4;}
                else if (s_tmp == "Hexahedra")       { type = GeoHexahedron8; numV=8;}
                else if (s_tmp == "End" || s_tmp == "END") elements_left = false;
                else {assert(false);}
                
                nodes_tmp.resize(numV);
                
                //Read the elements
                if (elements_left)
                {
                    file >> n;
                    elements.reserve(elements.size()+n);
                    
                    for (Index i=0;i<n;++i)
                    {
                        for (Index j=0;j<numV;++j)
                        {
                            file >> Indexmp;   //read the number of the Coord (index base = 1 !!)
                            nodes_tmp[j]=nodes[Indexmp-1];
                        }
                        file >> i_tmp;          //read the element partition
                        
                        elements.push_back(Factory(type,nodes_tmp,i_tmp));
                    }
                }
            }
            
            //close the file
            file.close();
        }
        
        inline void print_MESH(const Mesh & mesh, std::ostream & flux)
        {
            //Write the begining of the file
            flux << "MeshVersionFormatted 2" << std::endl;
            flux << "Dimension" << std::endl;
            
            if (mesh.getDim() == 1) flux << 2 << std::endl;
            else flux << mesh.getDim() << std::endl;
            
            
            //Write the vertices
            flux << "Vertices" << std::endl;
            flux << mesh.getNumNodes() << std::endl;
            if (mesh.getDim() == 3)
                for (Index i=0;i<mesh.getNumNodes();++i)
                    flux << mesh.getNode(i).x() << " " << mesh.getNode(i).y() << " " << mesh.getNode(i).z() << " "<< 0 << std::endl;
            else for (Index i=0;i<mesh.getNumNodes();++i)
                flux << mesh.getNode(i).x() << " " << mesh.getNode(i).y() << " " << 0 << std::endl;
            
            //intermediate numbering (non optimal way of doing it, but guarrentes that the function is const)
            std::map<Index, Index> num;
            for (Index i=0;i<mesh.getNumNodes();++i)
                num[mesh.getNode(i).getIndex()] = i;
            
            //Write the supported elements
            Index Num_Edges=0;        //
            Index Num_Triangles=0;    // Various counter
            Index Num_Quadrangles=0;  //
            Index Num_Tetrahedra=0;  //
            Index Num_Hexahedra=0;  //
            Index Num_HexahedraQ2=0;  //
            
            std::ostringstream flux_Edges(std::ostringstream::out);       //
            std::ostringstream flux_Triangles(std::ostringstream::out);   // Temporary fluxes
            std::ostringstream flux_Quadrangles(std::ostringstream::out); //
            std::ostringstream flux_Tetrahedra(std::ostringstream::out); //
            std::ostringstream flux_Hexahedra(std::ostringstream::out); //
            std::ostringstream flux_HexahedraQ2(std::ostringstream::out); //
            
            //Construct the temporary fluxes
            for (Index el=0;el<mesh.getNumElements();++el)
            {
                const Element & elem = mesh.getElement(el);
                
                switch (elem.getType())
                {
                    case GeoSegment3:
                    case GeoSegment2:
                        ++Num_Edges;
                        flux_Edges << num[elem.getNode(0).getIndex()]+1 << " " << num[elem.getNode(1).getIndex()]+1 << " " << elem.getLabel() << "\n";
                        break;
                    case GeoTriangle6:
                    case GeoTriangle3:
                        ++Num_Triangles;
                        flux_Triangles << num[elem.getNode(0).getIndex()]+1 << " " << num[elem.getNode(1).getIndex()]+1 << " "
                        << num[elem.getNode(2).getIndex()]+1 << " " << elem.getLabel() << "\n";
                        break;
                    case GeoQuadrangle9:
                    case GeoQuadrangle4:
                        ++Num_Quadrangles;
                        flux_Quadrangles << num[elem.getNode(0).getIndex()]+1 << " " << num[elem.getNode(1).getIndex()]+1 << " "
                        << num[elem.getNode(2).getIndex()]+1 << " " << num[elem.getNode(3).getIndex()]+1 << " " << elem.getLabel() << "\n";
                        break;
                    case GeoTetrahedron4:
                        ++Num_Tetrahedra;
                        flux_Tetrahedra <<  num[elem.getNode(0).getIndex()]+1 << " " << num[elem.getNode(1).getIndex()]+1 << " "
                        << num[elem.getNode(2).getIndex()]+1 << " " << num[elem.getNode(3).getIndex()]+1 << " " << elem.getLabel() << "\n";
                        break;
                    case GeoHexahedron27:
                        ++Num_HexahedraQ2;
                        flux_HexahedraQ2 
                        << num[elem.getNode(0).getIndex()]+1 << " " << num[elem.getNode(1).getIndex()]+1 << " "
                        << num[elem.getNode(2).getIndex()]+1 << " " << num[elem.getNode(3).getIndex()]+1 << " "
                        << num[elem.getNode(4).getIndex()]+1 << " " << num[elem.getNode(5).getIndex()]+1 << " "
                        << num[elem.getNode(6).getIndex()]+1 << " " << num[elem.getNode(7).getIndex()]+1 << " " 
                        
                        << num[elem.getNode(8).getIndex()]+1 << " " << num[elem.getNode(11).getIndex()]+1 << " " << num[elem.getNode(13).getIndex()]+1 << " " << num[elem.getNode(9).getIndex()]+1 << " "
                        << num[elem.getNode(16).getIndex()]+1 << " " << num[elem.getNode(18).getIndex()]+1 << " " << num[elem.getNode(19).getIndex()]+1 << " " << num[elem.getNode(17).getIndex()]+1 << " "  
                        << num[elem.getNode(10).getIndex()]+1 << " " << num[elem.getNode(12).getIndex()]+1 << " " << num[elem.getNode(14).getIndex()]+1 << " " << num[elem.getNode(15).getIndex()]+1 << " "                
                        
                        << num[elem.getNode(20).getIndex()]+1 << " " << num[elem.getNode(25).getIndex()]+1 << " " << num[elem.getNode(21).getIndex()]+1 << " "
                        << num[elem.getNode(23).getIndex()]+1 << " " << num[elem.getNode(24).getIndex()]+1 << " " << num[elem.getNode(22).getIndex()]+1 << " "

                        << num[elem.getNode(26).getIndex()]+1 << " "
                        
                        << elem.getLabel() << "\n";
                        break;
                    case GeoHexahedron8:
                        ++Num_Hexahedra;
                        flux_Hexahedra << num[elem.getNode(0).getIndex()]+1 << " " << num[elem.getNode(1).getIndex()]+1 << " "
                        << num[elem.getNode(2).getIndex()]+1 << " " << num[elem.getNode(3).getIndex()]+1 << " "
                        << num[elem.getNode(4).getIndex()]+1 << " " << num[elem.getNode(5).getIndex()]+1 << " "
                        << num[elem.getNode(6).getIndex()]+1 << " " << num[elem.getNode(7).getIndex()]+1 << " " << elem.getLabel() << "\n";
                        break;
                    case GeoPoint:
                    case GeoNone:         {assert(false);} break;
                        
                }
            }
            
            //Write the temporary fluxes
            if (Num_Edges != 0) flux << "Edges\n"  << Num_Edges << "\n" << flux_Edges.str();
            if (Num_Triangles != 0) flux << "Triangles\n" << Num_Triangles << "\n" << flux_Triangles.str();
            if (Num_Quadrangles != 0) flux << "Quadrilaterals\n"<< Num_Quadrangles << "\n" << flux_Quadrangles.str();
            if (Num_Tetrahedra != 0) flux << "Tetrahedra\n" << Num_Tetrahedra << "\n" << flux_Tetrahedra.str();
            if (Num_Hexahedra != 0) flux << "Hexahedra\n" << Num_Hexahedra <<"\n" << flux_Hexahedra.str();
            if (Num_HexahedraQ2 != 0) flux << "HexahedraQ2\n" << Num_HexahedraQ2 << "\n" << flux_HexahedraQ2.str();
            
            flux << "End" << "\n";
        }
        
   
    
    
} // namespace OndoMathX




