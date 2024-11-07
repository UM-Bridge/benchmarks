

namespace OndoMathX {
    
 
        inline void read_MSH(const char *fileName, std::vector<std::shared_ptr<Point> > &nodes, 
                                std::vector<std::shared_ptr<Element> > &elements)
        {
            nodes.clear();
            elements.clear();

            std::string s_tmp;
            int i_tmp,j_tmp,t_tmp;
            Index Indexmp;
            Index ui_tmp;
            Real d_tmp1,d_tmp2,d_tmp3;
                       
            //Open the file
            std::ifstream file(fileName);
            if (!file) assert(false);
            
            //Read the format
            file >> s_tmp; if (s_tmp != "$MeshFormat") {assert(false);}
            file >> d_tmp1;if (d_tmp1 >= 2.3) {assert(false);}
            file >> i_tmp; if (i_tmp != 0) {assert(false);}
            file >> i_tmp; if (i_tmp != 8) {assert(false);}
            file >> s_tmp; if (s_tmp != "$EndMeshFormat") {assert(false);}
            
            //Read the physical names (if needed)
            file >> s_tmp;
            if (s_tmp == "$PhysicalNames")
            {
                while (s_tmp != "$EndPhysicalNames" && !file.eof()) file >> s_tmp;
                file >> s_tmp;
            }
            
            
            //Read the vertices
            if (s_tmp != "$Nodes") {assert(false);}
            file >> Indexmp; nodes.resize(Indexmp);
            for (Index i=0;i<Indexmp;++i)
            {
                file >> ui_tmp; //index (base 1) !!
                file >> d_tmp1;
                file >> d_tmp2;
                file >> d_tmp3;
                nodes[i] = std::make_shared<Point>(d_tmp1,d_tmp2,d_tmp3,ui_tmp-1);
            }
            file >> s_tmp; if (s_tmp != "$EndNodes") {assert(false);}
            
            //Read the elements
            
            GeoElement type = GeoNone;
            int itype;
            int partition;
            Index numV=0;      //gives the number of Coord to read
            
            //temporary vector of pointers on vertices (uses for construction)
            std::vector<std::shared_ptr<Point> > nodes_tmp;
            
            
            file >> s_tmp; if (s_tmp != "$Elements") {assert(false);}
            file >> Indexmp; //read the number of elements
            elements.reserve(Indexmp);
            
            for (Index i=0;i<Indexmp;i++)
            {
                file >> j_tmp; //index of the element
                file >> itype;  //type of the element
                file >> t_tmp; //number of tagsvariable
                
                if (t_tmp > 0)
                {
                    file >> partition; //partition in which the element belong
                    for (int j=1;j<t_tmp;j++) file >> i_tmp;
                }
                else
                {
                    partition = 0;
                }
               
                bool construct_element = true;
                
                //compute the number of vertices to read
                switch(itype)
                {
                    case GeoPoint       : numV = 1;  type = GeoPoint; construct_element = false;
                        break;
                    case GeoSegment2    : numV = 2;  type = GeoSegment2;
                        break;
                    case GeoTriangle3   : numV = 3;  type = GeoTriangle3;
                        break;
                    case GeoQuadrangle4 : numV = 4;  type = GeoQuadrangle4;
                        break;
                    case GeoTetrahedron4: numV = 4;  type = GeoTetrahedron4;
                        break;
                    case GeoHexahedron8 : numV = 8;  type = GeoHexahedron8;
                        break;
                    case GeoSegment3    : numV = 3;  type = GeoSegment3;
                        break;
                    case GeoQuadrangle9 : numV = 9;  type = GeoQuadrangle9;
                        break;
                    case GeoHexahedron27: numV = 27; type = GeoHexahedron27;
                        break;
                        
                    default : construct_element = false;
                }
                
                nodes_tmp.resize(numV);
                
                //read the vertices
                for (Index j=0;j<numV;j++)
                {
                    file >> ui_tmp;  //read the number of the Coord (index base = 1 !!)
                    nodes_tmp[j]=nodes[ui_tmp-1];
                }
                
                //create the element if the element is not a Coord;
                if (construct_element)
                    elements.push_back(Factory(type,nodes_tmp,partition));
                
            }
            
            file >> s_tmp;
            if (s_tmp != "$EndElements") {assert(false);}
            
            //close the file
            file.close();
        }
        
        inline void print_MSH(const Mesh & mesh, std::ostream & flux)
        {
            //intermediate numbering
            std::map<Index, Index> num;
            for (Index i=0;i<mesh.getNumNodes();++i)
                num[mesh.getNode(i).getIndex()] = i;
            
            //Write the format
            flux << "$MeshFormat" << std::endl;
            flux << 2.2 << " " << 0 << " " << 8 << std::endl;
            flux << "$EndMeshFormat" <<std:: endl;
            
            //Write the nodes
            flux << "$Nodes" << std::endl;
            flux << mesh.getNumNodes() << std::endl;
            for (Index i=0;i<mesh.getNumNodes();++i)
                flux << i+1 << " " << mesh.getNode(i).x() << " " << mesh.getNode(i).y() << " " << mesh.getNode(i).z() << std::endl;
            
            flux << "$EndNodes" << std::endl;
            
            
            flux << "$Elements" << std::endl;
            flux << mesh.getNumElements() << std::endl;
            
            for (Index el=0;el<mesh.getNumElements();++el)
            {
                const Element & elem = mesh.getElement(el);
                
                flux << el << " " << elem.getType() << " " << 2 << " " << elem.getLabel() << " " << elem.getLabel() << " ";
                
                for (Index j=0;j<elem.getNumNodes(); ++j)
                {
                    flux << num[elem.getNode(j).getIndex()]+1;
                    
                    if ( j == (elem.getNumNodes()-1) ) flux << std::endl;
                    else flux << " ";
                    
                }
            }
            
            flux << "$EndElements" << std::endl;
        }
 
    
} // namespace OndoMathX




