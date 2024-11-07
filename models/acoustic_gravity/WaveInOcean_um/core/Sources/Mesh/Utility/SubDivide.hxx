
namespace OndoMathX {

inline Mesh Mesh::subDivide() const
{
    Mesh copied_mesh(*this);

    Index Nnodes = _nodes.size();

    //Consistent numbering of the nodes 
    for (Index i=0;i<Nnodes;++i)
    {
        copied_mesh._nodes[i]->setIndex(i);
    }

    //boolean table indicating if the node of the mesh to subdivide is used
    std::vector<bool> used_nodes(Nnodes);

    //Shared pointer to a duplicated node
    std::vector<std::shared_ptr<Point>> nodes_duplicate;  



 
	std::vector<std::shared_ptr<Element>> new_elements;         //New list of elements
    
	std::vector<std::shared_ptr<Point>> mid_edges(_numEdges);   //contain pointer on added vertices on edges
 
	std::vector<std::shared_ptr<Point>> mid_faces(_numFaces);   //contain pointer on added vertices on faces
   
    std::vector<std::shared_ptr<Point>> mid_volumes;            //contain pointer on added vertices on the volumes

    std::vector<std::shared_ptr<Point>> new_Q2_nodes;           //new nodes for Q2 mesh description

    std::vector<std::shared_ptr<Point>> vertices;               //vertices of the original mesh

    std::vector<std::shared_ptr<Point>> new_nodes;

	//Creation of the new list
	for (Index el = 0; el < _elements.size(); ++el)
	{

		auto pe = _elements[el];

		//Recover informations
		Index label = pe->getLabel();

        std:vector<std::shared_ptr<Point>> v(pe->getNumVertices());

		//for (Index i=0;i<pe->getNumVertices();++i) v[i]=pe->getVertexPtr(i);

		switch(pe->getType())
		{
			case GeoSegment2:
            {
			 
				auto mid_e = mid_edges[_neighbourEdges[el][0]];

				//Compute the coordinates of the new vertex on the edge
				*mid_e = Point((v[0]->x()+v[1]->x())*0.5,
							   (v[0]->y()+v[1]->y())*0.5,
							   (v[0]->z()+v[1]->z())*0.5); 

				//new_elements.push_back(std::make_shared<Segment2>({v[0],mid_e},label));
				//new_elements.push_back(std::make_shared<Segment2>({mid_e,v[1]},label));
            }
			break;

			case GeoQuadrangle4:
			{

                std::vector<std::shared_ptr<Point>> mid_e(4);

				for (Index i=0;i<4;++i) 
                {

                    mid_e[i] = mid_edges[_neighbourEdges[el][i]];

					auto ve1 = pe->getEdge(i).getVertex(0);
					auto ve2 = pe->getEdge(i).getVertex(1);

					*mid_e[i] = Point((ve1.x()+ve2.x())*0.5,
                                      (ve1.y()+ve2.y())*0.5,
                                      (ve1.z()+ve2.z())*0.5);
				 
				}

				auto mid_f = mid_faces[_neighbourFaces[el][0]];

				//Compute the coordinates of the new vertex on the face
				*mid_f = Point((v[0]->x()+v[1]->x()+v[2]->x()+v[3]->x())*0.25,
							   (v[0]->y()+v[1]->y()+v[2]->y()+v[3]->y())*0.25,
							   (v[0]->z()+v[1]->z()+v[2]->z()+v[3]->z())*0.25);

				//Construct the new quadrangles
				//new_elements.push_back(std::make_shared<Quadrangle4>({v[0],mid_e[0],mid_f,mid_e[3]},label));
				//new_elements.push_back(std::make_shared<Quadrangle4>({mid_e[0],v[1],mid_e[1],mid_f},label));
				//new_elements.push_back(std::make_shared<Quadrangle4>({mid_f,mid_e[1],v[2],mid_e[2]},label));
				//new_elements.push_back(std::make_shared<Quadrangle4>({mid_e[3],mid_f,mid_e[2],v[3]},label));
            }
			break;

			case GeoHexahedron8:
            {
                std::vector<std::shared_ptr<Point>> mid_e(12);
                std::vector<std::shared_ptr<Point>> mid_f(6);

				//Compute the coordinates of the new vertex on the face
				mid_volumes.push_back(std::make_shared<Point>(
                            (v[0]->x()+v[1]->x()+v[2]->x()+v[3]->x()+v[4]->x()+v[5]->x()+v[6]->x()+v[7]->x())*0.125,
							(v[0]->y()+v[1]->y()+v[2]->y()+v[3]->y()+v[4]->y()+v[5]->y()+v[6]->y()+v[7]->y())*0.125,
							(v[0]->z()+v[1]->z()+v[2]->z()+v[3]->z()+v[4]->z()+v[5]->z()+v[6]->z()+v[7]->z())*0.125));

                auto mid  = mid_volumes.back();


				//Get the middle vertices on edges
                for (Index i=0;i<12;++i) 
                {

                    mid_e[i] = mid_edges[_neighbourEdges[el][i]];

					auto ve1 = pe->getEdge(i).getVertex(0);
					auto ve2 = pe->getEdge(i).getVertex(1);

					*mid_e[i]  = Point((ve1.x()+ve2.x())*0.5,
                                       (ve1.y()+ve2.y())*0.5,
                                       (ve1.z()+ve2.z())*0.5);
				 
				}

                //Get the middle vertices on faces
			    for (Index i=0;i<6;++i) 
                {

                    mid_f[i] = mid_faces[_neighbourFaces[el][i]];

					auto vf1 = pe->getFace(i).getVertex(0);
					auto vf2 = pe->getFace(i).getVertex(1);
					auto vf3 = pe->getFace(i).getVertex(2);
					auto vf4 = pe->getFace(i).getVertex(3);

					*mid_f[i]  = Point((vf1.x()+vf2.x()+vf3.x()+vf4.x())*0.25,
                                       (vf1.y()+vf2.y()+vf3.y()+vf4.y())*0.25,
                                       (vf1.z()+vf2.z()+vf3.z()+vf4.z())*0.25);
				 
				}
 

				//Construct the new hexahedra
                /*
				new_elements.push_back(std::make_shared<Hexahedron8>({v[0],mid_e[0],mid_f[0],mid_e[3],mid_e[4],mid_f[1],mid,mid_f[2]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({v[1],mid_e[1],mid_f[0],mid_e[0],mid_e[5],mid_f[3],mid,mid_f[1]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({v[2],mid_e[2],mid_f[0],mid_e[1],mid_e[6],mid_f[4],mid,mid_f[3]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({v[3],mid_e[2],mid_f[0],mid_e[3],mid_e[7],mid_f[4],mid,mid_f[2]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({mid_e[4],mid_f[1],mid,mid_f[2],v[4],mid_e[8],mid_f[5],mid_e[11]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({mid_e[5],mid_f[3],mid,mid_f[1],v[5],mid_e[9],mid_f[5],mid_e[8]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({mid,mid_f[3],mid_e[6],mid_f[4],mid_f[5],mid_e[9],v[6],mid_e[10]},label));
				new_elements.push_back(std::make_shared<Hexahedron8>({mid_f[2],mid,mid_f[4],mid_e[7],mid_e[11],mid_f[5],mid_e[10],v[7]},label));*/
            }
			break;

			default: assert(false);
		}

 
	}
 
    //new_nodes = mid_edges;
    //new_nodes.insert(new_nodes.end(), mid_faces.begin(), mid_faces.end());
    //new_nodes.insert(new_nodes.end(), mid_volumes.begin(), mid_volumes.end());
   
    return Mesh(new_nodes,new_elements);
}
    
} // namespace OndoMathX




