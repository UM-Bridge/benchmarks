
namespace OndoMathX {


    void LinearIntersection(Real L, Index N1, Index N2, std::vector<std::vector<Index>> &NI, std::vector<std::vector<Real>> &xI, Real Tolerance = 1e-14)
    {
        assert(N1<N2);

        NI.resize(N1);
        xI.resize(N1);

        Index n2=0;

        Real dx2 = L/N2;

        for (Index n1=0;n1<N1;++n1)
        {
            Real xl = n1*L/N1;
            Real xr = (n1+1)*L/N1;
            Real x_current;

            xI[n1].push_back(xl);
            

            //Add the first element
            {
                NI[n1].push_back(n2);
                n2++;
                xI[n1].push_back(n2*dx2);
            }

            //Add elements inside
            while (((n2+1)*dx2)<(xr-Tolerance))
            {
                NI[n1].push_back(n2);
                n2++;
                xI[n1].push_back(n2*dx2);
            }

            //Add the last element  
            {
                NI[n1].push_back(n2);
                xI[n1].push_back(xr);

                if (fabs((n2+1)*dx2-xr) < Tolerance) n2++;
            }
            

        }
    }




    
        // Based on the assumption that the first mesh is a sub-mesh of the second one
        // the function return a vector indicating in which element of the second mesh
        // each element of the first one is included in
        inline std::vector<Index>
        IntersectionMatrix(const Mesh &Mesh1, const Mesh &Mesh2)
        {
            std::vector<Index> intersection_matrix(Mesh1.getNumElements());
            
            for (Index el1 = 0; el1 < Mesh1.getNumElements(); el1++)
            {
                for (Index el2 = 0; el2 < Mesh2.getNumElements(); el2++)
                {
                                    
                    if (Mesh1.getElement(el1).getType() == GeoSegment2 && Mesh2.getElement(el2).getType() == GeoSegment2 )
                    {
                        RealVector barycenter = Mesh1.getElement(el1).getBarycenter();
                     
                        bool is_inside = Mesh2.getElement(el2).isInside(barycenter);
                        
                        if (is_inside)
                        {
                            intersection_matrix[el1]=el2;
                            break;
                        }
                    }
                    else
                    {
                        assert(false);
                    }
                }
            }
            
            return intersection_matrix;
        }



    
        inline void ElementIntersection(const Element &element1, const Element &element2,
                                 std::vector<Point> &intersection_points, std::vector<bool> &search_vector)
        {
            assert(element1.getDim()<=2 && element2.getDim()==2);
            
            
            //Firt step: edge intersection
            search_vector.assign(element2.getNumEdges(),false);
            intersection_points.clear();
            
            
            for (Index j=0;j<element2.getNumEdges();++j)
            {
                Edge e2 = element2.getEdge(j);
                
                for (Index i=0;i<element1.getNumEdges();++i)
                {
                    
                    Edge e1 = element1.getEdge(i);
                    
                    std::array<double,2> b;
                    std::array<double,2> r;
                    std::array< std::array<double, 2>, 2 > A;
                    
                    A[0][0]= e1.getVertex(0).x()-e1.getVertex(1).x();
                    A[0][1]=-e2.getVertex(0).x()+e2.getVertex(1).x();
                    A[1][0]= e1.getVertex(0).y()-e1.getVertex(1).y();
                    A[1][1]=-e2.getVertex(0).y()+e2.getVertex(1).y();
                    
                    
                    
                    if (fabs(ArrayAlgebra::Det(A))>REF_COORD_TOL*REF_COORD_TOL)
                    {
                        b[0]=-e1.getVertex(1).x()+e2.getVertex(1).x();
                        b[1]=-e1.getVertex(1).y()+e2.getVertex(1).y();
                        
                        RealMatrix2x2 invA;
                        ArrayAlgebra::Inv(A,invA);
                        ArrayAlgebra::MatMlt(invA,b,r);
                        
                        if (r[0]>=0 && r[0]<=1 && r[1]>=0 && r[1]<=1)
                        {
                            Point intersection_point(r[0]*e1.getVertex(0).x()+(1-r[0])*e1.getVertex(1).x(),
                                                     r[0]*e1.getVertex(0).y()+(1-r[0])*e1.getVertex(1).y(),0.0);
                            
                            intersection_points.push_back(intersection_point);
                            search_vector[j]=true;
                        }
                    }
                }
            }
            
            //Second step: check wether or not vertices are inside elements
            for (Index j=0;j<element2.getNumVertices();++j)
            {
                if (element1.isInside(element2.getVertex(j).coord()))
                {
                    intersection_points.push_back(element2.getVertex(j));
                    
                    //Add all the search direction corresponding to edges that include the vertex j
                    for (Index k=0;k<element2.getNumEdges();++k)
                    {
                        Edge e2 = element2.getEdge(k);
                        
                        if (e2.getVertex(0).getIndex() == element2.getVertex(j).getIndex()
                            || e2.getVertex(1).getIndex() == element2.getVertex(j).getIndex())
                            
                            search_vector[k]=true;
                        
                    }
                }
            }
            
            for (Index i=0;i<element1.getNumVertices();++i)
            {
                if (element2.isInside(element1.getVertex(i).coord()))
                {
                    intersection_points.push_back(element1.getVertex(i));
                }
            }
            
            //Third step: sort the intersection points counter-clockwise
            if (intersection_points.size() < 2) intersection_points.clear();
            else
            {
                struct lessCoordinatesXYZ cmp_XYZ;
                struct sameCoordinatesXYZ same_XYZ;
                struct lessCoordinatesAngle2D cmp_angle;
                
                std::sort(intersection_points.begin(),intersection_points.end(),cmp_XYZ);
                auto it = std::unique (intersection_points.begin(),intersection_points.end(),same_XYZ);
                intersection_points.erase(it, intersection_points.end());
                
                cmp_angle.Barycenter = getBarycenter(intersection_points);
                std::sort(intersection_points.begin(),intersection_points.end(),cmp_angle);
               
            }
            
        }
        
        //Compute the intersected mesh between mesh1 and mesh2. THe first given mesh can be a mesh of lines.
        inline Mesh meshIntersection(const Mesh &Mesh1, std::vector<std::vector<Index> > &IntersectedElem1,
                              const Mesh &Mesh2, std::vector<std::vector<Index> > &IntersectedElem2)
        {
            assert(Mesh1.getDim()==2);
            assert(Mesh2.getDim()==2);
            assert(Mesh1.getDimElement()<=2);
            assert(Mesh2.getDimElement()==2);
            
            Index dimElement1 = Mesh1.getDimElement();
            

            std::vector<std::shared_ptr<Point> > new_nodes;
            std::vector<std::shared_ptr<Element> > new_elements;
            
            Index global_index = 0;
            
            
            //vectors helping through the search of neigbours
            std::vector<bool> search_vector_tmp;
            std::vector<bool> search_vector;
            std::vector<Index> next_candidates_1;
            
            //intersections point list
            std::vector<Point> intersection_points;
            
            //Get a couple of indices corresponding to two elements that intersects
            Index begin_el_1=0;
            Index begin_el_2=0;
            
            for (Index el2:Mesh2.getElementList(2))
            {
                for (Index el1:Mesh1.getElementList(dimElement1))
                {
                    ElementIntersection(Mesh1.getElement(el1),Mesh2.getElement(el2),intersection_points,search_vector);
                    
                    auto it = std::find(search_vector.begin(),search_vector.end(),true);
                    if (it != search_vector.end())
                    {
                        begin_el_1=el1;
                        begin_el_2=el2;
                        break;
                    }
                }
                
                auto it = std::find(search_vector.begin(),search_vector.end(),true);
                if (it != search_vector.end()) break;
            }
            
            
            //list of elements that must be treated
            std::list<Index> treatment_list_1;
            std::list<Index> treatment_list_2;
            treatment_list_2.push_back(begin_el_2);
            std::list<Index> starting_treatment_list_1;
            starting_treatment_list_1.push_back(begin_el_1);
            
            //boolean flag that indicates if an element shoudl be treated or not
            std::vector<bool> treatment_flag_2(Mesh2.getNumElements());
            std::vector<bool> treatment_flag_1(Mesh1.getNumElements());
            treatment_flag_2[begin_el_2] = true;
            
            while (treatment_list_2.size() != 0)
            {
                Index current_element_2 = treatment_list_2.front();
                treatment_list_2.pop_front();
                
                treatment_list_1.clear();
                treatment_list_1.push_back(starting_treatment_list_1.front());
                starting_treatment_list_1.pop_front();
                
                treatment_flag_1.assign(Mesh1.getNumElements(),false);
                treatment_flag_1[treatment_list_1.front()] = true;

                next_candidates_1.assign(Mesh2.getElement(current_element_2).getNumEdges(),0);
                search_vector.assign(next_candidates_1.size(),false);
                
                while (treatment_list_1.size() != 0)
                {
                    Index current_element_1 = treatment_list_1.front();
                    treatment_list_1.pop_front();
                    
                    //Compute all the interection points and the search_vector
                    intersection_points.clear();
                    search_vector_tmp.clear();
                    
                    ElementIntersection(Mesh1.getElement(current_element_1),
                                        Mesh2.getElement(current_element_2),
                                        intersection_points,search_vector_tmp);
                    
                    if (intersection_points.size() != 0)
                    {
                        //Add the vertices
                        Index first_coord = new_nodes.size();
                        
                        for (Index i=0;i<intersection_points.size();++i)
                        {
                            new_nodes.push_back(std::make_shared<Point>(intersection_points[i][0],
                                                                        intersection_points[i][1],
                                                                        intersection_points[i][2],
                                                                        global_index));
                            
                            global_index++;
                        }
                        
                        
                        /*
                         Routine that compute the intersected mesh from the intersection_points
                         */
                        if ((Mesh1.getElement(current_element_1).getType() == GeoTriangle3
                             || Mesh1.getElement(current_element_1).getType() == GeoQuadrangle4) && intersection_points.size() > 2)
                        {
                            if (intersection_points.size() == 4)
                            {
                                std::vector< std::shared_ptr<Point> > v(4);
                                v[0] = new_nodes[first_coord];
                                v[1] = new_nodes[first_coord+1];
                                v[2] = new_nodes[first_coord+2];
                                v[3] = new_nodes[first_coord+3];
                                
                                new_elements.push_back(std::make_shared<Quadrangle4>(v));
                            }
                            else for (Index i=2;i<intersection_points.size();++i)
                            {
                                std::vector< std::shared_ptr<Point> > v(3);
                                v[0] = new_nodes[first_coord];
                                v[1] = new_nodes[first_coord+i-1];
                                v[2] = new_nodes[first_coord+i];
                                
                                new_elements.push_back(std::make_shared<Triangle3>(v));
                            }
                        }
                        else if(Mesh1.getElement(current_element_1).getType() == GeoSegment2 && intersection_points.size() == 2)
                        {
                            std::vector< std::shared_ptr<Point> > v(2);
                            v[0] = new_nodes[first_coord];
                            v[1] = new_nodes[first_coord+1];
                            
                            new_elements.push_back(std::make_shared<Segment2>(v));
                        }
                        else {}
                        
                        std::vector<Index> neighbour_elements_2 = Mesh2.getNeighbourElements(current_element_2);
                        std::vector<Index> neighbour_elements_1 = Mesh1.getNeighbourElements(current_element_1);
                        
                        
                        if (Mesh1.getElement(current_element_1).getType() == GeoSegment2)
                        {
                            //Check if the intersection points are all on a common edge; put the treatment flag to true for that neighbour element
                            for(Index el : neighbour_elements_2)
                            {
                                if ( Mesh2.getElement(el).isInside(intersection_points[0].coord())
                                    && Mesh2.getElement(el).isInside(intersection_points[1].coord()) )
                                {
                                    treatment_flag_2[el] = true;
                                    break;
                                }
                            }
                        }
                        
                        //Put the element in the treatment list for the mesh 1 (and the treatment flag to 1)
                        for(Index el : neighbour_elements_1)
                        {
                            if (treatment_flag_1[el]==false)
                            {
                                treatment_list_1.push_back(el);
                                treatment_flag_1[el]=true;
                            }
                        }
                                  
                        for(Index interface = 0; interface < next_candidates_1.size() ; interface++ )
                        {
                            if (search_vector_tmp[interface]==true)
                            {
                                //Add the neigbouring element to the corresponding interface to a temporary list
                                next_candidates_1[interface]=current_element_1;
                                search_vector[interface]=true;
                            }
                        }
                    }
                    
                    std::vector<Index> neighbour_elements_2_by_interface;
                    std::vector<bool> has_neighbour_2;
                    Mesh2.getNeighbourElementsByInterface(current_element_2,has_neighbour_2,neighbour_elements_2_by_interface);
                    
                    //Put the element in the treatment list for the mesh 2 if the search_direction flag and the treatment are ok
                    for(Index interface = 0; interface < search_vector.size() ; interface++ )
                    {
                        if (has_neighbour_2[interface])
                        {
                            //Add the neigbouring element to the corresponding interface to a temporary list
                            Index next_element_2 = neighbour_elements_2_by_interface[interface];
                            
                            if (treatment_flag_2[next_element_2] == false && search_vector[interface]==true)
                            {
                                treatment_list_2.push_back(next_element_2);
                                
                                starting_treatment_list_1.push_back(next_candidates_1[interface]);
                                
                                treatment_flag_2[next_element_2] = true;
                                
                            }
                            
                        }
                    }
                }
            }
            
            return Mesh(new_nodes,new_elements);
        }

    
} // namespace OndoMathX




