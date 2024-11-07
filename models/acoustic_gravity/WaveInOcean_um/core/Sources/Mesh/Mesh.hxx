
namespace OndoMathX
{

    inline void Mesh::_copy(const Mesh &mesh)
    {
        std::vector<std::shared_ptr<Point>> nodes_tmp;

        _nodes.resize(mesh._nodes.size());

        std::map<Index, Index> num2Node;

        for (Index i = 0; i < mesh._nodes.size(); ++i)
        {
            _nodes[i] = std::make_shared<Point>(*mesh._nodes[i]);

            num2Node[mesh._nodes[i]->getIndex()] = i;
        }

        for (Index i = 0; i < mesh._elements.size(); ++i)
        {
            nodes_tmp.resize(mesh._elements[i]->getNumNodes());

            for (Index j = 0; j < mesh._elements[i]->getNumNodes(); ++j)
            {
                nodes_tmp[j] = _nodes[num2Node[mesh._elements[i]->getNode(j).getIndex()]];
            }

            _elements.push_back(Factory(mesh._elements[i]->getType(),
                                        nodes_tmp,
                                        mesh._elements[i]->getLabel()));
        }

        _dim = mesh._dim;
        _dimElement = mesh._dimElement;

        _neighbourEdges = mesh._neighbourEdges;
        _numEdges = mesh._numEdges;

        _neighbourFaces = mesh._neighbourFaces;
        _numFaces = mesh._numFaces;

        _neighbourVertices = mesh._neighbourVertices;
        _numVertices = mesh._numVertices;

        _neighbours = mesh._neighbours;

        _elementsByDim = mesh._elementsByDim;
        _elementsByDimAndLabel = mesh._elementsByDimAndLabel;
        _labels = mesh._labels;
        _labelsByDim = mesh._labelsByDim;

        _emptyVector = mesh._emptyVector;
    }

    // Mesh constructor//////////////////////////////////////////////////////////////////////
    inline Mesh::Mesh() : _dim(0), _dimElement(0), _numFaces(0), _numEdges(0) {}

    inline Mesh::~Mesh() { _clean(); }

    inline Mesh::Mesh(const char *fileName, MeshFormat format)
    {
        switch (format)
        {
        case FormatMesh:
            read_MESH(fileName, _nodes, _elements);
            break;
        case FormatMsh:
            read_MSH(fileName, _nodes, _elements);
            break;
        case FormatVtk:
        {
            assert(false);
        }
        break;
        }

        _initDim();
        _initElementSubList();
        _link();
    }

    inline Mesh::Mesh(const std::vector<std::shared_ptr<Point>> &nodes,
                      const std::vector<std::shared_ptr<Element>> &elements)
    {
        _nodes = nodes;
        _elements = elements;

        _initDim();
        _initElementSubList();
        _link();
    }

    inline Mesh::Mesh(const Mesh &mesh)
    {
        _copy(mesh);
    }

    inline Mesh &Mesh::operator=(const Mesh &mesh)
    {
        if (this == &mesh)
            return *this;

        _clean();

        _copy(mesh);

        return *this;
    }
    //////////////////////////////////////////////////////////////////////////////////////////

    inline void Mesh::_initDim()
    {
        _dim = 1;
        _dimElement = 1;

        for (Index i = 0; i < _nodes.size(); ++i)
        {
            if ((*_nodes[i])[1] != 0.0 & _dim < 2)
                _dim = 2;
            if ((*_nodes[i])[2] != 0.0 & _dim < 3)
                _dim = 3;
            if (_dim == 3)
                break;
        }

        for (Index i = 0; i < _elements.size(); ++i)
        {
            if (_dimElement < _elements[i]->getDim())
                _dimElement = _elements[i]->getDim();
            if (_dimElement == 3)
                break;
        }
    }

    inline void Mesh::_initElementSubList()
    {
        // I told ya ! It is just an empty vector
        _emptyVector.resize(0);

        _elementsByDim.clear();
        _elementsByDimAndLabel.clear();
        _labelsByDim.clear();

        _elementsByDim.resize(4);

        for (Index el = 0; el < _elements.size(); ++el)
        {
            _elementsByDim[_elements[el]->getDim()].push_back(el);
        }

        _elementsByDimAndLabel.resize(4);
        _labelsByDim.resize(4);
        _labels.resize(_elements.size());

        for (Index el = 0; el < _elements.size(); ++el)
        {
            _elementsByDimAndLabel[_elements[el]->getDim()][_elements[el]->getLabel()].push_back(el);

            _labelsByDim[_elements[el]->getDim()].push_back(_elements[el]->getLabel());

            _labels[el] = _elements[el]->getLabel();
        }

        // Store only one occurence of the label in _labelsByDim
        for (Index d = 0; d <= 3; ++d)
        {
            std::sort(_labelsByDim[d].begin(), _labelsByDim[d].end());
            auto last = std::unique(_labelsByDim[d].begin(), _labelsByDim[d].end());
            _labelsByDim[d].erase(last, _labelsByDim[d].end());
        }
    }

    inline const std::vector<Index> &Mesh::getElementList(Index dim) const
    {
        return _elementsByDim[dim];
    }

    inline const std::vector<Index> &Mesh::getElementList(Index dim, Index label) const
    {
        auto it = _elementsByDimAndLabel[dim].find(label);

        if (it == _elementsByDimAndLabel[dim].end())
            return _emptyVector;
        else
            return it->second;
    }

    inline const std::vector<Index> &Mesh::getLabels(Index dim) const
    {
        return _labelsByDim[dim];
    }

    inline const std::vector<Index> &Mesh::getLabels() const
    {
        return _labels;
    }

    inline void Mesh::setLabels(const std::vector<Index> &new_label)
    {
        assert(new_label.size() == _elements.size());

        for (Index el = 0; el < _elements.size(); ++el)
            _elements[el]->setLabel(new_label[el]);

        _initElementSubList();
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    inline void Mesh::_clean()
    {
        _elements.clear();
        _nodes.clear();
        _neighbourEdges.clear();
        _neighbourFaces.clear();
        _neighbourVertices.clear();
        _neighbours.clear();
        _elementsByDim.clear();
        _elementsByDimAndLabel.clear();
        _emptyVector.clear();
    }

    inline void Mesh::_clearLink()
    {
        _neighbourEdges.clear();
        _neighbourFaces.clear();
        _neighbourVertices.clear();
        _neighbours.clear();
        _numEdges = 0;
        _numFaces = 0;
    }

    //////////////////////////////////////////////////////////////////////////////////////////

    // the algorithm uses map, efficiency O(n)
    inline void Mesh::_link()
    {
        _clearLink();

        std::vector<Index> _tmp_numbering(_nodes.size());
        for (Index i = 0; i < _nodes.size(); ++i)
            _tmp_numbering[i] = _nodes[i]->getIndex();

        {
            std::vector<bool> verticesInserted(getNumNodes(), false);
            std::vector<Index> verticesNum(getNumNodes());

            _neighbourVertices.resize(_elements.size());

            // We use an intermediate consistent numbering to speed up the algorithm
            for (Index i = 0; i < _nodes.size(); ++i)
                _nodes[i]->setIndex(i);

            _numVertices = 0;

            for (Index el = 0; el < _elements.size(); ++el)
            {
                const Element &pe = *_elements[el];

                _neighbourVertices[el].resize(pe.getNumVertices());

                for (Index v = 0; v < pe.getNumVertices(); ++v)
                {
                    Index index = pe.getVertex(v).getIndex();

                    // Number the vertex if it nas not been numbered yet
                    if (verticesInserted[index] == false)
                    {
                        verticesNum[index] = _numVertices;
                        _numVertices++;
                        verticesInserted[index] = true;
                    }

                    // Add the vertices in the list of the element
                    _neighbourVertices[el][v] = verticesNum[index];
                }
            }

            // Recover the true numbering
            //(it must be done here because of the following computations of the orientation)
            for (Index i = 0; i < _nodes.size(); ++i)
                _nodes[i]->setIndex(_tmp_numbering[i]);

            // Construct the array Neighbours
            if (_dimElement == 1)
            {
                _neighbours.resize(_numVertices);

                for (auto el : getElementList(1))
                {

                    Index nvMax = _neighbourVertices[el].size();

                    for (Index nv = 0; nv < nvMax; ++nv)
                    {
                        _neighbours[_neighbourVertices[el][nv]].push_back(std::make_tuple(el, nv, 0));
                    }
                }
            }
        }

        {
            std::vector<std::map<Edge, Index>> mapEdges(getNumNodes());
            std::pair<std::map<Edge, Index>::iterator, bool> mapEdgesIt;

            _neighbourEdges.resize(_elements.size());

            // We use an antermediate numbering to speed up the algorithm
            for (Index i = 0; i < _nodes.size(); ++i)
                _nodes[i]->setIndex(i);

            for (Index el = 0; el < _elements.size(); ++el)
            {
                const Element &pe = *_elements[el];

                _neighbourEdges[el].resize(pe.getNumEdges());

                // Add edges if they do not already exist
                // 2D elements are referenced
                for (Index i = 0; i < pe.getNumEdges(); ++i)
                {
                    Edge edge = pe.getEdge(i);
                    Index num_Coord = edge.getSortedVertex(0).getIndex(); // sorted by their first Coord

                    // Insert the edge in the map
                    mapEdgesIt = mapEdges[num_Coord].insert(std::make_pair(edge, -1));
                    // The edge is numeroted if it has effectivly been inserted
                    if (mapEdgesIt.second == true)
                        mapEdgesIt.first->second = _numEdges++;
                    // Reference the neighbour relation
                    _neighbourEdges[el][i] = mapEdgesIt.first->second;
                }
            }

            // Recover the true numbering
            //(it must be done here because of the following computations of the orientation)
            for (Index i = 0; i < _nodes.size(); ++i)
                _nodes[i]->setIndex(_tmp_numbering[i]);

            // Construct the array Neighbours
            if (_dimElement == 2)
            {
                _neighbours.resize(_numEdges);

                for (auto el : getElementList(2))
                {
                    assert(_neighbourEdges[el].size() < std::numeric_limits<Index>::max());

                    Index neMax = (Index)(_neighbourEdges[el].size());

                    for (Index ne = 0; ne < neMax; ++ne)
                    {
                        _neighbours[_neighbourEdges[el][ne]].push_back(std::make_tuple(el, ne, 0));
                    }
                }

                // const auto & list = getElementList(1);
                for (auto el : getElementList(1))
                {
                    Index edge = _neighbourEdges[el][0];

                    for (auto &tup : _neighbours[edge])
                    {
                        Index face = std::get<0>(tup);
                        Index ne = std::get<1>(tup);

                        getOrientation(_elements[el]->getEdge(0), _elements[face]->getEdge(ne), std::get<2>(tup));
                    }
                }
            }
        }

        {
            std::vector<std::map<Face, Index>> mapFaces(getNumNodes());
            std::pair<std::map<Face, Index>::iterator, bool> mapFacesIt;

            _neighbourFaces.resize(_elements.size());

            // We use an antermediate numbering to speed up the algorithm
            for (Index i = 0; i < _nodes.size(); ++i)
                _nodes[i]->setIndex(i);

            for (Index el = 0; el < _elements.size(); ++el)
            {
                const Element &pe = *_elements[el];

                _neighbourFaces[el].resize(pe.getNumFaces());

                // Add faces if they do not already exist
                // 3D elements are referenced
                for (Index i = 0; i < pe.getNumFaces(); ++i)
                {
                    Face face = pe.getFace(i);
                    Index num_Coord = face.getSortedVertex(0).getIndex(); // sorted by their first Coord

                    // Insert the face in the map
                    mapFacesIt = mapFaces[num_Coord].insert(std::make_pair(face, -1));
                    // The edge is numeroted if it has effectivly been inserted
                    if (mapFacesIt.second == true)
                        mapFacesIt.first->second = _numFaces++;
                    // Reference the neighbour relation
                    _neighbourFaces[el][i] = mapFacesIt.first->second;
                }
            }

            // Recover the true numbering
            //(it must be done here because of the following computations of the orientation)
            for (Index i = 0; i < _nodes.size(); ++i)
                _nodes[i]->setIndex(_tmp_numbering[i]);

            // Construct the array Neighbours
            if (_dimElement == 3)
            {
                _neighbours.resize(_numFaces);
                for (auto el : getElementList(3))
                {
                    assert(_neighbourFaces[el].size() < std::numeric_limits<Index>::max());

                    Index nfMax = (Index)(_neighbourFaces[el].size());

                    for (Index nf = 0; nf < nfMax; ++nf)
                    {
                        _neighbours[_neighbourFaces[el][nf]].push_back(std::make_tuple(el, nf, 0));
                    }
                }

                for (auto el : getElementList(2))
                {
                    Index face = _neighbourFaces[el][0];

                    for (auto &tup : _neighbours[face])
                    {
                        Index vol = std::get<0>(tup);
                        Index nf = std::get<1>(tup);

                        getOrientation(_elements[el]->getFace(0), _elements[vol]->getFace(nf), std::get<2>(tup));
                    }
                }
            }
        }
    }

    inline const vector<Index> &Mesh::getNeighbourInterfaces(Index el) const
    {
        if (_dimElement == 3)
        {
            assert(el < _neighbourFaces.size());
            return _neighbourFaces[el];
        }

        if (_dimElement == 2)
        {
            assert(el < _neighbourEdges.size());
            return _neighbourEdges[el];
        }

        assert(el < _neighbourVertices.size());
        return _neighbourVertices[el];
    }

    inline vector<Index> Mesh::getNeighbourElements(Index el) const
    {
        assert(_elements[el]->getDim() == _dimElement);

        vector<Index> neighbourElements;

        vector<Index> interfaces = getNeighbourInterfaces(el);

        for (Index i = 0; i < interfaces.size(); ++i)
        {
            auto neighbours = getNeighbours(interfaces[i]);

            for (Index n = 0; n < neighbours.size(); ++n)
                neighbourElements.push_back(std::get<0>(neighbours[n]));
        }

        sort(neighbourElements.begin(), neighbourElements.end());
        auto it = unique(neighbourElements.begin(), neighbourElements.end());
        neighbourElements.erase(it, neighbourElements.end());
        it = remove(neighbourElements.begin(), neighbourElements.end(), el);
        neighbourElements.erase(it, neighbourElements.end());

        return neighbourElements;
    }

    inline void Mesh::getNeighbourElementsByInterface(Index el, std::vector<bool> &hasNeighbour,
                                                      std::vector<Index> &neighbourElement) const
    {
        assert(_elements[el]->getDim() == _dimElement);

        vector<Index> neighbourElements;
        vector<Index> interfaces = getNeighbourInterfaces(el);

        if (_elements[el]->getDim() == 3)
        {
            hasNeighbour.assign(_elements[el]->getNumFaces(), false);
            neighbourElement.assign(_elements[el]->getNumFaces(), 0);
        }
        else if (_elements[el]->getDim() == 2)
        {
            hasNeighbour.assign(_elements[el]->getNumEdges(), false);
            neighbourElement.assign(_elements[el]->getNumEdges(), 0);
        }
        else if (_elements[el]->getDim() == 1)
        {
            hasNeighbour.assign(2, false);
            neighbourElement.assign(2, 0);
        }

        for (Index i = 0; i < interfaces.size(); ++i)
        {
            auto neighbours = getNeighbours(interfaces[i]);

            if (neighbours.size() != 1)
            {
                hasNeighbour[i] = true;

                if (get<0>(neighbours[0]) == el)
                    neighbourElement[i] = std::get<0>(neighbours[1]);
                else
                    neighbourElement[i] = std::get<0>(neighbours[0]);
            }
        }
    }

    inline Mesh Mesh::extract(const std::vector<Index> &elems) const
    {
        // array of nodes
        std::vector<std::shared_ptr<Point>> new_nodes;

        // Equals true if the Coord has been added to the liste
        std::vector<bool> added(_nodes.size(), false);
        std::vector<Index> num(_nodes.size());
        Index counter = 0;

        // array of elements
        std::vector<std::shared_ptr<Element>> new_elements;

        // Temporary vectors of points
        std::vector<std::shared_ptr<Point>> vertices_tmp;

        for (Index i = 0; i < elems.size(); ++i)
        {
            std::shared_ptr<Element> elements_tmp = _elements[elems[i]];

            vertices_tmp.resize(elements_tmp->getNumNodes());

            for (Index j = 0; j < elements_tmp->getNumNodes(); ++j)
            {
                if (added[elements_tmp->getNode(j).getIndex()] == false)
                {
                    Real x = elements_tmp->getNode(j).x();
                    Real y = elements_tmp->getNode(j).y();
                    Real z = elements_tmp->getNode(j).z();

                    std::shared_ptr<Point> coord =
                        std::make_shared<Point>(x, y, z, elements_tmp->getNode(j).getIndex());

                    num[elements_tmp->getNode(j).getIndex()] = counter;

                    new_nodes.push_back(coord);
                    counter++;

                    vertices_tmp[j] = coord;

                    added[elements_tmp->getNode(j).getIndex()] = true;
                }
                else
                {
                    vertices_tmp[j] = new_nodes[num[elements_tmp->getNode(j).getIndex()]];
                }
            }

            new_elements.push_back(Factory(elements_tmp->getType(), vertices_tmp, elements_tmp->getLabel()));
        }

        Mesh new_mesh(new_nodes, new_elements);

        return new_mesh;
    }

} // namespace OndoMathX
