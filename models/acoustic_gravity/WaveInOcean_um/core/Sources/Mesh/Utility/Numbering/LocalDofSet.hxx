

namespace OndoMathX
{

    template <class ReferenceElement>
    LocalDofSet<ReferenceElement>::LocalDofSet(vector<RealVector> DofList)
    {

        // Initialisation of basic informations/////////////////////////////////////////////
        // We get the number of interfaces
        Index numFaces = ReferenceElement::getNumFaces();
        Index numEdges = ReferenceElement::getNumEdges();
        Index numVert = ReferenceElement::getNumVertices();
        Index numOrientationEdge = 2;
        std::vector<Index> numOrientationFace(ReferenceElement::getNumFaces(), 1);

        _dofList.resize(DofList.size());

        // First we construct the list of dofs for the specific finite element
        for (Index i = 0; i < _dofList.size(); ++i)
        {
            // Set the index into increasing order in the list
            _dofList[i] = Point(DofList[i], i);

            // We check if the element is inside the element
            assert(ReferenceElement::isInside(_dofList[i].coord()));

            // create the map coord to index.
            _dofSet.insert(_dofList[i]);
        }

        // Compute the number of Orientations of the interface (used only for 3D)
        for (Index i = 0; i < numFaces; ++i)
        {
            if (ReferenceElement::getFaceRefType(i) == Quadrangle)
                numOrientationFace[i] = 8;
            else if (ReferenceElement::getFaceRefType(i) == Triangle)
                numOrientationFace[i] = 6;
        }

        {
            _interiorFaceDofs.resize(numFaces);
            _faceDofs.resize(numFaces);

            for (Index f = 0; f < numFaces; f++)
            {
                _interiorFaceDofs[f].resize(numOrientationFace[f]);
                _faceDofs[f].resize(numOrientationFace[f]);
            }

            _interiorEdgeDofs.resize(numEdges);
            _edgeDofs.resize(numEdges);

            for (Index e = 0; e < numEdges; e++)
            {
                _interiorEdgeDofs[e].resize(numOrientationEdge);
                _edgeDofs[e].resize(numOrientationEdge);
            }

            _vertexDofs.resize(numVert);
        }
        //////////////////////////////////////////////////////////////////////////////////

        // Sort the nodes so that they are sorted with (x,y,z) increasing when stored in the array
        std::sort(_dofList.begin(), _dofList.end());

        // Insertion of the nodes in the appropiate arrays
        for (Index i = 0; i < _dofList.size(); i++)
        {
            RealVector pos;

            // If the node is shared and is on a vertex
            for (Index v = 0; v < numVert; v++)
                if (ReferenceElement::isOnVertex(_dofList[i].coord(), v))
                {
                    _vertexDofs[v].push_back(_dofList[i].getIndex());
                }

            // If the node is shared and is on the interior of an edge.
            for (Index e = 0; e < numEdges; ++e)
            {
                if (ReferenceElement::isOnEdge(_dofList[i].coord(), e))
                {
                    for (Index t = 0; t < numOrientationEdge; ++t)
                    {
                        ReferenceElement::posEdgeTransformed(_dofList[i].coord(), pos, e, t);

                        Index nodeIndex = 0;

                        if (_pos2LocalDof(pos, nodeIndex))
                        {
                            _edgeDofs[e][t].push_back(nodeIndex);

                            if (!ReferenceElement::isOnVertices(_dofList[i].coord()))
                                _interiorEdgeDofs[e][t].push_back(nodeIndex);
                        }
                        else
                        {
                            assert(false);
                        }
                    }
                }
            }

            // If the node is shared and is on the interior of a face
            for (Index f = 0; f < numFaces; f++)
            {
                if (ReferenceElement::isOnFace(_dofList[i].coord(), f))
                {
                    for (Index t = 0; t < numOrientationFace[f]; t++)
                    {
                        ReferenceElement::posFaceTransformed(_dofList[i].coord(), pos, f, t);

                        Index nodeIndex = 0;

                        if (_pos2LocalDof(pos, nodeIndex))
                        {
                            _faceDofs[f][t].push_back(nodeIndex);

                            if (!ReferenceElement::isOnEdges(_dofList[i].coord()))
                                _interiorFaceDofs[f][t].push_back(nodeIndex);
                        }
                        else
                        {
                            assert(false);
                        }
                    }
                }
            }

            // If the node is not Shared
            if (ReferenceElement::getDim() == 3 && !ReferenceElement::isOnFaces(_dofList[i].coord()))
            {
                _volumeDofs.push_back(_dofList[i].getIndex());
            }
        }
    }

    template <class ReferenceElement>
    bool LocalDofSet<ReferenceElement>::_pos2LocalDof(const RealVector &pos, Index &dofIndex)
    {
        Point coord(pos[0], pos[1], pos[2], 0);

        auto it_coord = _dofSet.find(coord);

        if (it_coord != _dofSet.end())
        {
            dofIndex = (*it_coord).getIndex();
            return true;
        }
        else
        {
            dofIndex = 0;
            return false;
        }
    }

    //! Print function.
    template <class ReferenceElement>
    void LocalDofSet<ReferenceElement>::
        print(std::ostream &out)
    {

        for (Index v = 0; v < _vertexDofs.size(); ++v)
        {
            if (_vertexDofs[v].size() != 0)
            {
                out << "vertex " << v << ": (";
                for (Index n = 0; n < _vertexDofs[v].size(); ++n)
                {
                    out << _vertexDofs[v][n];
                    if (n != _vertexDofs[v].size() - 1)
                        out << ",";
                }
                out << ")" << std::endl;
            }
        }

        for (Index e = 0; e < _interiorEdgeDofs.size(); ++e)
        {
            for (Index t = 0; t < _interiorEdgeDofs[e].size(); ++t)
            {
                if (_interiorEdgeDofs[e][t].size() != 0)
                {
                    out << "edge " << e << ", Orientation " << t << ": (";
                    for (Index n = 0; n < _interiorEdgeDofs[e][t].size(); ++n)
                    {
                        out << _interiorEdgeDofs[e][t][n];
                        if (n != _interiorEdgeDofs[e][t].size() - 1)
                            out << ",";
                    }
                    out << ")" << std::endl;
                }
            }
        }

        for (Index f = 0; f < _interiorFaceDofs.size(); ++f)
        {
            for (Index t = 0; t < _interiorFaceDofs[f].size(); ++t)
            {
                if (_interiorFaceDofs[f][t].size() != 0)
                {
                    out << "face " << f << ", Orientation " << t << ": (";
                    for (Index n = 0; n < _interiorFaceDofs[f][t].size(); ++n)
                    {
                        out << _interiorFaceDofs[f][t][n];
                        if (n != _interiorFaceDofs[f][t].size() - 1)
                            out << ",";
                    }
                    out << ")" << std::endl;
                }
            }
        }

        if (_volumeDofs.size() != 0)
        {
            out << "volume : (";
            for (Index n = 0; n < _volumeDofs.size(); ++n)
            {
                out << _volumeDofs[n];
                if (n != _volumeDofs.size() - 1)
                    out << ",";
            }
            out << ")" << std::endl;
        }
    }

} // OndoMathX
