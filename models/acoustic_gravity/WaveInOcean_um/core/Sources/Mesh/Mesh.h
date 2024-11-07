#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <string>
#include <limits>
#include <vector>
#include <list>
#include <queue>
#include <map>
#include <tuple>
#include <algorithm>

#include "Element/Element.h"
#include "Element/Factory.h"
 
#include <cassert>

namespace OndoMathX {
    
 
		using namespace std;
	   
        /**
         * \class Mesh
         * Describing the sets of methods and attributes for a mesh.
         */
        
        class Mesh
        {
        private:
            
			
			//Data : vertices and elements ///////////////////////////
            //!array of vertices
            vector<shared_ptr<Point>> _nodes;
            //!array of elements
            vector<shared_ptr<Element>> _elements; 
            // access first element of _elements:
            // _elements[0] == shared_ptr<Element>
            // *_elements[0] == Element

            ///////////////////////////////////////////////////////////
			
			
			//Topology ///////////////////////////////////////////////
            /**If _dim_element == 2 : => 2D interfaces : each cell correspond to an edge and contain its neighbour elements.
             * If _dim_element == 3 : => 3D interfaces : each cell correspond to a face and contain its neighbour elements.
             * The Index corresponds to the number of the element
             * The first Index value corresponds to the number of the interface from the point of view of the element.
             * The second Index value corresponds to the Orientation mapping the interface to its element interface.
             */
            vector<vector<tuple<Index,Index,Index>>> _neighbours;
            //!Each cell corresponds to a volume and contains the neighbour faces
            vector<vector<Index> > _neighbourFaces;
            Index _numFaces;
            //!Each cell corresponds to a face and contains the neighbour edges
            vector<vector<Index> > _neighbourEdges;
            Index _numEdges;
            //!Each cell corresponds to an edge and contains the neighbour vertices
            std::vector< std::vector<Index> > _neighbourVertices;
            Index _numVertices;
			//////////////////////////////////////////////////////////
			
			
			//Dimensions ////////////////////////////////////////////
            //!Dimension of the vertices (if all vertices has a Z-component to 0 then _dim < 3, etc..)
            Index _dim;
            //!Highest dimension among the elements
            Index _dimElement;
            //////////////////////////////////////////////////////////

            
			//Lists of elements sorted by dimension and label/////////
            //!Element indices sorted by dimension
            vector<vector<Index>> _elementsByDim;
            //!Elements indices sorted by dimension and label
            vector<map<Index,vector<Index>>> _elementsByDimAndLabel;
            //!Labels
            vector<Index> _labels;
            vector<vector<Index>> _labelsByDim;
            //!Empty vector used when returning a void list
            vector<Index> _emptyVector;
            /////////////////////////////////////////////////////////
		
            
            //!Clear the link between the neighbours
            void _clearLink();
            
            //!Create the link between the neighbours
            void _link();
            
            //!Function that copies the _nodes and _elements of an other mesh
            void _copy(const Mesh &);
            
            //!Create elements sublists
            void _initElementSubList();
            
            //! Init the dimension variable
            void _initDim();
            
            //!Clean the mesh (deallocate)
            void _clean();
            
        public:
            
            //!A constructor based upon .msh or .mesh files
            /**
             * @param fileName a pointer to the name of the file in which is contained the mesh
             * @param format MeshFormat
             */
            explicit Mesh(const char *fileName, MeshFormat format);
            
            //! Default constructor
            explicit Mesh();
            
            //! Explicit constructor of the mesh, from the set of vertices, elements and the dimension.
            explicit Mesh(const vector<shared_ptr<Point> > &nodes,
                          const vector<shared_ptr<Element> > &elements);
            
            //! Constructor for the duplication of a mesh.
            Mesh(const Mesh &);
            
            Mesh & operator=(const Mesh &);
            
            ~Mesh();
            
            //! Returns a list of element of the corresponding dimension
            /**
             *@param dim dimension of the elements
             *@return a list of indices on the elements
             */
            const vector<Index> & getElementList(Index dim) const;
            
            //! Returns a list of element of the corresponding dimension and label
            /**
             *@param dim dimension of the elements
             *@param label label of the elements
             *@return a list of indices on the elements
             */
            const vector<Index> & getElementList(Index dim, Index label) const;
            
            //! Returns the list of all labels for a given dimension
            /**
             *@param dim dimension of the elements
             *@return a list of labels
             */
            const vector<Index> & getLabels(Index dim) const;
            
            const vector<Index> & getLabels() const;
        
            void setLabels(const std::vector<Index> &);
            
            //!Create a sub-mesh from a list of elements of the mesh
            /**
             *@param elems the list of elements of the sub-mesh
             */
            Mesh extract(const std::vector<Index> &elems) const;
            
            
            //!Create a mesh where each elements of the actual mesh are subdivided
            Mesh subDivide() const;



            //! Transform the mesh using the a function
            /**
             * @param f function that modifies its parameter coordinates
             */
            void transform(std::function<void(Point &)> f) {for (auto v : _nodes) f(*v);}
            
            // Set all z coordinates to 0
            void makePlane()
            {
                assert(_dimElement < 3);
                
                for (auto & coord : _nodes) (*coord)[2]=0.0;
                
                _dim = 2;
            };
            
            //! Returns the number of elements of the mesh
            /**
             *@return the number of elements of the mesh
             */
            Index getNumElements() const {return _elements.size();}
            
            //! Access the attribute _dim_element
            /**
             *@return the attribute \a _dim_element
             */
            Index getDimElement() const {return _dimElement;}
            
            //! Access the attribute _dim
            /**
             *@return the attribute \a _dim
             */
            Index getDim() const {return _dim;}
            
            //! Returns the element number \a el
            /**
             *@param el index of an element
             *@return the  element \a el
             */
            const Element & getElement(Index el) const
            {
                assert (el<_elements.size());
                return *_elements[el];
            }
            
            //! Returns the neighbour vertices of the element \a el
            /**
             *@param el index of an element
             *@return the neighbour vertices of the element \a el
             */
            const vector<Index>& getNeighbourVertices(Index el) const
            {
                assert (el<_neighbourVertices.size());
                return _neighbourVertices[el];
            }
            
            //! Returns the neighbour edges of the element \a el
            /**
             *@param el index of an element
             *@return the neighbour edges of the element \a el
             */
            const vector<Index>& getNeighbourEdges(Index el) const
            {
                assert (el<_neighbourEdges.size());
                return _neighbourEdges[el];
            }
            
            //! Returns the neighbour faces of the element \a el
            /**
             *@param el index of an element
             *@return the neighbour faces of the element \a el
             */
            const vector<Index>& getNeighbourFaces(Index el) const
            {
                assert (el<_neighbourFaces.size());
                return _neighbourFaces[el];
            }
            
            //! Returns the neighbour interface of the element \a el
            /**
             *@param el index of an element
             *@return the neighbour interfaces of the element \a el
             */
            const vector<Index>& getNeighbourInterfaces(Index el) const;

            
            //! Returns the neighbour elements of the same dimension \a el
            /**
             *@param el index of an element
             *@return the neighbour elements of the element \a el
             */
            vector<Index> getNeighbourElements(Index el) const;
            
            //! Returns the neighbour elements of the same dimension \a el
            /**
             *@param el index of an element
             *@return the neighbour elements of the element \a el
             */
            void getNeighbourElementsByInterface(Index el, std::vector<bool> & hasNeighbour,
                                                 std::vector<Index> & neighbourElement) const;
            
            
            //! Returns the number of nodes (vertices + points for the geometric description) of the mesh
            /**
             *@return the number of nodes of the mesh
             */
            Index getNumNodes() const {return static_cast<Index>(_nodes.size());}
            
            
            //! Returns the Coord number \a v
            /**
             *@param v index of a Coord
             *@return the Coord number \a v
             */
            const Point & getNode(Index v) const
            {
                assert (v<_nodes.size());
                return *_nodes[v];
            }

            //! Returns the Coord number \a v
            /**
             *@param v index of a Coord
             *@return the Coord number \a v
             */
            Point & getNode(Index v)
            {
                assert (v<_nodes.size());
                return *_nodes[v];
            }

            //! Returns the number of edges of the mesh
            /**
             *@return the number of edges of the mesh
             */
            Index getNumEdges() const {return _numEdges;}
            
            //! Returns the number of faces of the mesh
            /**
             *@return the number of faces of the mesh
             */
            Index getNumFaces() const {return _numFaces;}
            
            //! Returns the number of faces of the mesh
            /**
             *@return the number of faces of the mesh
             */
            Index getNumVertices() const {return _numVertices;}
            
            
            const vector<tuple<Index,Index,Index>> &
            getNeighbours(Index elem) const
            {
                assert (elem<_neighbours.size());
                return _neighbours[elem];
            }

        };
        
} //namespace Ondomatic

#include "Utility/Intersection.hxx"
#include "Utility/Rasterization.hxx"
#include "Utility/Coloring.hxx"
#include "Utility/SubDivide.hxx"
#include "Utility/Numbering/Numbering.hxx"
#include "Utility/Other.hxx"

#include "IO/IO_MESH.hxx"
#include "IO/IO_MSH.hxx"
#include "IO/IO_VTK.hxx"

#include "Mesh.hxx"

