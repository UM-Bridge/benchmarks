#pragma once

#include <set>
#include <vector>
#include <algorithm>
#include <cassert>
#include <map>

 
namespace OndoMathX {
    

        
        template<class ReferenceElement>
        class LocalDofSet
        {
        protected:
        
            std::set<Point,lessCoordinatesXYZ> _dofSet;
            std::vector<Point> _dofList;
            
            //3D vectors : (vertex,edge,face) numbers x Orientation numbers x list of dof indexes
            std::vector<std::vector<Index> > _vertexDofs;
            
            std::vector<std::vector<std::vector<Index> > > _interiorEdgeDofs;
            
            std::vector<std::vector<std::vector<Index> > > _interiorFaceDofs;
            
            std::vector<std::vector<std::vector<Index> > > _edgeDofs;
            
            std::vector<std::vector<std::vector<Index> > > _faceDofs;
            
            
            std::vector<Index> _volumeDofs;
            
            ReferenceElement _refElement;
            
            bool _pos2LocalDof(const RealVector & pos, Index & dofIndex);
            
            
        public:
            
            LocalDofSet(vector<RealVector> DofSet);
            
        
            //////////////////////////////////////////////////////////////////////////////////
            const std::vector<Index> &
            getVertexDofs(Index v) const
            {
                return _vertexDofs[v];
            }
            
            const std::vector<Index> &
            getInteriorEdgeDofs(Index e, Index t) const
            {
                return _interiorEdgeDofs[e][t];
            }
            
            const std::vector<Index> &
            getInteriorFaceDofs(Index f, Index t) const
            {
                return _interiorFaceDofs[f][t];
            }
            
            const std::vector<Index> &
            getEdgeDofs(Index e, Index t) const
            {
                return _edgeDofs[e][t];
            }
            
            const std::vector<Index> &
            getFaceDofs(Index f, Index t) const
            {
                return _faceDofs[f][t];
            }
            
            
            const std::vector<Index> &
            getVolumeDofs() const
            {
                return _volumeDofs;
            }
            
            //////////////////////////////////////////////////////////////////////////////////
             
            
            //////////////////////////////////////////////////////////////////////////////////
            Index getNumFaces() const {return _refElement.getNumFaces();}
            Index getNumEdges() const {return _refElement.getNumEdges();}
            Index getNumVertices() const {return _refElement.getNumVertices();}
            Index getDim() const {return _refElement.getDim();}
            RefElement getRefType() const {return _refElement.getRefType();}
            //////////////////////////////////////////////////////////////////////////////////
            
            //! Print function.
            void print(std::ostream& out);
            
        };
        
    
} //OndoMathX



#include "LocalDofSet.hxx"


 




