#pragma once

#include "LocalDofSet.h"

namespace OndoMathX
{

template<class ReferenceElement>  Index Numbering(const Mesh &mesh, std::vector<RealVector> DofSet, std::vector<std::vector<Index>> & loc2glob)
{
    LocalDofSet<ReferenceElement> localDofSet(DofSet);
    
    
    //Clear and resize the previous numbering of each operators
    loc2glob.resize(mesh.getNumElements());
    
    
    //Init of the nodes list for each numbering
    //The access uses the infos (nodeId)-> array of dof number
    std::vector< std::vector <Index> > verticesNumbers;
    std::vector< std::vector <Index> > edgesNumbers;
    std::vector< std::vector <Index> > facesNumbers;
    
    Index NDoF = 0;
    
 
    
    //Loop over the element to be numbered
    for (Index el = 0; el < mesh.getNumElements()  ; el++ )
    {
        //Recover a reference on the local-to-global array that is being built
        const Element & curElem =  mesh.getElement(el);
        
        if (curElem.getDim() == ReferenceElement::getDim())
        {
            loc2glob[el].resize(DofSet.size());
            
            //vertices
            for (Index v=0;v<ReferenceElement::getNumVertices();++v)
            {
                //recover the number of the node
                Index k = mesh.getNeighbourVertices(el)[v];
                
                //Allocate enough space to deal with the nodes infos.
                if (verticesNumbers.size()==0)
                    verticesNumbers.resize(mesh.getNumVertices());
                
                const std::vector<Index> & loc2loc = localDofSet.getVertexDofs(v);
                Index localNumDofs = loc2loc.size();
                
                if (localNumDofs>0)
                {
                    
                    //create the numbering on the node (if needed)
                    if ( verticesNumbers[k].size()==0)
                    {
                        verticesNumbers[k].resize(localNumDofs);
                        for (Index n = 0;n<localNumDofs;++n)
                        {
                            verticesNumbers[k][n]=NDoF;
                            NDoF++;
                        }
                    }
                    
                    //compute the local2global
                    for (Index n = 0;n<localNumDofs;++n)
                    loc2glob[el][loc2loc[n]]=verticesNumbers[k][n];
                }
            } // end managment of vertices
            
            //interior of the edges
            if (ReferenceElement::getDim() >= 1)
                for (Index e=0;e<ReferenceElement::getNumEdges();++e)
            {
                //recover the number of the node
                Index k = mesh.getNeighbourEdges(el)[e];
                
                //recover the orientation
                Index t = curElem.getOrientationEdge(e);
                
                //For all variable that must be numbered together
                
                //Allocate enough space to deal with the nodes infos.
                if (edgesNumbers.size()==0)
                    edgesNumbers.resize(mesh.getNumEdges());
                
                const std::vector<Index> & loc2loc = localDofSet.getInteriorEdgeDofs(e,t);
                
                Index localNumDofs = loc2loc.size();
                
                if (localNumDofs>0)
                {
                    //create the numbering on the node (if needed)
                    if (edgesNumbers[k].size()==0)
                    {
                        edgesNumbers[k].resize(localNumDofs);
                        
                        for (Index n = 0;n<localNumDofs;n++)
                        {
                            edgesNumbers[k][n]=NDoF;
                            NDoF++;
                        }
                    }
                    
                    //compute the local2global
                    for (Index n = 0;n<localNumDofs;++n)
                    loc2glob[el][loc2loc[n]] = edgesNumbers[k][n];
                }
            }// end managment of the interior of the edges
            
            //interior of the faces
            if (ReferenceElement::getDim() >= 2)
                for (Index f=0;f<ReferenceElement::getNumFaces();++f)
            {
                //recover the number of the node
                Index k = mesh.getNeighbourFaces(el)[f];
                
                //recover the orientation
                Index t = curElem.getOrientationFace(f);
                
                //Allocate enough space to deal with the nodes infos.
                if (facesNumbers.size()==0)
                    facesNumbers.resize(mesh.getNumFaces());
                
                const std::vector<Index> & loc2loc = localDofSet.getInteriorFaceDofs(f,t);
                
                Index localNumDofs = static_cast<Index> (loc2loc.size());
                
                if (localNumDofs>0)
                {
                    //create the numbering on the node (if needed)
                    if (facesNumbers[k].size()==0)
                    {
                        facesNumbers[k].resize(localNumDofs);
                        
                        for (Index n = 0;n<localNumDofs;n++)
                        {
                            facesNumbers[k][n]=NDoF;
                            NDoF++;
                        }
                    }
                    
                    for (Index n = 0;n<localNumDofs;++n)
                    loc2glob[el][loc2loc[n]]=facesNumbers[k][n];
                }
                
            } // end managment of the interior of the faces
            
            
            //interior of the volume
            if (ReferenceElement::getDim() == 3)
            {
                auto & loc2loc = localDofSet.getVolumeDofs();
                
                for (Index n = 0;n<loc2loc.size();++n)
                {
                    loc2glob[el][loc2loc[n]]=NDoF;
                    NDoF++;
                }
            }
            
        }
    } // end of the loop of the element
    
   

    //Numbering of the elements on the boundary
    for (Index el = 0; el < mesh.getNumElements()  ; el++ )
    {
        //Recover a reference on the local-to-global array that is being built
        const Element & curElem =  mesh.getElement(el);
        
        if (curElem.getDim() == (ReferenceElement::getDim()-1))
        {
            
            Index face;
            
            if (curElem.getDim() == 1)
                face = mesh.getNeighbourEdges(el)[0];
            else face = mesh.getNeighbourFaces(el)[0];
            
            auto & neighbours = mesh.getNeighbours(face);
            
            //Make sure that there is at least one neighbouring element
            assert(neighbours.size() > 0);
            
            Index volume = std::get<0>(neighbours[0]);
            Index local_face = std::get<1>(neighbours[0]);
            Index orientation = std::get<2>(neighbours[0]);
            
            if (curElem.getDim() == 1)
            {
                auto & loc2loc = localDofSet.getEdgeDofs(local_face,orientation);
                
                loc2glob[el].resize(loc2loc.size());
                
                for (Index n = 0;n<loc2loc.size();++n)
                loc2glob[el][n] = loc2glob[volume][loc2loc[n]];
            }
            else
            {
                auto & loc2loc = localDofSet.getFaceDofs(local_face,orientation);
                
                loc2glob[el].resize(loc2loc.size());
                
                for (Index n = 0;n<loc2loc.size();++n)
                loc2glob[el][n] = loc2glob[volume][loc2loc[n]];
            }
            
        }
    }
    
    return NDoF;
}

} // namespace OndoMathX




