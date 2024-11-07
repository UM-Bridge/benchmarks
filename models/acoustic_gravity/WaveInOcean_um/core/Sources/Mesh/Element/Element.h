#pragma once

#include <array>
#include <cassert>
#include "../../Utility/Defines.h"
#include "../Topology/Edge.h"
#include "../Topology/Face.h"
#include "../Topology/Point.h"



namespace OndoMathX {
        
        /**
         * \class Element
         * Parent class for the reference elements.
         * Virtual class describing the set of methods
         * for a reference element.
         */
        class Element
        {
     
        public:
            
            
            virtual ~Element();
            
            // get the edges
            virtual Index getNumEdges() const = 0;
            virtual Edge getEdge(Index num) const = 0;
            
            // get the faces
            virtual Index getNumFaces() const = 0;
            virtual Face getFace(Index num) const = 0;
            
            // get the vertices (used for connectivity)
            virtual Index getNumVertices() const = 0;
            virtual const Point & getVertex(Index num) const = 0;

            // get all the nodes (used to describe the deformation of the element)
            virtual Index getNumNodes() const = 0;
            virtual const Point & getNode(Index num) const = 0;
            
            // returns the geometrical dimension of the element
            virtual Index getDim() const = 0;
            
            // get/set the partition to which the element belongs
            virtual Index getLabel() const = 0;
            virtual void setLabel(Index num) = 0;
            

            // Returns the interpolating nodal shape function associated with
            // node num, evaluated at point (u,v,w) in parametric coordinates
            virtual void getShapeFunction(Index num,const RealVector&,Real &) const = 0;
            
            // Returns the gradient of of the nodal shape function associated
            // with node num, evaluated at point (u,v,w) in parametric
            // coordinates
            virtual void getGradShapeFunction(Index num,
                                              const RealVector&,
                                              RealVector&) const = 0;
            
            // Returns the Jacobian of the element evaluated at point (u,v,w) in
            // parametric coordinates
            virtual void getGradDef(const RealVector&, RealMatrix3x3 &) const = 0;
            virtual void getGradDef(const RealVector&, RealMatrix2x2 &) const = 0;
            virtual void getGradDef(const RealVector&, Real &) const = 0;
            
            // Interpolate the nodal data at point (u,v,w) in parametric coordinates
            virtual void uvw2xyz(const RealVector&, RealVector&) const = 0;
            
            // invert the parametrisation
            virtual bool xyz2uvw(const RealVector&, RealVector&) const = 0;
            
            //Access the type and Topology of the element (using the corresponding enum)
            virtual GeoElement getType() const = 0;
            virtual RefElement getRefType() const = 0;
            
            //Return the barycenter of the element
            virtual RealVector getBarycenter() const = 0;
            
            //Return the radius of the circonscript circle
            virtual Real getRadius() const = 0;
            
            //Get the number of the Orientation for the edges or the faces
            virtual Index getOrientationEdge(Index e) const = 0;
            virtual Index getOrientationFace(Index f) const = 0;
            
            virtual bool isAffine() const = 0;
            
            //Check if a coordinate is inside the element
            virtual bool isInside(const RealVector &pos) const = 0;
            
            
            Index getOrientationInterface(Index i) const
            {
                if (getDim()==3) return getOrientationFace(i);
                if (getDim()==2) return getOrientationEdge(i);
                return 0;
            }
            
            Index getOrientation() const
            {
                if (getDim()==2) return getOrientationFace(0);
                if (getDim()==1) return getOrientationEdge(0);
                return 0;
            }
            
        };
        
        inline Element::~Element() = default;
 
}

 
