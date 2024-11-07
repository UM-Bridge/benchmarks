#pragma once

#include <array>
#include <cassert>
#include <memory>
#include <cmath>
#include "../../Utility/Defines.h"
#include "../Topology/Edge.h"
#include "../Topology/Face.h"
#include "../Topology/Point.h"
#include "../../Utility/ArrayAlgebra.h"
#include "Element.h"

namespace OndoMathX {
    
        //Implementation of the CRTP pattern
        template<class DerivedT,class ReferenceElement>
        class ElementT : public Element
        {
        protected:
            
            Index _label;
            
        public:
            
            explicit ElementT(Index label=0) : _label(label) {}
            virtual ~ElementT() {}
            
            
            Index getLabel() const override { return _label; }
            void setLabel(Index num) override { _label = num; }
            
            void getGradDef(const RealVector&, RealMatrix3x3&) const override;
            void getGradDef(const RealVector&, RealMatrix2x2&) const override;
            void getGradDef(const RealVector&, Real&) const override;
            
            // Interpolate the nodal data at point (u,v,w) in parametric coordinates
            void uvw2xyz(const RealVector&, RealVector&) const override;
            
            // invert the parametrisation
            bool xyz2uvw(const RealVector&, RealVector&) const override;
            
            //Return the barycenter of the element
            RealVector getBarycenter() const override;
            
            //Return the radius of the circonscript circle
            Real getRadius() const override;
            
            //Check if a coordinate is inside the element
            bool isInside(const RealVector &pos) const override;
            
            RefElement getRefType() const override {return ReferenceElement::getRefType();}
            Index getNumEdges() const override {return ReferenceElement::getNumEdges();}
            Index getNumFaces() const override {return ReferenceElement::getNumFaces();}
            Index getNumVertices() const override {return ReferenceElement::getNumVertices();}
            Index getDim() const override {return ReferenceElement::getDim();}
            
        };
        
} // OndoMathX



#include "ElementT.hxx"
 
