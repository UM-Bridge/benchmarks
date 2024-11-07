#pragma once

#include "ElementT.h"
#include "../ReferenceElement/RSegment.h"

namespace OndoMathX
{

class Segment2 : public ElementT<class Segment2,RSegment>
{
    
protected:
    
    std::array< std::shared_ptr<Point>, 2> _v;
    
    public :
    
    explicit Segment2(std::vector< std::shared_ptr<Point> > &v, Index label=0)
    : ElementT<class Segment2,RSegment>(label)
    {
        _v[0] = v[0];
        _v[1] = v[1];
    }
    
    Segment2() = delete;
    explicit Segment2(const Segment2& rhs) = default;
    explicit Segment2(Segment2&& rhs) = default;
    Segment2& operator=(const Segment2& rhs) = default;
    virtual ~Segment2();
    
    inline GeoElement getType() const { return GeoSegment2; }
    
    inline Index getNumVertices() const { return 2; }
    inline const Point & getVertex(Index num) const
    {
        assert(num<2);
        return *_v[num];
    }
    /*
    inline std::shared_ptr<Point> getVertexPtr(Index num)
    {
        assert(num < 2);
        return _v[num];
    }*/
    
    inline Index getNumNodes() const { return getNumVertices(); };
    inline const Point & getNode(Index num) const
    {
        return getVertex(num);
    };
    
    
    inline void getShapeFunction(Index num,const RealVector& uvw, Real &s) const
    {
        assert(num<2);
        
        if (num == 0)
        {
            s = (1. - uvw[0]);
        }
        else
        {
            s = uvw[0];
        }
        
    }
    
    inline void getGradShapeFunction(Index num,const RealVector& /*uvw*/, RealVector &s) const
    {
        assert(num<2);
        
        if (num == 0)
        {
            s[0] = -1.0;
            s[1] = 0.;
            s[2] = 0.;
        }
        else
        {
            s[0] = 1.0;
            s[1] = 0.;
            s[2] = 0.;
        }
    }
    
    inline Edge getEdge(Index num) const
    {
        assert(num==0);
        return Edge(_v[0],_v[1]);
    }
    
    inline Face getFace(Index /*num*/) const
    {
        assert(false);
        return Face(_v[0],_v[0],_v[0]);
    }
    
    inline Index getOrientationEdge(Index e) const
    {
        assert(e==0);
        
        if ( _v[0]->getIndex()>_v[1]->getIndex())
            return 1;
        else  return 0;
    }
    
    inline Index getOrientationFace(Index /*f*/) const { assert(false); return 0;}
    
    
    inline bool isAffine() const
    {
        return true;
    }
    
    
    bool xyz2uvw(const RealVector& xyz, RealVector& uvw) const;
    
    /*
     virtual void refine(int option, const vector<Real> & ref_pos,
     vector<MVertex*> > 		 & points_nodes,
     vector<vector<MVertex*> > & edges_nodes,
     vector<vector<MVertex*> > & faces_nodes,
     vector<MElement*> & new_elements,
     vector<MVertex*>  & new_nodes);
     */
    
};


inline Segment2::~Segment2() = default;



// invert the parametrisation
inline bool Segment2::xyz2uvw(const RealVector& xyz, RealVector& uvw) const
{
    uvw[1]=0.0;
    uvw[2]=0.0;
    
    //Construct a system of the form : S u = F
    RealVector c;
    RealVector S;
    RealVector F;
    
    F[0]=xyz[0];
    F[1]=xyz[1];
    F[2]=xyz[2];
    
    //First shape function = 1 - u
    c[0] = _v[0]->x();
    c[1] = _v[0]->y();
    c[2] = _v[0]->z();
    
    for (int i=0;i<3;++i)
    {
        F[i]-=c[i];
        S[i]=-c[i];
    }
    
    //Second shape function = u
    c[0] = _v[1]->x();
    c[1] = _v[1]->y();
    c[2] = _v[1]->z();
    
    for (int i=0;i<3;++i)
    {
        S[i]+=c[i];
    }
    
    //Solve the linear problem
    Real STF=ArrayAlgebra::ScalarProduct(S,F);
    Real STS=ArrayAlgebra::ScalarProduct(S,S);
    
    uvw[0]=STF/STS;
    
    return true;
}

} // OndoMathX





