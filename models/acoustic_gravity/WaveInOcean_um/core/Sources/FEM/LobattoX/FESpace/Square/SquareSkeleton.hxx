#pragma once
 

// -----------------------------------------------------------------------------------------//
// Class definition
namespace OndoMathX {

 



template <typename GeoTransform> class SquareSkeleton
{

private:

    const Square<GeoTransform> & _Square;

    Index _NX;
    Index _NY;
    Index _P;

    PeriodicityType _PerType;

    GaussLobattoElement _GLE;

public:
   
    /*! \brief Definition of the dimension as a static variable.*/
    static const Index Dim = 2;
   
    SquareSkeleton() = delete; 

    SquareSkeleton(const Square<GeoTransform> &aSquare, PeriodicityType PerType = None) : _Square(aSquare)
    {
        assert(aSquare._PX == aSquare._PY);

        _P = aSquare._PX;
        _NX = aSquare._NX;
        _NY = aSquare._NY;

        _PerType = PerType;

        _GLE = GaussLobattoElement(_P+1,0,0);
    }

    // ---------------------------------------------------------------------------------//
    Index getNumDoFs() const {return _Square.getNumDoFs();}
    Index getNumLocDoFs() const { return _P+1;}
    Index getNumColors() const {return 1;}
    Index getNumElements(Index iColor) const {return getNumElements();}
    Index getNumElements() const 
    {
        if (_PerType == PeriodicityType::None) 
            return (_NX-1)*_NY+(_NY-1)*_NX;
        else if (_PerType == PeriodicityType::XY)
             return 2*_NX*_NY;

        return 0;
    }
    const GaussLobattoElement & getFE() const { return _GLE;}
    Index EltColorNum2EltGlobalNum(Index c, Index e) const {return e;}

    const Square<GeoTransform> & getFESpace() const {return _Square;}

    // ---------------------------------------------------------------------------------//

    void getNeighbours(Index iEltSkeleton, Index &iEltK, Index &iEltL) const
    {
        if (_PerType == PeriodicityType::None) 
        {
            if (iEltSkeleton < (_NX-1)*_NY)
            {
                Index I = iEltSkeleton % (_NX-1);
                Index J = iEltSkeleton / (_NX-1);

                iEltK = I + _NX*J;
                iEltL = iEltK + 1;
            }
            else
            {
                iEltSkeleton -= (_NX-1)*_NY;

                Index I = iEltSkeleton % _NX;
                Index J = iEltSkeleton / _NX;
                
                iEltK = I + _NX*J;
                iEltL = iEltK + _NX;
            }
        }
        else if (_PerType == PeriodicityType::XY) 
        {
            if (iEltSkeleton < _NX*_NY)
            {
                Index I = iEltSkeleton % _NX;
                Index J = iEltSkeleton / _NX;

                iEltK = I + _NX*J;

                if (I == (_NX-1)) iEltL = iEltK - (_NX-1);
                else
                    iEltL = iEltK + 1;

            }
            else
            {
                iEltSkeleton -= _NX*_NY;

                Index I = iEltSkeleton % _NX;
                Index J = iEltSkeleton / _NX;
                
                iEltK = I + _NX*J;
                
                if (J == (_NY-1)) iEltL = iEltK - (_NY-1)*_NX;
                else
                    iEltL = iEltK + _NX;
            }
        }
        

    }

    Index Loc2Loc(Index iEltSkeleton, Index iLoc, bool K = true) const
    {
        
        if (_PerType == PeriodicityType::None) 
        {
            if (iEltSkeleton < (_NX-1)*_NY)
            {
                if (K)
                {
                    return (_P+1)*iLoc + _P;
                }
                else
                {
                    return (_P+1)*iLoc;
                }
            }
            else
            {
                if (K)
                {
                    return (_P+1)*_P + iLoc;
                }
                else
                {
                    return iLoc;
                }
            } 
        }    
        else if (_PerType == PeriodicityType::XY) 
        {
            if (iEltSkeleton < _NX*_NY)
            {
                if (K)
                {
                    return (_P+1)*iLoc + _P;
                }
                else
                {
                    return (_P+1)*iLoc;
                }
            }
            else
            {
                if (K)
                {
                    return (_P+1)*_P + iLoc;
                }
                else
                {
                    return iLoc;
                }
            } 
        } 

        return 0;
    }

    void getNormal(Index iEltSkeleton, Index iLoc, std::array<Real, 2> &n) const
    {
        Index iEltK;
        Index iEltL;
        Index iLocK;

        getNeighbours(iEltSkeleton,iEltK,iEltL);
        iLocK = Loc2Loc(iEltSkeleton,iLoc);

        if (_PerType == PeriodicityType::None) 
        {
            if (iEltSkeleton < (_NX-1)*_NY)
            {
                n[0]=1.0;
                n[1]=0.0;
            }
            else
            {
                n[0]=0.0;
                n[1]=1.0;
            }
        }
        else if (_PerType == PeriodicityType::XY) 
        {
            if (iEltSkeleton < _NX*_NY)
            {
                n[0]=1.0;
                n[1]=0.0;
            }
            else
            {
                n[0]=0.0;
                n[1]=1.0;
            }
        }
 
        RealMatrix2x2 F;
        RealMatrix2x2 CoMatF;

        _Square.getGradDef(iEltK,iLocK,F);

        ArrayAlgebra::CoMat(F,CoMatF);
        ArrayAlgebra::MatMlt(CoMatF,n);
    }

    void getDoFCoordinate(const Index iEltGlob, const Index iLoc, RealVector& xyz) const
    {
        Index iEltK;
        Index iEltL;
        Index iLocK;

        getNeighbours(iEltGlob,iEltK,iEltL);
        iLocK = Loc2Loc(iEltGlob,iLoc);

        _Square.getDoFCoordinate(iEltK,iLocK,xyz);  
    }
    
    
};


 
}//OndoMathX
 
