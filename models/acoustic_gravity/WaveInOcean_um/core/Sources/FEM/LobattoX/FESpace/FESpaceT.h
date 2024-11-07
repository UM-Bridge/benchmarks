#pragma once


// Class definition
namespace OndoMathX {


    template<class FESpace, Index Dim> class FESpaceTSubdomain;

    template<class FESpace,Index Dim>
    class FESpaceT {

    public: 

        typedef std::array<std::array<Real,Dim>,Dim> RealMatrix;

        typedef GaussLobattoInterpolator Interpolator;

    protected:

        //buffer for pre-computations of the geometric informations and loc2glob
        std::vector<std::vector<RealMatrix>> _GradDef;
        std::vector<std::vector<Index>>      _loc2glob;

        bool _GeomInfoPreComputed;
        bool _Loc2GlobPreComputed;
        
    public: 

        // ---------------------------------------------------------------------------------//

        FESpaceT() 
        {
            _GeomInfoPreComputed = false;
            _Loc2GlobPreComputed = false;
        }


        void ComputeGeometricInfo() 
        {
            const FESpace & self = static_cast<const FESpace &>(*this);

            Index NElem = self.getNumElements();
            Index NLocDoF= self.getNumLocDoFs();

            _GradDef.resize(NElem);

            for (Index iEltGlob=0;iEltGlob<NElem;++iEltGlob)
            {
                _GradDef[iEltGlob].resize(NLocDoF);
    
                for (Index iLoc=0;iLoc<NLocDoF;++iLoc)
                {
                    self._getGradDef(iEltGlob,iLoc,_GradDef[iEltGlob][iLoc]);
                } 
            }
            
            _GeomInfoPreComputed = true;
        }
        
        void getGradDef(Index iEltGlob, Index iLoc,
			    RealMatrix& GradDef) const
        {          
            if (_GeomInfoPreComputed)
            {
                GradDef = _GradDef[iEltGlob][iLoc];
            }
            else
            {
                const FESpace & self = static_cast<const FESpace &>(*this);

                self._getGradDef(iEltGlob,iLoc,GradDef);
            }
        }

        void ComputeLoc2Glob()
        {
            const FESpace & self = static_cast<const FESpace &>(*this);

            Index NElem = self.getNumElements();
            Index NLocDoF= self.getNumLocDoFs();
         
            _loc2glob.resize(NElem);
            
            for (Index iEltGlob=0;iEltGlob<NElem;++iEltGlob)
            {
                _loc2glob[iEltGlob].resize(NLocDoF);
    
                for (Index iLoc=0;iLoc<NLocDoF;++iLoc)
                    _loc2glob[iEltGlob][iLoc] = self._Loc2Glob(iEltGlob,iLoc);
                
            }
            
            _Loc2GlobPreComputed = true;
        }

 
        Index Loc2Glob(Index iEltGlob, Index iLoc) const
        {
            if (_Loc2GlobPreComputed)
            {
                return _loc2glob[iEltGlob][iLoc];
            }
            else
            {
                const FESpace & self = static_cast<const FESpace &>(*this); 

                return self._Loc2Glob(iEltGlob,iLoc);
            }
        } 

        void Loc2Glob(Index iEltGlob, std::vector<Index> & loc2glob) const
        {
            if (_Loc2GlobPreComputed)
            {
                loc2glob = _loc2glob[iEltGlob];
            }
            else
            {
                const FESpace & self = static_cast<const FESpace &>(*this); 

                Index NLocDoF= self.getNumLocDoFs();
                
                loc2glob.resize(NLocDoF);
            
                for (Index iLoc=0; iLoc < NLocDoF; iLoc++)
                    loc2glob[iLoc] = self._Loc2Glob(iEltGlob,iLoc);
            }
        }


        Index Loc2GlobDisc(Index iEltGlob, Index iLoc) const
        {
            const FESpace & self = static_cast<const FESpace &>(*this); 
            Index NLocDoF= self.getNumLocDoFs();

            return iEltGlob*NLocDoF+iLoc;
        }
        

        void Loc2GlobDisc(Index iEltGlob, std::vector<Index> & loc2glob) const
        {
            const FESpace & self = static_cast<const FESpace &>(*this); 

            Index NLocDoF= self.getNumLocDoFs();
                
            loc2glob.resize(NLocDoF);
            
            for (Index iLoc=0; iLoc < NLocDoF; iLoc++)
                loc2glob[iLoc] = iEltGlob*NLocDoF+iLoc;
        }

        // ---------------------------------------------------------------------------------//

		Index getNumDoFsDisc() const
		{
            const FESpace & self = static_cast<const FESpace &>(*this); 

            Index NElem = self.getNumElements();
            Index NLocDoF= self.getNumLocDoFs();

            return NElem*NLocDoF;
		}

        Real getJacobian(Index iEltGlob, Index iLoc) const
        {
            RealMatrix GradDef;
            getGradDef(iEltGlob,iLoc,GradDef);
            return ArrayAlgebra::Det(GradDef);
        }

        void MltCoMatGradDef(Index iEltGlob, Index iLoc, std::array<Real,Dim> &U) const
		{
            RealMatrix GradDef;
            RealMatrix CoMatGradDef;
            getGradDef(iEltGlob,iLoc,GradDef);
            ArrayAlgebra::CoMat(GradDef,CoMatGradDef);
            ArrayAlgebra::MatMlt(CoMatGradDef,U);
        }

        void TransposeMltCoMatGradDef(Index iEltGlob, Index iLoc, std::array<Real,Dim> &U) const
		{
            RealMatrix GradDef;
            RealMatrix CoMatGradDef;
            getGradDef(iEltGlob,iLoc,GradDef);
            ArrayAlgebra::CoMat(GradDef,CoMatGradDef);
            ArrayAlgebra::TransposeMatMlt(CoMatGradDef,U);
        }

        void getDoFCoordinateDisc(const Index iGlob, RealVector& xyz) const
        {
            const FESpace & self = static_cast<const FESpace &>(*this); 
            Index NLocDoF= self.getNumLocDoFs();
            
            Index NElem = iGlob/NLocDoF;
            Index iLoc = iGlob - NElem*NLocDoF;

            self.getDoFCoordinate(NElem, iLoc, xyz);
        }
           
        Mesh getMesh()
        {
            assert(Dim!=1);

            const FESpace & self = static_cast<const FESpace &>(*this);

            Index NumVertices;

            if constexpr(Dim == 2) NumVertices = 4;
            if constexpr(Dim == 3) NumVertices = 8;

            std::vector<std::shared_ptr<Point>> coords;
            std::vector<std::shared_ptr<Point>> local_coords(NumVertices);
            std::vector<std::shared_ptr<Element>> elements;
                
            RealVector xyz;
            for (Index v = 0; v<self.getNumVertices(); v++)
            {
                self.getVertexCoordinate(v,xyz);
                coords.push_back(std::make_shared<Point>(xyz[0],xyz[1],xyz[2],v));
            }
                
            for (Index e = 0; e<self.getNumElements(); e++)
            {
                std::vector<Index> numbering;
                self.getElementVertices(e,numbering);

                for (Index v = 0; v<NumVertices; v++)
                {
                    local_coords[v]=coords[numbering[v]];
                }

                if constexpr(Dim == 2) elements.push_back(std::make_shared<Quadrangle4>(local_coords,0));
                if constexpr(Dim == 3) elements.push_back(std::make_shared<Hexahedron8>(local_coords,0));
            }
                
            return Mesh(coords,elements);
        }

        Mesh getMeshQ2()
        {
            assert(Dim!=1);

            const FESpace & self = static_cast<const FESpace &>(*this);

            Index NumVerticesQ2;

            if constexpr(Dim == 2) NumVerticesQ2 = 9;
            if constexpr(Dim == 3) NumVerticesQ2 = 27;

            std::vector<std::shared_ptr<Point>> coords;
            std::vector<std::shared_ptr<Point>> local_coords(NumVerticesQ2);
            std::vector<std::shared_ptr<Element>> elements;
                
            RealVector xyz;
            for (Index v = 0; v<self.getNumVerticesQ2(); v++)
            {
                self.getVertexCoordinateQ2(v,xyz);
                coords.push_back(std::make_shared<Point>(xyz[0],xyz[1],xyz[2],v));
            }
                
            for (Index e = 0; e<self.getNumElements(); e++)
            {
                std::vector<Index> numbering;
                self.getElementVerticesQ2(e,numbering);

                for (Index v = 0; v<NumVerticesQ2; v++)
                {
                    local_coords[v]=coords[numbering[v]];
                }

                if constexpr(Dim == 2) elements.push_back(std::make_shared<Quadrangle9>(local_coords,0));
                if constexpr(Dim == 3) elements.push_back(std::make_shared<Hexahedron27>(local_coords,0));
            }
                
            return Mesh(coords,elements);
        }

        FESpaceTSubdomain<FESpace,Dim> getSubDomainFESpace(std::function<bool(RealVector)> isInside)
        {
            std::vector<Index> elems;
            const FESpace & self = static_cast<const FESpace &>(*this);

            Index NElem = self.getNumElements();

            for (Index e=0;e<NElem;++e)
            {
                RealVector barycenter;
                self.getBarycenter(e,barycenter);

                if (isInside(barycenter)) elems.push_back(e);
            }

            return FESpaceTSubdomain<FESpace,Dim>(self,elems);
        }
    
    };
 
}//OndoMathX


#include "FESpaceTSubdomain.hxx"



