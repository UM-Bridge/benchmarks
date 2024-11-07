 
namespace OndoMathX {
        
        template<class DerivedT,class ReferenceElement> void ElementT<DerivedT, ReferenceElement>::
        getGradDef(const RealVector& uvw , RealMatrix3x3 & graddef) const
        {
            graddef[0][0] = 0.0;
            graddef[0][1] = 0.0;
            graddef[0][2] = 0.0;
            graddef[1][0] = 0.0;
            graddef[1][1] = 0.0;
            graddef[1][2] = 0.0;
            graddef[2][0] = 0.0;
            graddef[2][1] = 0.0;
            graddef[2][2] = 0.0;
            
            std::array<Real, 3> s;
            
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            switch(getDim())
            {
                case 3 :
                    
                    for(Index i = 0; i < self.getNumNodes(); i++)
                    {
                        self.getGradShapeFunction(i, uvw, s);
                        const Point &p = self.getNode(i);
                        
                        graddef[0][0] += p[0] * s[0]; graddef[0][1] += p[0] * s[1]; graddef[0][2] += p[0] * s[2];
                        graddef[1][0] += p[1] * s[0]; graddef[1][1] += p[1] * s[1]; graddef[1][2] += p[1] * s[2];
                        graddef[2][0] += p[2] * s[0]; graddef[2][1] += p[2] * s[1]; graddef[2][2] += p[2] * s[2];
                    }
                    
                    break;
                    
                case 2 :
                    
                    for(Index i = 0; i <  self.getNumNodes(); i++)
                    {
                        self.getGradShapeFunction(i, uvw, s);
                        const Point &p = self.getNode(i);
                        
                        graddef[0][0] += p[0] * s[0]; graddef[0][1] += p[0] * s[1];
                        graddef[1][0] += p[1] * s[0]; graddef[1][1] += p[1] * s[1];
                        graddef[2][0] += p[2] * s[0]; graddef[2][1] += p[2] * s[1];
                    }
                    
                    break;
                    
                case 1:
                    
                    for(Index i = 0; i <  self.getNumNodes(); i++)
                    {
                        self.getGradShapeFunction(i, uvw, s);
                        const Point &p = self.getNode(i);
                        
                        graddef[0][0] += p[0] * s[0];
                        graddef[1][0] += p[1] * s[0];
                        graddef[2][0] += p[2] * s[0];
                    }
                    
                    break;
                    
                default: {assert(false);}
            }
        }
        
        
        template<class DerivedT,class ReferenceElement> void ElementT<DerivedT, ReferenceElement>::
        getGradDef(const RealVector& uv, RealMatrix2x2 & graddef) const
        {
            graddef[0][0] = graddef[0][1] =  0.;
            graddef[1][0] = graddef[1][1] =  0.;
            
            RealVector s;
            
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            //Check that the element is indeed embeded in R^2 (z = 0)
            for(Index i=0;i< self.getNumNodes();++i)
            {
                assert(self.getVertex(i)[2]==0.0);
            }
    
            
            switch(getDim())
            {
                    
                case 2 :
                    
                    for(Index i=0;i< self.getNumNodes();++i)
                    {
                        self.getGradShapeFunction(i, uv, s);
                        const Point &p = self.getNode(i);
                        
                        graddef[0][0] += p[0] * s[0]; graddef[0][1] += p[0] * s[1];
                        graddef[1][0] += p[1] * s[0]; graddef[1][1] += p[1] * s[1];
                    }
                    
                    break;
                case 1 :
                    
                    for(Index i=0;i< self.getNumNodes();++i)
                    {
                        self.getGradShapeFunction(i, uv, s);
                        const Point &p = self.getNode(i);
                        
                        graddef[0][0] += p[0] * s[0];
                        graddef[1][0] += p[1] * s[0];
                    }
                    
                    break;
                    
                default: {assert(false);}
            }
        }
        
        
        template<class DerivedT,class ReferenceElement> void ElementT<DerivedT, ReferenceElement>::
        getGradDef(const RealVector& u , Real & graddef) const
        {
            graddef=0;
            
            RealVector s;
            
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            //Check that the element is indeed embeded in R (y = z = 0)
            for(Index i=0;i< self.getNumNodes();++i)
            {
                assert(self.getVertex(i)[1]==0.0);
                assert(self.getVertex(i)[2]==0.0);
            }
            
            switch(getDim())
            {
                    
                case 1 :
                    
                    for(Index i=0;i< self.getNumNodes();++i)
                    {
                        self.getGradShapeFunction(i, u, s);
                        const Point &p = self.getNode(i);
                        
                        graddef += p[0] * s[0];
                    }
                    
                    break;
                    
                default: {assert(false);}
            }
        }
        
        
        template<class DerivedT,class ReferenceElement> void ElementT<DerivedT, ReferenceElement>::
        uvw2xyz(const RealVector & uvw,  RealVector & xyz) const
        {
            xyz[0] = 0.0;
            xyz[1] = 0.0;
            xyz[2] = 0.0;
            
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            for(Index i = 0; i <   self.getNumNodes(); ++i)
            {
                Real s;
                self.getShapeFunction(i, uvw , s);
                const Point &p = self.getNode(i);
                
                xyz[0] += p[0] * s;
                xyz[1] += p[1] * s;
                xyz[2] += p[2] * s;
            }
        }
        
        
        
        template<class DerivedT,class ReferenceElement> RealVector ElementT<DerivedT, ReferenceElement>::
        getBarycenter() const
        {
            RealVector xyz = {0.0, 0.0, 0.0};
            
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            Index n =  self.getNumNodes();
            for (Index i=0;i<n;i++)
            {
                const Point &p = self.getNode(i);
                
                xyz[0]+=p[0]/n;
                xyz[1]+=p[1]/n;
                xyz[2]+=p[2]/n;
            }
            
            return xyz;
        }
        
        
        template<class DerivedT,class ReferenceElement> Real ElementT<DerivedT, ReferenceElement>::
        getRadius() const
        {
            RealVector barycenter = getBarycenter();
            
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            Real dist=0;
            
            for (Index i=0;i< self.getNumNodes();i++)
            {
                Real tmp = ArrayAlgebra::Dist(self.getNode(i).coord(),barycenter);
                if (tmp>dist) dist=tmp;
            }
            
            return dist;
        }
        
        template<class DerivedT,class ReferenceElement> bool ElementT<DerivedT, ReferenceElement>::
        isInside(const RealVector &pos) const
        {
            const DerivedT & self = static_cast<const DerivedT &>(*this);
        
            RealVector barycenter = self.getBarycenter();
            Real radius = self.getRadius();
            
            //Check if we are far from the element
            if (ArrayAlgebra::Dist(barycenter,pos) > radius)
                return false;
            
            //If we are close we invert the parametrization  to check if we are inside in the reference config.
            RealVector uvw;
            bool invert_parametrization = self.xyz2uvw(pos,uvw);
            
            if (invert_parametrization==false)
            {
                return false;
                //The default behavior could be : assert(invert_parametrization==true);
                //but we assume that the parametrization is invertible by the algorithm if
                //the point is inside
            }
            
             
            
            if (ReferenceElement::isInside(uvw))
            {
                RealVector xyz;
                self.uvw2xyz(uvw,xyz);
                
                if (ArrayAlgebra::Dist(pos,xyz) > COORD_TOL)
                    return false;
                else return true;
            }
            else return false;
        }
        
        // invert the parametrisation in the least square sense ,
        // minimize the functionnal  1/2| F(uvw) - xyz |_2 with gauss-newton
        // or solving F(uvw) - xyz  (with newton algorithm)
        template<class DerivedT,class ReferenceElement> bool ElementT<DerivedT, ReferenceElement>::
        xyz2uvw(const RealVector & xyz,  RealVector & uvw) const
        {
            const DerivedT & self = static_cast<const DerivedT &>(*this);
            
            
            //Starting value
            uvw =  {0.0, 0.0, 0.0};
            
            Index iter = 1;
            Index maxiter = 100;
            Real tol = REF_COORD_TOL;
            Real diff = tol+1.0;
            
            if (getDim()==3)
            {
                //General Newton routine, solve  F(uvw) = xyz

                RealMatrix3x3 graddef;
                RealMatrix3x3 inv;
                RealVector uvw_new;
                RealVector xyz_tmp = {0.0, 0.0, 0.0};
                
                while (diff > tol && iter < maxiter)
                {
                    uvw2xyz(uvw,xyz_tmp);
                    
                    xyz_tmp[0]=xyz[0]-xyz_tmp[0];
                    xyz_tmp[1]=xyz[1]-xyz_tmp[1];
                    xyz_tmp[2]=xyz[2]-xyz_tmp[2];
                    
                    getGradDef(uvw, graddef);
                    ArrayAlgebra::Inv(graddef, inv);
                    
                    ArrayAlgebra::MatMlt(inv, xyz_tmp, uvw_new);
                    
                    uvw_new[0]+=uvw[0];
                    uvw_new[1]+=uvw[1];
                    uvw_new[2]+=uvw[2];
                    
                    diff = ArrayAlgebra::Dist(uvw_new,uvw);
                    
                    uvw[0] = uvw_new[0];
                    uvw[1] = uvw_new[1];
                    uvw[2] = uvw_new[2];
                    
                    iter++ ;
                }
              
            }
            else if (getDim()==2)
            {
                Index GeometricDim = 2;
                
                if (xyz[2]!=0.0)
                {
                    GeometricDim = 3;
                }
                else
                {
                    for(Index i = 0; i <   self.getNumNodes(); ++i)
                    {
                        const Point &p = self.getNode(i);
                    
                        if(p[2]!=0.0)
                        {
                            GeometricDim=3;
                            break;
                        }
                    }
                }
                
                if (GeometricDim==2)
                {
                    //General Newton routine, solve  F(uvw) = xyz
                    
                    RealMatrix2x2 graddef;
                    RealMatrix2x2 inv;
                    std::array<Real, 2> uv_new = {0.0, 0.0};
                    RealVector xyz_tmp = {0.0, 0.0, 0.0};
                    
                    while (diff > tol && iter < maxiter)
                    {
                        uvw2xyz(uvw,xyz_tmp);
                        
                        uv_new[0]=xyz[0]-xyz_tmp[0];
                        uv_new[1]=xyz[1]-xyz_tmp[1];
                        
                        
                        getGradDef(uvw, graddef);
                        ArrayAlgebra::Inv(graddef, inv);
                        ArrayAlgebra::MatMlt(inv,uv_new);
                        
                        uv_new[0]+=uvw[0];
                        uv_new[1]+=uvw[1];
                        
                        diff = sqrt((uv_new[0] - uvw[0])*(uv_new[0] - uvw[0])
                                    + (uv_new[1] - uvw[1])*(uv_new[1] - uvw[1]));
                        
                        uvw[0] = uv_new[0];
                        uvw[1] = uv_new[1];
                        
                        iter++ ;
                    }
                }
                else
                {
                    //gauss-newton routine, minimize 1/2 | F(uvw) - xyz |^2
                    
                    RealMatrix3x3 graddef;
                    RealMatrix2x2 hess;
                    RealMatrix2x2 inv;
                    
                    std::array<Real, 2> grad;
                    std::array<Real, 2> uv;
                     
                    RealVector xyz_tmp = {0.0, 0.0, 0.0};
                    
                    while (diff > tol && iter < maxiter)
                    {
                        uvw2xyz(uvw,xyz_tmp);
                        
                        xyz_tmp[0]=xyz[0]-xyz_tmp[0];
                        xyz_tmp[1]=xyz[1]-xyz_tmp[1];
                        xyz_tmp[2]=xyz[2]-xyz_tmp[2];
                        
                        getGradDef(uvw, graddef);
                        
                        grad[0] = xyz_tmp[0]*graddef[0][0]+xyz_tmp[1]*graddef[1][0]+xyz_tmp[2]*graddef[2][0];
                        grad[1] = xyz_tmp[0]*graddef[0][1]+xyz_tmp[1]*graddef[1][1]+xyz_tmp[2]*graddef[2][1];
                        
                        hess[0][0]=graddef[0][0]*graddef[0][0]+graddef[1][0]*graddef[1][0]+graddef[2][0]*graddef[2][0];
                        hess[0][1]=graddef[0][0]*graddef[0][1]+graddef[1][0]*graddef[1][1]+graddef[2][0]*graddef[2][1];
                        hess[1][0]=hess[0][1];
                        hess[1][1]=graddef[0][1]*graddef[0][1]+graddef[1][1]*graddef[1][1]+graddef[2][1]*graddef[2][1];
                        
                        ArrayAlgebra::Inv(hess, inv);
                        ArrayAlgebra::MatMlt(inv, grad, uv);
                        
                        uv[0] += uvw[0];
                        uv[1] += uvw[1];
                        
                        diff = sqrt((uv[0] - uvw[0])*(uv[0] - uvw[0]) + (uv[1] - uvw[1])*(uv[1] - uvw[1]));
                        
                        uvw[0] = uv[0];
                        uvw[1] = uv[1];
                        
                        iter++ ;
                    }
                    
                }
            }
            else //getDim() == 1
            {
                Index GeometricDim = 1;
                
                if (xyz[1]!=0.0)
                {
                    if (xyz[2]!=0.0) GeometricDim = 3;
                    else GeometricDim = 2;
                }
                else
                {
                    for(Index i = 0; i <   self.getNumNodes(); ++i)
                    {
                        const Point &p = self.getNode(i);
                    
                        if (p[2]!=0.0) {GeometricDim = 3; break;}
                        else {GeometricDim = 2; break;}
                    }
                }
                
                if (GeometricDim==1)
                {
                    Real graddef;
                    Real inv;
                    RealVector xyz_tmp = {0.0,0.0,0.0};
                    
                    while (diff > tol && iter < maxiter)
                    {
                        uvw2xyz(uvw,xyz_tmp);
                        
                        getGradDef(uvw, graddef);
                        inv=1.0/graddef;
                        
                        Real u = uvw[0] - inv * (xyz_tmp[0]-xyz[0]);
                        
                        diff = fabs(u - uvw[0]);
                        
                        uvw[0] = u;
                        
                        iter++;
                    }
                }
                else if (GeometricDim==2)
                {
                    //gauss-newton routine, minimize | F(uvw) - xyz |^2
                    
                    RealMatrix2x2 graddef;
                    Real hess;
                    Real grad;
                    Real u;
                    RealVector xyz_tmp = {0.0,0.0,0.0};
                    
                    
                    while (diff > tol && iter < maxiter)
                    {
                        uvw2xyz(uvw,xyz_tmp);
                        
                        xyz_tmp[0]=xyz[0]-xyz_tmp[0];
                        xyz_tmp[1]=xyz[1]-xyz_tmp[1];
                        
                        getGradDef(uvw, graddef);
                        
                        grad=xyz_tmp[0]*graddef[0][0]+xyz_tmp[1]*graddef[1][0];
                        hess=graddef[0][0]*graddef[0][0]+graddef[1][0]*graddef[1][0];
                    
                        u = uvw[0]+grad/hess;
                       
                        diff = fabs(u - uvw[0]);
                        
                        uvw[0] = u;
                      
                        iter++;
                    }
                    
                }
                else // GeometricDim==3
                {
                    //gauss-newton routine, minimize | F(uvw) - xyz |^2
                    
                    RealMatrix3x3 graddef;
                    Real hess;
                    Real grad;
                    Real u;
                    RealVector xyz_tmp;
                    
                    while (diff > tol && iter < maxiter)
                    {
                        uvw2xyz(uvw,xyz_tmp);
                        
                        xyz_tmp[0]=xyz[0]-xyz_tmp[0];
                        xyz_tmp[1]=xyz[1]-xyz_tmp[1];
                        xyz_tmp[2]=xyz[2]-xyz_tmp[2];
                        
                        getGradDef(uvw, graddef);
                        
                        grad=xyz_tmp[0]*graddef[0][0]+xyz_tmp[1]*graddef[1][0]+xyz_tmp[2]*graddef[2][0];
                        hess=graddef[0][0]*graddef[0][0]+graddef[1][0]*graddef[1][0]+graddef[2][0]*graddef[2][0];
                    
                        u = uvw[0]+grad/hess;
                       
                        diff = fabs(u - uvw[0]);
                        
                        uvw[0] = u;
                      
                        iter++;
                    }
                }
            }
          
            if (iter<maxiter) return true;
            else return false;
        }
        

} // OndoMathX

 
