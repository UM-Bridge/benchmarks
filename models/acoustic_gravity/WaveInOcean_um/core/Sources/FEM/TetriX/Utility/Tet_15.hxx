#pragma once


// -----------------------------------------------------------------------------------------//
namespace OndoMathX
{

    namespace TetriX
    {
 
        namespace Tet_15
        {
        
            const  RealVector Point[15] = 
            {{0.0,0.0,0.0},
             {1.0,0.0,0.0},
             {0.0,1.0,0.0},
             {0.0,0.0,1.0},
             
             {0.5,0.0,0.0},
             {0.5,0.5,0.0},
             {0.0,0.5,0.0},
             {0.0,0.0,0.5},
             {0.5,0.0,0.5},
             {0.0,0.5,0.5},
             
             {1.0/3.0,1.0/3.0,0.0},
             {1.0/3.0,0.0,1.0/3.0},
             {0.0,1.0/3.0,1.0/3.0},
             {1.0/3.0,1.0/3.0,1.0/3.0},

             {0.25,0.25,0.25}
             };

            const Real _wv = 17.0/5040.0;
            const Real _we = 2.0/315.0;
            const Real _wf = 9.0/560.0;
            const Real _wb = 16.0/315.0;


            const Real Weight[15] = 
            {_wv,
             _wv,
             _wv,
             _wv,
             
             _we,
             _we,
             _we,
             _we,
             _we,
             _we,
             
             _wf,
             _wf,
             _wf,
             _wf,

             _wb
             };


        inline const std::function<Real(Real &u,Real &v,Real &w,Real &z)> _ShapeFunctions[15] =
        {
            [](Real &u,Real &v,Real &w,Real &z){return z*(1 - 2*u - 2*v - 2*w + 3*u*v + 3*u*w + 3*v*w - 4*u*v*w);},
            [](Real &u,Real &v,Real &w,Real &z){return u*(1 - 2*v - 2*w - 2*z + 3*v*w + 3*v*z + 3*w*z - 4*v*w*z);},
            [](Real &u,Real &v,Real &w,Real &z){return v*(1 - 2*u - 2*w - 2*z + 3*u*w + 3*u*z + 3*w*z - 4*u*w*z);},
            [](Real &u,Real &v,Real &w,Real &z){return w*(1 - 2*u - 2*v - 2*z + 3*u*v + 3*u*z + 3*v*z - 4*u*v*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*z*(1 - 3*v - 3*w + 8*v*w);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*v*(1 - 3*w - 3*z + 8*w*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*v*z*(1 - 3*u - 3*w + 8*u*w);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*w*z*(1 - 3*u - 3*v + 8*u*v);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*w*(1 - 3*v - 3*z + 8*v*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*v*w*(1 - 3*u - 3*z + 8*u*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*z*u*v*(1-4*w);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*z*u*w*(1-4*v);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*z*v*w*(1-4*u);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*v*w*(1-4*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 256*u*v*w*z;}
        };

        inline const std::function<Real(Real &u,Real &v,Real &w,Real &z)> _DuShapeFunctions[15] =
        {
            [](Real &u,Real &v,Real &w,Real &z){return 2*u + 2*v + 2*w - 2*z - 3*u*v - 3*u*w - 3*v*w + 3*v*z + 3*w*z + 4*u*v*w - 4*v*w*z - 1;},
            [](Real &u,Real &v,Real &w,Real &z){return 2*u - 2*v - 2*w - 2*z - 3*u*v - 3*u*w + 3*v*w + 3*v*z + 3*w*z + 4*u*v*w - 4*v*w*z + 1;},
            [](Real &u,Real &v,Real &w,Real &z){return v*(4*w - 3)*(u - z);},
            [](Real &u,Real &v,Real &w,Real &z){return w*(4*v - 3)*(u - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*(u - z)*(3*v + 3*w - 8*v*w - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*v*(3*u - 3*w - 3*z - 8*u*w + 8*w*z + 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*v*(3*u + 3*w - 3*z - 8*u*w + 8*w*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*w*(3*u + 3*v - 3*z - 8*u*v + 8*v*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*w*(3*u - 3*v - 3*z - 8*u*v + 8*v*z + 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 32*v*w*(z-u);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*v*(4*w - 1)*(u - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*w*(4*v - 1)*(u - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*v*w*(4*u - 4*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*v*w*(4*u - 4*z + 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 256*v*w*(z-u);}
        };

        inline const std::function<Real(Real &u,Real &v,Real &w,Real &z)> _DvShapeFunctions[15] =
        {
            [](Real &u,Real &v,Real &w,Real &z){return 2*u + 2*v + 2*w - 2*z - 3*u*v - 3*u*w - 3*v*w + 3*u*z + 3*w*z + 4*u*v*w - 4*u*w*z - 1;},
            [](Real &u,Real &v,Real &w,Real &z){return u*(4*w - 3)*(v - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 2*v - 2*u - 2*w - 2*z - 3*u*v + 3*u*w - 3*v*w + 3*u*z + 3*w*z + 4*u*v*w - 4*u*w*z + 1;},
            [](Real &u,Real &v,Real &w,Real &z){return w*(4*u - 3)*(v - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*(3*v + 3*w - 3*z - 8*v*w + 8*w*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*(3*v - 3*w - 3*z - 8*v*w + 8*w*z + 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*(v - z)*(3*u + 3*w - 8*u*w - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*w*(3*u + 3*v - 3*z - 8*u*v + 8*u*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 32*u*w*(z - v);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*w*(1-3*u + 3*v - 3*z - 8*u*v + 8*u*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*(4*w - 1)*(v - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*w*(4*v - 4*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*w*(4*u - 1)*(v - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*w*(4*v - 4*z + 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 256*u*w*(z - v);}
        };

        inline const std::function<Real(Real &u,Real &v,Real &w,Real &z)> _DwShapeFunctions[15] =
        {
            [](Real &u,Real &v,Real &w,Real &z){return 2*u + 2*v + 2*w - 2*z - 3*u*v - 3*u*w - 3*v*w + 3*u*z + 3*v*z + 4*u*v*w - 4*u*v*z - 1;},
            [](Real &u,Real &v,Real &w,Real &z){return u*(4*v - 3)*(w - z);},
            [](Real &u,Real &v,Real &w,Real &z){return v*(4*u - 3)*(w - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 2*w - 2*v - 2*u - 2*z + 3*u*v - 3*u*w - 3*v*w + 3*u*z + 3*v*z + 4*u*v*w - 4*u*v*z + 1;},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*(3*v + 3*w - 3*z - 8*v*w + 8*v*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 32*u*v*(z - w);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*v*(3*u + 3*w - 3*z - 8*u*w + 8*u*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*(w - z)*(3*u + 3*v - 8*u*v - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*u*(1 - 3*v + 3*w - 3*z - 8*v*w + 8*v*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 4*v*(1 - 3*u + 3*w - 3*z - 8*u*w + 8*u*z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*v*(4*w -4*z - 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*(4*v - 1)*(w - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*v*(4*u - 1)*(w - z);},
            [](Real &u,Real &v,Real &w,Real &z){return 27*u*v*(4*w - 4*z + 1);},
            [](Real &u,Real &v,Real &w,Real &z){return 256*u*v*(z - w);}
        };

        Real ShapeFunction(const Index &i, const RealVector & uvw) 
        {
            Real u = uvw[0];
            Real v = uvw[1];
            Real w = uvw[2];
            Real z = 1-u-v-w;
        
            return _ShapeFunctions[i](u,v,w,z);
    
        }

        Real DShapeFunction(const Index &i, const Index &d,  const RealVector & uvw) 
        {
            Real u = uvw[0];
            Real v = uvw[1];
            Real w = uvw[2];
            Real z = 1-u-v-w;
        
            if (d == 0)
                return _DuShapeFunctions[i](u,v,w,z);
            if (d == 1)
                return _DvShapeFunctions[i](u,v,w,z);
            if (d == 2)
                return _DwShapeFunctions[i](u,v,w,z);

            return 0.0;
        }
 

  
        }

    }
        
}
// -----------------------------------------------------------------------------------------//

 