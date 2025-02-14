#pragma once

#include "Field_ScalarToVector.h"

namespace OndoMathX
{
// First part of the VTK file name (second part is in the .cpp file)
std::string NameSimu = "um";

// Adimensionned number delta = Ma/Fr
// Take H=1500, C=1500, g=10. Then delta^2 = Ma^2/Fr^2 = gH/C^2 = 0.006
Real delta2 = 0.006;  

//// Final time and output frequency
Real T_end = 10. ; //50.0; // Final time in seconds
Real OutputDeltat = 0.1; // write pressure output every outputDeltat (in seconds)

// Finite element order in the x respectively the z direction 
const Index FEOrderX = 5; 
const Index FEOrderZ = 5; 

//std::string NameSimu="Dual" +std::to_string(FEOrderX)+"_"+std::to_string(FEOrderZ)+"_"+ std::to_string(Nx)+"_"+std::to_string(Nz);
Real Lx = 43;// scaling in the x direction -> 64.5 km
Real Lz = 1; // scaling in the z direction -> 1 km
Index Nx = 43; // number of subdivisions in the x direction  
Index Nz = 10;  // number of subdivisions in the z direction  

//// Recordings at one point
const Index nPointX = 1;
const Index nPointZ = 1;
// Captors coordinates given in [0,Lx] and [0,Lz]
Real listX[nPointX] = {28}; 
Real listZ[nPointZ] = {0.9};

// parameters for the function g(t)
Real timeWidth = 25;
Real timeFreq  = 20;
Real timeStart = 0.5;

Real functionTime(Real t){
    //return Ricker(t, 2, 4);
    return DoubleSigmoid(t, timeFreq, timeStart, timeWidth);
}

// parameters for the function f(x)
// Real spaceWidth = 5; 
// Real spaceFreq = 10; 
// RealVector Center = {0.5, 0.5};
//
//Real functionSpace(const RealVector &p){
//    // Input: p = (x,z),  x \in [0,Lx], z \in [0,Lz]
//    return DoubleSigmoid(p[0], spaceFreq, Center[0] - spaceWidth/2, spaceWidth);
//}   

//  Physical parameters (density, velocity)
Real _Rho(const RealVector &p)
{   
    // Input: (x,z) in the transformed domain (= not the reference square)
    return 1.;
}

Real _c(const RealVector &p)
{   
    // Input: (x,z) in the transformed domain (= not the reference square)
    return 1.;
}

Real _Rho_cm2(const RealVector &p)
{   
    // Input: (x,z) in the transformed domain (= not the reference square)
    Real Rho = _Rho(p); 
    Real c = _c(p);
    
    return Rho / (c*c);
}


std::array<Real,2> _field_ez(const RealVector &p)
{
    return  std::array<Real,2>({0,1});
}


// Parameters for the PML
Real LxPml_L = 1.0;
Real LxPml_R = 1.0;
Real Dissip_coeff = 10;

std::array<Real, 2> ex = {1.0, 0.0};
std::array<Real, 2> ez = {0.0, 1.0};

Real rhoSurface = _Rho({0,1,0});


}