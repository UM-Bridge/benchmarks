#pragma once


namespace OndoMathX
{
double function_interpolate( std::vector<double> &xData, std::vector<double> &yData, double x)
{
    int size = xData.size();
    //std::cout<< x << " " << xData[0] << " " <<  xData[size-1] << endl;
    assert(x>= xData[0] && x <= xData[size-1]);

    float dx = xData[1] - xData[0];
    int i = 0;
    if ( x == xData[size - 1] )  {return yData[size-1];} // special case: beyond right end
    i = max((int)floor(x/dx),0);
    
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];    // points on either side (unless beyond ends)
    
    double dydx = ( yR - yL ) / ( xR - xL );                                  // gradient

    return yL + dydx * ( x - xL );                                            // linear interpolation
}

Real function_interpolate2D( std::vector<Real> &xData, std::vector<Real> &tData, std::vector<Real> &fData, Real x, Real t)
{
    // Fit with a bilinear polynomial
    // f(x,y) = a00 + a10 x + a01 t + a11 xt

    // We assume that the grids xData, tData are regular 
    // and that the data fData are ordered as follows
    // [f(x=0, t=0), f(x=1, t=0), ..., f(x=0, t=1), .... f(x=xMax, t=tMax)]

    Index size_x = xData.size();
    Index size_t = tData.size();
    Index size_f = fData.size();
    assert(x>= xData[0] && x <= xData[size_x-1]);
    assert(t>= tData[0] && t <= tData[size_t-1]);

    Real dx = xData[1] - xData[0];
    Real dt = tData[1] - tData[0]; 

    // Find the four points (x1, x2, t1, t2) surrounding the point (x,t) 
    // and the corresponding values f(x1,t1), f(x1, t2), f(x2, t1), f(x2,t2)
    Index i = max((int)floor(x/dx),0);
    Index j = max((int)floor(t/dt),0);
    Index ip1, jp1; 

    // Special cases: last point
    if ( i == size_x && j == size_t )  {return fData[size_f-1];} 

     // Special case: x=xMax 
    if ( i == size_x ) { ip1=i; } 
    else { ip1=i+1; }
    
    // Special case: t=tMax
    if ( j == size_t ) { jp1=j; } 
    else { jp1=j+1; }

    Real  x1 = xData[i];
    Real  x2 = xData[ip1];
    Real  t1 = tData[j];
    Real  t2 = tData[jp1];
    //std::cout << "x1, x2: "<< x1 << " " << x2;
    //std::cout << " ;  t1, t2: "<< t1 << " " << t2 << endl;;

    Real f11 = fData[i + j*size_x];
    Real f12 = fData[i + jp1*size_x];
    Real f21 = fData[ip1 + j*size_x];
    Real f22 = fData[ip1 + jp1*size_x];

    // Compute the coefficients of the polynomial
    Real a00 = 1/(dx*dt) * ( (x2*t2)*f11 - x2*t1*f12 - x1*t2*f21 + x1*t1*f22 ) ; 
    Real a10 = 1/(dx*dt) * (  -   t2*f11 +    t1*f12 +    t2*f21 -    t1*f22 ) ; 
    Real a01 = 1/(dx*dt) * (  -   x2*f11 +    x2*f12 +    x1*f21 -    x1*f22 ) ; 
    Real a11 = 1/(dx*dt) * (         f11 -       f12 -       f21 +       f22 ) ; 

    return a00 + a10*x + a01*t + a11*x*t; 
 
}


void read_data(const std::string& fileName, std::vector<double> &data){
    string line;
    ifstream myfile (fileName);
    if (myfile.is_open())
    {
        while ( getline (myfile,line) ) {
            data.push_back(std::stod(line));
        }
        myfile.close();
    }
    else std::cout << "Unable to open file" << std::endl;
}




template<class FESpace> Index findGlobLocIndex(const RealVector& q, FESpace & fespace)
{   
    // From a point q given in the reference square, get the closest global index 
    Index Nx = fespace.getNX();
    Index Nz = fespace.getNY();
    Index FEOrderX = fespace.getPX();
    Index FEOrderZ = fespace.getPY();
    const GaussLobattoElement GLE = fespace.getFE();

    // We assume that if the function geoTransformTopo_Inverse returns something not in [0,1]x[0,1] 
    // it means that the target point p is not in the transformed domain
    // So we return an index larger that the max index
    if(q[0]<0 or q[0]>1 or q[1]<0 or q[1]>1) return fespace.getNumDoFs() + 1  ; 

    // Get global element number
    Index iEltGlobX = Index(floor(q[0] * Nx));
    Index iEltGlobZ = Index(floor(q[1] * Nz)); 
    Index iEltGlob = iEltGlobZ * Nx + iEltGlobX; 
    
    // Get local coordinates of the point inside the element
    Real dx = Nx * q[0] - iEltGlobX ;
    Real dz = Nz * q[1] - iEltGlobZ ;

    // Get the local index: get the closest quadrature point to the local coordinates
    std::vector<Real> pointsI = GLE.getPointsI();
    for (Real& d : pointsI ) {  d = abs(d - dx); }
    std::vector<Real>::iterator resultI = std::min_element(pointsI.begin(), pointsI.end());
    Index i = std::distance(pointsI.begin(), resultI);

    std::vector<Real> pointsJ =  GLE.getPointsJ();
    for (Real& d : pointsJ ) { d = abs(d - dz); }
    std::vector<Real>::iterator resultJ = std::min_element(pointsJ.begin(), pointsJ.end());
    Index j = std::distance(pointsJ.begin(), resultJ);

    // Get the local index from i and j
    // In _GLE: index = j*NI + i with NI = FEOrderX + 1
    Index indexLoc = j*(FEOrderX+1) + i;

    // Get the global index 
    Index iGlob = fespace._Loc2Glob(iEltGlob, indexLoc);

    return iGlob;
}





////////// Various functions of space and time

Real Sigmoid (Real t, Real k, Real t0)
{
    return 1/(1+exp( -k * ( t - t0 ))) ;
}

Real DoubleSigmoid(Real t, Real k, Real t0, Real deltaT)
{
    return Sigmoid(t,k,t0) - Sigmoid(t,k, t0+deltaT);
}

Real dSigmoid (Real t, Real k, Real t0)
{
    Real den = 1/(1+exp( -k * ( t - t0 )));
    return k * exp( -k * ( t - t0 )) * (den * den) ;
}

Real dDoubleSigmoid(Real t, Real k, Real t0, Real deltaT)
{
    return dSigmoid(t,k,t0) - dSigmoid(t,k, t0+deltaT);
}


Real Ramp(Real t, Real k, Real t0)
{
    return t + 1/k * log( (1+exp(-k*(t-t0)))/(1+exp(k*t0))) ;

}

Real DoubleRamp(Real t, Real k, Real t0, Real deltaT)
{
    return Ramp(t,k,t0) - Ramp(t,k, t0+deltaT);
}

Real Sinc(Real t, Real t0, Real s0)
{
    if(abs(t-t0) < Zero_Machine) return 1;
    return sin(s0*M_PI*(t-t0)) / (s0*M_PI*(t-t0));
}

Real CInfty0_1D (Real t, Real t0, Real r0 , Real s0)
{
    Real r2 = (t-t0)*(t-t0);
    Real r0_2 = r0*r0;
    
    if (r2>=r0_2) return 0.0;
    
    Real den = (r2 - r0_2);
    Real func = exp(s0 * r0_2/den + s0);
    
    return func;
}

Real dCInfty0_1D(Real t, Real t0, Real r0, Real s0)
{
    t = t-t0;
    Real r2 = t*t;
    Real r0_2 = r0*r0;
    
    if (r2>=r0_2) return 0.0;
    
    Real den = (r2 - r0_2);
    Real num =  - 2 * s0 * t * r0_2 /(den*den);
    Real func = num*exp(s0 * r0_2/den + s0);
    
    return func;
}

Real Munk(const Real z)
{
    // Munk profile: taken from Jensen, Chap. 5.6.
    // adim version of Munk: divided by 1500.

    Real zMunk = 5000 * (1 - z);
    //std::cout<< "in Munk, z= " << zMunk <<  " " ;
    zMunk = 2*(zMunk-1300)/1300;
    double res = (1+0.00737 * (zMunk-1+exp(-zMunk)));
    return res*res ;
    
}

}