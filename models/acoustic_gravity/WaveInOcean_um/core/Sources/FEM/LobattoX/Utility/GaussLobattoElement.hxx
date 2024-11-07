#pragma once

// File Regrouping various functions used to define spectral finite elements.


// -----------------------------------------------------------------------------------------//
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>
#include <algorithm>
// -----------------------------------------------------------------------------------------//


// -----------------------------------------------------------------------------------------//
namespace OndoMathX {
    
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*! \brief Intermediary function used for the computation of Gauss-Lobatto points
 */
inline void ComputeFourierCoef(Index n, Real & cz, std::vector<Real> & cp,
                         std::vector<Real> & dcp, std::vector<Real> & ddcp)
{
    //
    //     computes the fourier coefficients of the legendre
    //     polynomial p_n0 and its derivative.
    //     n is the degree and n/2 or (n+1)/2
    //     coefficients are returned in cp depending on whether
    //     n is even or odd. The same number of coefficients
    //     are returned in dcp. For n even the constant
    //     coefficient is returned in cz.
    
    Real t1 = -1.0;
    Real t2 = n + 1.0;
    Real t3 = 0.0;
    Real t4 = 2*n + 1.0;

    Index ncp = (n + 1) / 2;

    if (n % 2 == 0)
    {
        cp[ncp] = 1.0;
        for (Index j = ncp; j >= 2; j--)
        {
            t1 += 2.0;
            t2 += -1.0;
            t3 += 1.0;
            t4 += -2.0;
            cp[j - 1] = (t1*t2) / (t3*t4)*cp[j];
        }
        
        t1 = t1 + 2.0;
        t2 = t2 - 1.0;
        t3 = t3 + 1.0;
        t4 = t4 - 2.0;
        cz = (t1*t2) / (t3*t4)*cp[1];
        
        for (Index j = 1; j <= ncp; j++)
        {
            dcp[j] = (j + j)*cp[j];
            ddcp[j] = (j + j)*dcp[j];
        }
    }
    else
    {
        cp[ncp] = 1.0;
        for (Index j = ncp - 1; j >= 1; j--)
        {
            t1 = t1 + 2.0;
            t2 = t2 - 1.0;
            t3 = t3 + 1.0;
            t4 = t4 - 2.0;
            cp[j] = (t1*t2) / (t3*t4)*cp[j + 1];
        }
        for (Index j = 1; j <= ncp; j++)
        {
            dcp[j] = (j + j - 1)*cp[j];
            ddcp[j] = (j + j - 1)*dcp[j];
        }
    }
    
    return;
}





///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

/*! \brief Intermediary function used for the computation of Gauss points
 */
inline void ComputeLegendrePolAndDerivative(Index n, Real theta, const Real & cz,
                                     const std::vector<Real> & cp,
                                     const std::vector<Real> & dcp,
                                     const std::vector<Real>& ddcp,
                                     Real & pb, Real & dpb, Real & ddpb)
{
    //     computes pn(theta) and its derivative dpb(theta) with
    //     respect to theta
    //
    
    
    Real cdt = cos(2 * theta);
    Real sdt = sin(2 * theta);
    Real cth = 0.0;
    Real sth = 0.0;
    Real chh = 0.0;
    
    if (n % 2 == 0)
    {
        
        // n even
        Index kdo = n / 2;
        pb = cz / 2;
        dpb = Real(0.0);
        ddpb = Real(0.0);
        if (n > 0)
        {
            cth = cdt;
            sth = sdt;
            for (Index k = 1; k <= kdo; k++)
            {
                //      pb = pb+cp(k)*cos(2*k*theta)
                pb = pb + cp[k] * cth;
                //      dpb = dpb-(k+k)*cp(k)*sin(2*k*theta)
                dpb = dpb - dcp[k] * sth;
                ddpb = ddpb - ddcp[k] * cth;
                chh = cdt*cth - sdt*sth;
                sth = sdt*cth + cdt*sth;
                cth = chh;
            }
        }
    }
    else
    {
        
        //  n odd
        Index kdo = (n + 1) / 2;
        pb = 0.0;
        dpb = 0.0;
        ddpb = 0.0;
        cth = cos(theta);
        sth = sin(theta);
        
        for (Index k = 1; k <= kdo; k++)
        {
            //     pb = pb+cp(k)*cos((2*k-1)*theta)
            pb = pb + cp[k] * cth;
            // dpb = dpb-(k+k-1)*cp(k)*sin((2*k-1)*theta)
            dpb = dpb - dcp[k] * sth;
            ddpb = ddpb - ddcp[k] * cth;
            chh = cdt*cth - sdt*sth;
            sth = sdt*cth + cdt*sth;
            cth = chh;
        }
    }
    
    return;
}


    
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    /* Compute points and weight of a gauss quadrature forumal from the number of wanted points.
     */
    
    inline void _ComputeGaussFormulas(
                              Index nb_points_gauss,
                              std::vector<Real> & points,
                              std::vector<Real> & weights,
                              bool gauss_lobatto)
    {
        // For details about that function, take a look at
        //
        //     P. N. Swarztrauber, Computing the points and weights for
        //     Gauss-Legendre quadrature, SIAM J. Sci. Comput.,
        //     24(2002) pp. 945-954.
        //
        //
        
        const Real eps = numeric_limits<Real>::epsilon();
        
        const Real pi = M_PI;
        const Real pis2 = 0.5*pi;
        
        Index nb_points_quadrature = nb_points_gauss;
        
        // points(i) will be equal to cos(theta(i))
        vector<Real> theta;
        
        // for Gauss-Lobatto, we compute internal nodes
        if (gauss_lobatto) nb_points_quadrature--;
        
        theta.resize(nb_points_gauss);
        points.resize(nb_points_gauss);
        weights.resize(nb_points_gauss);
        
        // coefficient of Fourier transform of Legendre polynom and derivative
        // we need the second derivative for Gauss-Lobatto points
        vector<Real> coef(nb_points_quadrature, 0.0);
        vector<Real> deriv_coef(nb_points_quadrature, 0.0);
        vector<Real> second_deriv_coef(nb_points_quadrature, 0.0);
        
        // exact expressions are used for order 1, 2, 3
        if (gauss_lobatto)
        {
            if (nb_points_quadrature == 1)
            {
                // [0,1] with weights 1/2, 1/2
                points[0] = 0.0;
                points[1] = 1.0;
                
                weights[0] = 0.5;
                weights[1] = 0.5;
                return;
            }
            else if (nb_points_quadrature == 2)
            {
                // [0, 1/2, 1] with weights 1/6, 2/3, 1/6
                points[0] = 0.0;
                points[1] = 0.5;
                points[2] = 1.0;
                
                weights[0] = 1.0 / 6.0;
                weights[1] = 2.0 / 3.0;
                weights[2] = weights[0];
                
                return;
            }
            else if (nb_points_quadrature == 3)
            {
                points[0] = 0.0;
                points[1] = (5.0 - sqrt(5.0)) / 10.0;
                points[2] = 1.0 - points[1];
                points[3] = 1.0;
                
                weights[0] = 1.0 / 12.0;
                weights[1] = 5.0 / 12.0;
                weights[2] = weights[1];
                weights[3] = weights[0];
            }
        }
        else
        {
            if (nb_points_quadrature == 1)
            {
                // 1/2 with weight 1
                points[0] = 0.5;
                weights[0] = 1.0;
                return;
            }
            else if (nb_points_quadrature == 2)
            {
                // [(3 - sqrt(3))/6, (3 + sqrt(3))/6] with weights 1/2, 1/2
                points[0] = (3.0 - sqrt(3.0)) / 6.0;
                points[1] = 1.0 - points[0];
                
                weights[0] = 0.5;
                weights[1] = 0.5;
                
                return;
            }
            else if (nb_points_quadrature == 3)
            {
                points[0] = 0.5 - sqrt(15.0) / 10.0;
                points[1] = 0.5;
                points[2] = 1.0 - points[0];
                
                weights[0] = 5.0 / 18.0;
                weights[1] = 4.0 / 9.0;
                weights[2] = weights[0];
                
                return;
            }
        }
        
        // for order greater than 3, numerical computation
        Index parity_nb_points = nb_points_quadrature % 2;
        Index half_nb_points = nb_points_gauss / 2;
        Index nhalf = (nb_points_quadrature + 1) / 2;
        
        Real zero = 0.0;
        Real cz = zero;
        ComputeFourierCoef(nb_points_quadrature, cz, coef, deriv_coef, second_deriv_coef);
        
        Real dtheta = zero;
        if (gauss_lobatto && (parity_nb_points == 1))
            dtheta = pis2 / (nhalf - 1);
        else
            dtheta = pis2 / nhalf;
        
        Real dthalf = dtheta / 2.0;
        Real cmax = dtheta / 5.0;
        // the zeros of Legendre polynom
        Real zprev = zero, zlast = zero, zhold = zero;
        Real pb = zero, dpb = zero, ddpb = zero, dcor = zero, sgnd = zero;
        Index nix;
        //
        //     estimate first point next to theta = pi/2
        //
        if (parity_nb_points != 0)
        {
            // nb_points_quadrature = 2 nhalf-1
            // if odd the first zero is at the middle pi/2
            // and the following pi/2-pi/(2*nhalf)
            if (gauss_lobatto)
            {
                zero = pis2 - dthalf;
                zprev = pis2 + dthalf;
                nix = nhalf - 1;
            }
            else
            {
                zero = pis2 - dtheta;
                zprev = pis2;
                nix = nhalf - 1; // index of the point
            }
        }
        else
        {
            // if even, no zero on the middle, the first zero is on pi/2-pi/(4*nhalf)
            if (gauss_lobatto)
            {
                zero = pis2 - dtheta;
                zprev = pis2;
                nix = nhalf - 1; // index of the point
            }
            else
            {
                zero = pis2 - dthalf;
                zprev = pis2 + dthalf;
                nix = nhalf;
            }
        }
        
        bool each_point_not_computed = true;
        while (each_point_not_computed)
        {
            int nb_iter = 0;
            bool test_convergence_newton = true;
            Real residu(1.0), residu_prec(1.0);
            while ((test_convergence_newton) && (nb_iter < 100))
            {
                nb_iter++;
                zlast = zero;
                //
                //     newton iterations
                //
                ComputeLegendrePolAndDerivative(nb_points_quadrature, zero, cz, coef, deriv_coef, second_deriv_coef, pb, dpb, ddpb);
                
                if (gauss_lobatto)
                    dcor = dpb / ddpb;
                else
                    dcor = pb / dpb;
                
                sgnd = 1.0;
                if (dcor != 0.0) sgnd = dcor / fabs(dcor);
                
                // we don't move the point further than 0.2*delta_theta
                Real tmp = fabs(dcor);
                dcor = sgnd*min(tmp, cmax);
                zero = zero - dcor;
                residu_prec = residu;
                residu = fabs(zero - zlast);
                // we check if the stopping criteria are reached
                if ((fabs(zero - zlast) <= eps*fabs(zero)) || ((nb_iter > 5) && (residu_prec < residu))) test_convergence_newton = false;
            }
            
            theta[nix - 1] = zero;
            zhold = zero;
            //      weights(nix) = (nb_points_quadrature+nb_points_quadrature+1)/(dpb*dpb)
            //
            //     yakimiw's formula permits using old pb and dpb
            //
            if (gauss_lobatto)
            {
                Real tmp = pb - dcor*dpb + 0.5*dcor*dcor*ddpb;
                tmp *= tmp;
                weights[nix - 1] = 1.0 / tmp;
            }
            else
            {
                Real tmp = dpb + pb*cos(zlast) / sin(zlast); tmp *= tmp;
                weights[nix - 1] = Real(nb_points_quadrature + nb_points_quadrature + 1) / tmp;
            }
            nix--;
            
            if (nix == 0)
                each_point_not_computed = false;
            else if (nix <= (nhalf - 1))
                zero += zero - zprev;
            
            zprev = zhold;
        }
        
        //
        //     extend points and weights via symmetries
        //
        if ((!gauss_lobatto) && (parity_nb_points != 0))
        {
            theta[nhalf - 1] = pis2;
            ComputeLegendrePolAndDerivative(nb_points_quadrature, pis2, cz, coef, deriv_coef, second_deriv_coef, pb, dpb, ddpb);
            
            weights[nhalf - 1] = (nb_points_quadrature + nb_points_quadrature + 1) / (dpb*dpb);
        }
        
        if ((gauss_lobatto) && (parity_nb_points == 0))
        {
            theta[nhalf - 1] = pis2;
            ComputeLegendrePolAndDerivative(nb_points_quadrature, pis2, cz, coef, deriv_coef, second_deriv_coef, pb, dpb, ddpb);
            weights[nhalf - 1] = 1.0 / (pb*pb);
        }
        
        // DISP(nhalf); DISP(theta);DISP(weights);
        if (gauss_lobatto)
        {
            for (int i = (int(nhalf) - 1); i >= 0; i--)
            {
                theta[(Index)(i + 1)] = theta[(Index)(i)];
                weights[(Index)(i + 1)] = weights[(Index)(i)];
            }
            
            
            theta[0] = 0.0;
            ComputeLegendrePolAndDerivative(nb_points_quadrature, theta[0], cz, coef, deriv_coef, second_deriv_coef, pb, dpb, ddpb);
            
            // DISP(cz); DISP(pb);
            weights[0] = 1.0 / (pb*pb);
        }
        
        // DISP(theta);
        for (Index i = 0; i < half_nb_points; i++)
        {
            weights[nb_points_gauss - i - 1] = weights[i];
            theta[nb_points_gauss - i - 1] = pi - theta[i];
        }
        
        // DISP(theta); DISP(weights);
        Real sum = 0.0;
        for (Index i = 0; i < nb_points_gauss; i++)
        {
            sum += weights[i];
        }
        
        for (Index i = 0; i < nb_points_gauss; i++)
        {
            weights[i] = weights[i] / sum;
            points[i] = 0.5 + 0.5*cos(theta[i]);
        }
        
        //reverse the arrays
        for (Index i = 0; i < weights.size() / 2; i++)
        {
            Real tmp = weights[i];
            weights[i] = weights[weights.size() - i - 1];
            weights[weights.size() - i - 1] = tmp;
        }
        
        for (Index i = 0; i < points.size() / 2; i++)
        {
            Real tmp = points[i];
            points[i] = points[points.size() - i - 1];
            points[points.size() - i - 1] = tmp;
        }
    }
    
    
    
    
    
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    /*! \brief 1D Interpolation at 'p' of the basis function 'phi'. The lagrangian basis are built upon
     the points 'pos'. 'phi' is between 0 and n-1. No identical position are allowed.
     */
    template <class Vector> inline Real _Interpolation(const Vector &pos, Index n, Real p, Index phi)
    {
        Real s = 1.0;
        
        for (Index r = 0; r < n; ++r) if (r != phi) s *= (p - pos[r]) / (pos[phi] - pos[r]);
        return s;
    }
    
   
    template <class Vector> inline Real _InterpolationWeight(const Vector &pos, Index n, Index phi)
    {
        Real s = 1.0;

        for (Index r = 0; r < n; ++r) if (r != phi) s /= (pos[phi] - pos[r]);
         
        return s;
    }
    
    

    
    
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
    
    /*! \brief 1D Interpolation at 'p' of the derivative of a basis function 'phi'.
     The lagrangian basis are built upon the points 'pos'.
     'phi' is between 0 and n-1. No identical position are allowed.
     */
    template <class Vector> inline Real _DerivativeInterpolation(Vector &pos, Index n, Real p, Index phi)
    {
        Real c;
        Real s = 0.0;
        Real d = 1.0;
        
        for (Index r = 0; r < n; ++r) if (r != phi) d *= (pos[phi] - pos[r]);
        
        for (Index r1 = 0; r1 < n; ++r1){
            if (r1 != phi){
                c = 1.0;
                for (Index r2 = 0; r2 < n; ++r2) if (r2 != r1 && r2 != phi) c *= (p - pos[r2]);
                s += c;
            }
        }
        return s / d;
    }
    
    
    
 

} // OndoMathX


 
