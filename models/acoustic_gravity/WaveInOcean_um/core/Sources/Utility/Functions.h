#pragma once

namespace OndoMathX
{

    // Smooth function with compact support
    template <Index Dim = 3>
    inline Real CInfty0(const RealVector &xyz, const RealVector &Center, Real r0, Real s0)
    {
        Real r_2 = 0;
        Real r0_2 = r0 * r0;

        for (Index d = 0; d < Dim; ++d)
            r_2 += (xyz[d] - Center[d]) * (xyz[d] - Center[d]);

        if (r_2 >= r0_2)
            return 0.0;

        Real den = (r_2 - r0_2);
        Real func = exp(s0 * r0_2 / den + s0);

        return func;
    }

    // Gradient of the smooth function with compact support
    template <Index Dim = 3>
    inline std::array<Real, Dim> dCInfty0(const RealVector &xyz, const RealVector &Center, Real r0, Real s0)
    {
        Real r = 0;
        Real r0_2 = r0 * r0;

        std::array<Real, Dim> grad;

        for (Index d = 0; d < Dim; ++d)
            grad[d] = 0.0;

        for (Index d = 0; d < Dim; ++d)
            r += (xyz[d] - Center[d]) * (xyz[d] - Center[d]);

        if (r >= r0_2)
            return grad;

        Real den = (r - r0_2);
        Real func = exp(s0 * r0 / den + s0);
        Real d_func = -s0 * r0 * func / (den * den);

        for (Index d = 0; d < Dim; ++d)
        {
            Real dr = 2 * (xyz[d] - Center[d]);
            grad[d] = d_func * dr;
        }

        return grad;
    }

    // Ricker function
    inline Real Ricker(Real t, Real t0, Real s0)
    {
        t = t - t0;
        return (4 * s0 * s0 * t * t - 2 * s0) * exp(-t * t * s0);
    }

    // Gaussian function
    inline Real Gaussian(Real t, Real t0, Real s0)
    {
        t = t - t0;
        return exp(-t * t * s0);
    }

    // Gaussian function
    inline Real dGaussian(Real t, Real t0, Real s0)
    {
        t = t - t0;
        return -2 * s0 * t * exp(-t * t * s0);
    }

    // Solution u of the wave equation in Rˆ2 with
    // zero source term, u(t=0)=0 and u'(t=0) = exp(-alpha |x|^2)
    inline Real wave_solution(const RealVector &P, Real t, Real alpha)
    {
        constexpr Index N_t = 100;
        constexpr Index N_I = 100;

        Real Rx = ArrayAlgebra::Norm(P);

        auto _bess_func = [=](Real Ry)
        {
            // Gauss quadrature with 3 points
            std::array<Real, 3> quad_pos;
            std::array<Real, 3> quad_weight;

            Real sum = 0.0;
            Real theta = 0;
            Real delta_theta = 2.0 * M_PI / N_I;

            quad_weight[0] = delta_theta * (5. / 9.) * 0.5;
            quad_weight[1] = delta_theta * (8. / 9.) * 0.5;
            quad_weight[2] = delta_theta * (5. / 9.) * 0.5;

            quad_pos[0] = delta_theta * (-sqrt(3. / 5.) + 1.0) / 2.0;
            quad_pos[1] = delta_theta * 0.5;
            quad_pos[2] = delta_theta * (sqrt(3. / 5.) + 1.0) / 2.0;

            Real tmp_1 = Rx * Rx + Ry * Ry;
            Real tmp_2 = 2 * Rx * Ry;

            for (Index I = 0; I < N_I; ++I)
            {
                sum += quad_weight[0] * exp(-alpha * (tmp_1 - tmp_2 * cos(theta + quad_pos[0])));
                sum += quad_weight[1] * exp(-alpha * (tmp_1 - tmp_2 * cos(theta + quad_pos[1])));
                sum += quad_weight[2] * exp(-alpha * (tmp_1 - tmp_2 * cos(theta + quad_pos[2])));

                theta += delta_theta;
            }

            return Ry * sum / (2 * M_PI * sqrt(t + Ry));
        };

        Real Ry;
        Real delta_t = t / N_t;
        Real sum = 0.0;

        // Gauss quadrature with 3 points
        std::array<Real, 3> quad_pos;
        std::array<Real, 3> quad_weight;

        quad_weight[0] = delta_t * (5. / 9.) * 0.5;
        quad_weight[1] = delta_t * (8. / 9.) * 0.5;
        quad_weight[2] = delta_t * (5. / 9.) * 0.5;

        quad_pos[0] = delta_t * (-sqrt(3. / 5.) + 1.0) / 2.0;
        quad_pos[1] = delta_t * 0.5;
        quad_pos[2] = delta_t * (sqrt(3. / 5.) + 1.0) / 2.0;

        Real tmp = _bess_func(t);

        for (Index n = 0; n < N_t; ++n)
        {
            Ry = n * delta_t;

            sum += quad_weight[0] * (_bess_func(t - Ry - quad_pos[0]) - tmp) / sqrt(Ry + quad_pos[0]);

            sum += quad_weight[1] * (_bess_func(t - Ry - quad_pos[1]) - tmp) / sqrt(Ry + quad_pos[1]);

            sum += quad_weight[2] * (_bess_func(t - Ry - quad_pos[2]) - tmp) / sqrt(Ry + quad_pos[2]);
        }

        return sum + 2.0 * sqrt(t) * tmp;
    }

} // OndoMathX
