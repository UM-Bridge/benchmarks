#pragma once

namespace OndoMathX
{

    namespace LAL
    {
        inline void SolvePCG(Vector &x,
                             const Vector &b,
                             std::function<void(const Vector &, Vector &)> Mlt,
                             std::function<void(const Vector &, Vector &)> PreCondSolve,
                             Index N,
                             Index &NIter, Real &Residue,
                             Real Tol = 10e-7, Index NItermax = 200, bool verbose = false)
        {
            NIter = 0;

            Vector r;
            Vector p;
            Vector q;
            Vector z;

            Allocate(r, N);
            Allocate(p, N);
            Allocate(q, N);
            Allocate(z, N);

            Real norme_b = EuclidianNorm(b);

            Mlt(x, r);
            r = b - r;
            PreCondSolve(r, z);
            p = z;

            do
            {
                Real alpha, beta, num, den;

                Mlt(p, q);

                num = ScalarProduct(r, z);
                den = ScalarProduct(q, p);

                if (den == 0.0)
                    break;

                alpha = num / den;

                x += alpha * p;
                r -= alpha * q;

                PreCondSolve(r, z);

                beta = ScalarProduct(r, z) / num;

                p = z + beta * p;

                NIter++;

                Residue = EuclidianNorm(r);

                if (verbose)
                    std::cout << "N. Iteration: " << NIter << " \t \t Relative residue: " << Residue / norme_b << std::endl;

            } while (Residue > Tol * norme_b && NIter < NItermax);
        }
        //---------------------------------------------------------------------------------//

        // Function to Check on GradJ --> to finish
        inline void CheckGradJ(const Vector &x, std::function<Real(const Vector &)> J, std::function<void(const Vector &, Vector &)> GradJ, Index N, Real eps = 10e-4, bool verbose = true, std::function<Real(const Vector &, const Vector &)> ScalProd = [](const LAL::Vector &u, const LAL::Vector &v) -> Real
                               { return LAL::ScalarProduct(u, v); })
        {

            Vector dJx;
            Vector xm;
            Vector xp;

            Allocate(dJx, N);
            Allocate(xm, N);
            Allocate(xp, N);

            xm = x;
            xp = x;

            GradJ(x, dJx);
            Real J_value = J(x);

            for (Index i = 0; i < N; i++)
            {
                xp[i] += eps;
                xm[i] -= eps;

                Real Jp = J(xp);
                Real Jm = J(xm);

                Real dJx_app = (Jp - Jm) / (2 * eps);

                if (verbose)
                    std::cout << i << ": \t \t" << J_value << "\t \t " << dJx[i] << " \t \t " << dJx_app << "\t  \t" << fabs(dJx[i] - dJx_app) / fabs(dJx[i]) << std::endl;

                xp[i] -= x[i];
                xm[i] += x[i];
            }
        }

        //---------------------------------------------------

        inline void SolveNLPCG(Vector &x, std::function<Real(const Vector &)> J, std::function<void(const Vector &, Vector &)> GradJ, std::function<void(const Vector &, Vector &)> PreCondSolve, Index N, Index &NIter, Real &J_value, Real MaxStep = 1, Real Tol = 10e-8, Index NItermax = 200, bool verbose = false, std::function<Real(const Vector &, const Vector &)> ScalProd = [](const LAL::Vector &u, const LAL::Vector &v) -> Real
                               { return LAL::ScalarProduct(u, v); })
        {
            NIter = 0;

            Vector d;
            Vector y;
            Vector g;
            Vector g_old;
            Vector y_old;
            Vector xd;

            Allocate(d, N);
            Allocate(y, N);
            Allocate(g, N);
            Allocate(g_old, N);
            Allocate(y_old, N);
            Allocate(xd, N);

            // Init
            GradJ(x, g);
            PreCondSolve(g, y);
            d = -y;

            do
            {
                Real alpha = MaxStep;
                Real beta, num, den;
                Real Jx = J(x);
                Real m = ScalProd(g, d);

                if (m >= 0.0)
                {
                    if (verbose)
                        std::cout << "No local decrease " << std::endl;
                    break;
                }

                // Backtracking line search -- Armijo–Goldstein criterion.
                {
                    xd = x + alpha * d;
                    Real Jalpha = J(xd);

                    //  if (verbose)
                    //      std::cout << "In iteration " << NIter << " use of a direction descent step of " << alpha << std::endl;

                    while ((Jx - Jalpha) < -0.5 * alpha * m && alpha > 1e-8 * MaxStep)
                    {
                        alpha *= 0.5;
                        xd = x + alpha * d;
                        Jalpha = J(xd);

                        //   if (verbose)
                        //      std::cout << "In iteration " << NIter << " use of a direction descent step of " << alpha << std::endl;
                    }
                }

                x = x + alpha * d;

                g_old = g; // Back up of the old g
                y_old = y;

                GradJ(x, g);
                PreCondSolve(g, y);

                // Polak-Ribiere-Polyak conjugate gradient
                num = ScalProd(g, y) - ScalProd(g_old, y);
                den = ScalProd(g_old, y_old);

                if (den == 0.0)
                    break;

                beta = num / den;

                d = -y + beta * d;

                NIter++;

                J_value = J(x);

                if (verbose)
                    std::cout << "N. Iteration: " << NIter << " \t \t Value of the functionnal: " << J_value << std::endl;

            } while (NIter < NItermax);
        }
        //---------------------------------------------------------------------------------//

        inline void SolveNLPCG(Vector &x, std::function<Real(const Vector &)> J, std::function<void(const Vector &, Vector &)> GradJ, std::function<void(const Vector &, Vector &)> PreCondSolve, Index &Converged, Index N, Index &NIter, Real &J_value, Real MaxStep = 1, Real Tol = 10e-8, Index NItermax = 200, bool verbose = false, std::function<Real(const Vector &, const Vector &)> ScalProd = [](const LAL::Vector &u, const LAL::Vector &v) -> Real
                               { return LAL::ScalarProduct(u, v); })
        {
            NIter = 0;

            Vector d;
            Vector y;
            Vector g;
            Vector g_old;
            Vector y_old;
            Vector xd;

            Allocate(d, N);
            Allocate(y, N);
            Allocate(g, N);
            Allocate(g_old, N);
            Allocate(y_old, N);
            Allocate(xd, N);

            // Init
            GradJ(x, g);
            PreCondSolve(g, y);
            d = -y;

            do
            {
                Real alpha = MaxStep;
                Real beta, num, den;
                Real Jx = J(x);
                Real m = ScalProd(g, d);

                if (m >= 0.0)
                {
                    if (verbose)
                        std::cout << "No local decrease " << std::endl;
                    Converged = 1;
                    break;
                }

                // Backtracking line search -- Armijo–Goldstein criterion.
                {
                    xd = x + alpha * d;
                    Real Jalpha = J(xd);

                    //  if (verbose)
                    //      std::cout << "In iteration " << NIter << " use of a direction descent step of " << alpha << std::endl;

                    while ((Jx - Jalpha) < -0.5 * alpha * m && alpha > 1e-8 * MaxStep)
                    {
                        alpha *= 0.5;
                        xd = x + alpha * d;
                        Jalpha = J(xd);

                        //   if (verbose)
                        //      std::cout << "In iteration " << NIter << " use of a direction descent step of " << alpha << std::endl;
                    }
                }

                x = x + alpha * d;

                g_old = g; // Back up of the old g
                y_old = y;

                GradJ(x, g);
                PreCondSolve(g, y);

                // Polak-Ribiere-Polyak conjugate gradient
                num = ScalProd(g, y) - ScalProd(g_old, y);
                den = ScalProd(g_old, y_old);

                if (den == 0.0)
                    break;

                beta = num / den;

                d = -y + beta * d;

                NIter++;

                J_value = J(x);

                if (verbose)
                    std::cout << "N. Iteration: " << NIter << " \t \t Value of the functionnal: " << J_value << std::endl;

            } while (NIter < NItermax);
        }
        //---------------------------------------------------------------------------------//

    } // LAL

} // OndoMathX
