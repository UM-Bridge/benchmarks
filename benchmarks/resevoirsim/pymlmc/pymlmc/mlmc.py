import numpy
import numpy.linalg

class WeakConvergenceFailure(Exception):
    pass

def mlmc(Lmin, Lmax, N0, eps, mlmc_fn, alpha_0, beta_0, gamma_0, *args, **kwargs):
    """
    Multilevel Monte Carlo estimation.

    (P, Nl, Cl) = mlmc(...)

    Inputs:
      N0:   initial number of samples    >  0
      eps:  desired accuracy (rms error) >  0
      Lmin: minimum level of refinement  >= 2
      Lmax: maximum level of refinement  >= Lmin

      mlmc_fn: the user low-level routine for level l estimator. Its interface is

        (sums, cost) = mlmc_fn(l, N, *args, **kwargs)

        Inputs:  l: level
                 N: number of paths
                 *args, **kwargs: optional additional user variables

        Outputs: sums[0]: sum(Y)
                 sums[1]: sum(Y**2)
                    where Y are iid samples with expected value
                        E[P_0]            on level 0
                        E[P_l - P_{l-1}]  on level l > 0
                 cost: cost of N samples

      alpha ->  weak error is  O(2^{-alpha*l})
      beta  ->  variance is    O(2^{-beta*l})
      gamma ->  sample cost is O(2^{ gamma*l})

      If alpha, beta are not positive then they will be estimated.

      *args, **kwargs = optional additional user variables to be passed to mlmc_fn

    Outputs:
      P:  value
      Nl: number of samples at each level
      Cl: cost of samples at each level
    """

    # Check arguments

    if Lmin < 2:
        raise ValueError("Need Lmin >= 2")
    if Lmax < Lmin:
        raise ValueError("Need Lmax >= Lmin")
    if N0 <= 0 or eps <= 0:
        raise ValueError("Need N0 > 0, eps > 0")

    # Initialisation

    alpha = max(0, alpha_0)
    beta  = max(0, beta_0)
    gamma = max(0, gamma_0)

    theta = 0.25

    L = Lmin

    Nl   = numpy.zeros(L)
    suml = numpy.zeros((2, L))
    costl = numpy.zeros(L)
    dNl  = N0*numpy.ones(L)

    while sum(dNl) > 0:

        # update sample sums

        for l in range(0, L):
            if dNl[l] > 0:
                (sums, cost) = mlmc_fn(l, int(dNl[l]), *args, **kwargs)
                Nl[l]        = Nl[l] + dNl[l]
                suml[0, l]   = suml[0, l] + sums[0]
                suml[1, l]   = suml[1, l] + sums[1]
                costl[l]     = costl[l] + cost

        # compute absolute average, variance and cost

        ml = numpy.abs(       suml[0, :]/Nl)
        Vl = numpy.maximum(0, suml[1, :]/Nl - ml**2)
        Cl = costl/Nl

        # fix to cope with possible zero values for ml and Vl
        # (can happen in some applications when there are few samples)

        for l in range(2, L):
            ml[l-1] = max(ml[l-1], 0.5*ml[l-2]/2**alpha)
            Vl[l-1] = max(Vl[l-1], 0.5*Vl[l-2]/2**beta)

        # use linear regression to estimate alpha, beta, gamma if not given
        if alpha_0 <= 0:
            A = numpy.ones((L, 2)); A[:, 0] = range(1, L+1)
            x = numpy.linalg.lstsq(A, numpy.log2(ml[1:]))[0]
            alpha = max(0.5, -x[0])

        if beta_0 <= 0:
            A = numpy.ones((L, 2)); A[:, 0] = range(1, L+1)
            x = numpy.linalg.lstsq(A, numpy.log2(Vl[1:]))[0]
            beta = max(0.5, -x[0])

        if gamma_0 <= 0:
            A = numpy.ones((L, 2)); A[:, 0] = range(1, L+1)
            x = numpy.linalg.lstsq(A, numpy.log2(Cl[1:]))[0]
            gamma = max(0.5, x[0])

        # set optimal number of additional samples

        Ns = numpy.ceil( numpy.sqrt(Vl/Cl) * sum(numpy.sqrt(Vl*Cl)) / ((1-theta)*eps**2) )
        dNl = numpy.maximum(0, Ns-Nl)

        # if (almost) converged, estimate remaining error and decide
        # whether a new level is required

        if sum(dNl > 0.01*Nl) == 0:
            # 23/3/18 this change copes with cases with erratic ml values
            rang = list(range(min(3, L)))
            rem = ( numpy.amax(ml[[L - 1 - x for x in rang]] / 2.0**(numpy.array(rang)*alpha))
                    / (2.0**alpha - 1.0) )
            # rem = ml[L] / (2.0**alpha - 1.0)

            if rem > numpy.sqrt(theta)*eps:
                if L == Lmax:
                    raise WeakConvergenceFailure("Failed to achieve weak convergence")
                else:
                    L = L + 1
                    Vl = numpy.append(Vl, Vl[-1] / 2.0**beta)
                    Nl = numpy.append(Nl, 0.0)
                    suml = numpy.column_stack([suml, [0, 0]])
                    Cl = numpy.append(Cl, Cl[-1]*2**gamma)
                    costl = numpy.append(costl, 0.0)

                    Ns = numpy.ceil( numpy.sqrt(Vl/Cl) * sum(numpy.sqrt(Vl*Cl))
                            / ((1-theta)*eps**2) )
                    dNl = numpy.maximum(0, Ns-Nl)

    # finally, evaluate the multilevel estimator
    P = sum(suml[0,:]/Nl)

    return (P, Nl, Cl)
