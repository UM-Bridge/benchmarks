import numpy as np
import multiprocessing as mp
import scipy.stats as stats
import time
import os
import random

class ParticleAprx:
    def __init__(self, prior, M, correction_steps, ess_ratio, min_increase = 0.001, init=True):
        """
        CLASS OBJECT: ParticleAprx = SMC approximation of the parameter distributions


        >>> INPUT to initialize
        prior               prior object, containing the prior information, for sampling and its logpdf in RWMH
        M                   number of particles
        correction_steps    number of applications of the MCMC kernel (per mutation step)
        ess_ratio           ratio of M, at which particles should be resampled
        min_increase        (see TEMPERING MODIFICATIONS)


        >>> ATTRIBUTES
        prior               (see INPUT)
        M                   (see INPUT)
        d                   dimension of parameter space

        correction_steps    (see INPUT)
        ess_ratio           (see INPUT)

        acceptance_rate     avg. acceptance rate of the particles
        scaling             scaling factor for the MCMC kernel's variance
        ess                 effective sample size
        
        particles           array of particles, shape: [1:d,1:M]
        weights             associated weights
        means               array of marginal means
        variances           array of marginal variances
        MAPs                array of marginal MAPs (maximum a posteriori estimates)
        logLL_particles     array, in which the current loglikelihoods of the particles are stored 

        >> for model comparison:
        BIC_evol            array, in which the BIC values per SMC step should be stored
        BIC                 BIC value of the current approximation
        evidence            array, in which the model evidence values per SMC step should be stored
        Z                   model evidence of the current approximation
        
        >> for TEMPERING :
        min_increase        minimum exponent increase for tempering
        former_T            former tempering exponent
        current_T           current tempering exponent
        temperingfinished   boolean which is 1 when the tempering is finished (i.e. when the tempering exponent is 1)
        diff_tar_impLL      array of differences of loglikelihoods (for updating each filtering step), tempering modification
        
        > attributes to save/store information about the tempering steps:
        ess_temp
        tempering_steps
        list_essTemp
        list_exponents
        list_unique
        list_MCMC
        list_highest 
        list_scaling 
        temperingfinished 
        reweightingfinished 
        lastMCMCdone 
        total_accepted 
        total_accepted_tmp 
        
        """
        self.prior = prior
        self.M = M;  self.d = prior.d

        ### NEW: v3 MCMC
        self.correction_steps = correction_steps  # now does not matter as it's chosen adaptively
        self.corrsteps_fct = 1.  # adaptivity modification
        self.ess_ratio = ess_ratio

        ### create an initial Monte Carlo approximation of the prior
        if init:
            self.step = 0
            self.acceptance_rate = 0.25 # should be 0.5 for parameter dimension d<=2
                            
            self.scaling = 1.
            self.ess = None

            ### NEW: v1 Tempering
            self.min_increase = min_increase
            self.former_T = 0.
            self.current_T = 1.
            self.temperingfinished = False

            ### NEW: v3 MCMC
            self.nsteps_max = 64 # upper bound on n. MCMC steps per SMC iteration
            self.nsteps_min = 4  # lower bound

            self.BIC_evol = None; self.BIC = 0.
            self.evidence = None; self.Z = 0. # in scale log10(Z)

            self.particles = self.prior.rvs(self.M) # initial sample from prior
            self.weights = np.full(self.M, 1. / self.M) # uniform weights
            self.means = np.array([np.dot(self.particles[i], self.weights) for i in range(self.d)]) # sample means
            self.variances = np.array([np.dot((self.particles[i]-self.means[i])**2, self.weights) for i in range(self.d)]) # sample variances
            self.MAPs = self.means
            self.logLL_particles = np.zeros(self.M) # Was vlogPDF and I changed, it was probably a bug

            ### NEW: Tempering
            self.diff_tar_impLL = np.zeros(self.M) # initialization will be overwritten, so it does not really matter, tempering modification
            self.ess_temp = 0.
            self.tempering_steps = -1
            self.list_essTemp = []
            self.list_exponents = []
            self.list_unique = []
            self.list_MCMC = []
            self.list_highest = []
            self.list_scaling = []
            self.temperingfinished = True
            self.reweightingfinished = False
            self.lastMCMCdone = 0
            self.total_accepted = 0
    
    def save(self, filename):
        """
        saves the particle approximation in a file for later use
        
        >> INPUT        filename:  name (string) to use for the file
        >> no OUTPUT    (but generates .npz file)
        """
        np.savez(filename, prior=self.prior, M=self.M, d=self.d, correction_steps=self.correction_steps,
                 ess_ratio=self.ess_ratio, acceptance_rate=self.acceptance_rate, scaling=self.scaling, ess=self.ess,
                 particles=self.particles, weights=self.weights, means=self.means, variances=self.variances,
                 evidence=self.evidence, Z=self.Z, BIC=self.BIC, BICs=self.BIC_evol, step=self.step, MAPs=self.MAPs,
                 logLLs=self.logLL_particles, minIncr=self.min_increase, formerT=self.former_T, currentT=self.current_T,
                 tempDone=self.temperingfinished, diffLL=self.diff_tar_impLL, corrstepsFct=self.corrsteps_fct,
                 nstepsMax=self.nsteps_max, nstepsMin=self.nsteps_min, ESStemp=self.ess_temp,tempStep=self.tempering_steps,
                 listESS=self.list_essTemp, listEXP=self.list_exponents, listUNI=self.list_unique, listMCMC=self.list_MCMC,
                 listHIGH=self.list_highest, listScale=self.list_scaling, rewDone=self.reweightingfinished, lastMCMC=self.lastMCMCdone,
                 totAcc=self.total_accepted)

    def update_mean_var(self):
        """
        updating the 'means' and 'variances' attributes according to current particles and weights
        >> no INPUT / no OUTPUT
        """
        self.means = np.array([np.dot(self.particles[i], self.weights) for i in range(self.d)])
        self.variances = np.array([np.dot((self.particles[i]-self.means[i])**2, self.weights) for i in range(self.d)])

    def resample(self): ### systematic resampling
        """
        resamples the particles from the current approximation using the systematic resampling algorithm,
        see e.g. A. Doucet and A. M. Johansen. A Tutorial on Particle filtering and smoothing: Fiteen years later. :
        https://warwick.ac.uk/fac/sci/statistics/staff/academic-research/johansen/publications/DJ11.pdf
        >> no INPUT / no OUTPUT
        """
        new_samples = []
        new_logLL = []
        new_diffLL = []

        r = random.uniform(0, 1./self.M)
        c = self.weights[0]
        i = 0
        for part in range(self.M):
            U = r + part / self.M
            while U > c:
                i += 1
                c += self.weights[i]
            new_samples += [self.particles[...,i]]
            new_logLL += [self.logLL_particles[i]]
            new_diffLL += [self.diff_tar_impLL[i]]

        self.particles = np.array(new_samples).T
        self.weights = np.full(self.M, 1.0/self.M)
        self.logLL_particles = np.array(new_logLL)
        self.diff_tar_impLL = np.array(new_diffLL)

    # def resample(self): ### random resampling (OLD method)
    #     """
    #     resamples the particles from the current approximation of the target distribution
    #      and reset the weighting of the particles to an uniform weighting
    #     >> no INPUT / no OUTPUT
    #     """
    #     indizes = np.random.choice(self.M, size=self.M, p=self.weights, replace=True)
    #     self.particles = self.particles[...,indizes]
    #     self.weights = np.full(self.M, 1.0/self.M)
    #     self.logLL_particles = self.logLL_particles[indizes]
    #     self.diff_tar_impLL = self.diff_tar_impLL[indizes] # tempering modification (addition) ### NEW v2
   
    def effective_sample_size(self):
        """
        computes effective sample size to approximate the effective number of samples of this
         approximation based of the deviation of the approximate weights variance

        >> no INPUT
        >> OUTPUT   effective sample size (float) of the current particle approximation
        """
        return 1. / np.dot(self.weights, self.weights)
        
    def vlogPDF(self, pars): ### pars = list of parameter vectors
        """
        paralellized computation of the vectorial logarithmic probability density funtion
         (see definition of function "logpdf" of class "parsPriors")
        >> INPUT    pars: array of shape [1:M,1:d], containing M d-dimensional parameter vectors
        >> OUTPUT   array of shape [1:M], containing the log-pdf values at the resp. parameter vectors in pars
        """
        pool = mp.Pool()
        tmp = pool.map(self.prior.logpdf, pars)
        pool.close()
        return np.array(tmp)

    # Tempering modification ### v2.5: remove 'step'
    def compute_diffpotentials(self, LL_fct):
        """
        computes difference potentials (updates for weights) for filtering steps (i.e. no multiplication by tempering exponent)

        >> INPUT
        step     filtering step
        LL_fct   potential function of the target distribution (same as target_potential)
                 (np.ndarray) -> np.ndarray

        >> no OUTPUT
        """

        ### compute importance-weights: tar_LL - imp_LL => LL_fct called with True
        self.diff_tar_impLL = LL_fct(self.particles)

    def reweight(self, num_data):
        """
        updates 'weight' attribute to approximate the target distribution:
         given the potential functions of a importance and target distribution,
         and assuming the current ParticleAprx approximates the importance distribution

        >> INPUT
        step                    filtering step - tempering modification
        importance_potential    potential function of the importance distribution
                                 (np.ndarray) -> np.ndarray
        target_potential        potential function of the target distribution
                                 (np.ndarray) -> np.ndarray

        >> no OUTPUT
        """
        ### for tempering, computation of diff_tar_impLL moved to routine compute_diffpotentials,
        ### which is called directly by smc_update, tempering modification

        # This is for computing model evidence
        if self.current_T == 1: # update evidence and BIC only at end of filtering steps, tempering modification
            self.Z += np.log10(np.dot(np.exp(self.diff_tar_impLL), self.weights)) ### log10(Z)
            self.BIC = -2. * np.max(self.diff_tar_impLL + self.logLL_particles) + self.d * np.log(num_data)

        ### update weights, tempering modification
        self.weights = np.exp(np.log(self.weights)+(self.current_T - self.former_T)*self.diff_tar_impLL)
        self.weights /= self.weights.sum() ### normalize

        if self.current_T == 1:
            ### MAP of unsmoothed, discrete particle approxmation (just for tracking)
            unique_particles, indices, counts = np.unique(self.particles, axis=1, return_counts=True, return_index=True)
            unique_weights = self.weights[indices]*counts
            tmp = list(zip(unique_particles.T, unique_weights)); tmp.sort(key=lambda x: x[-1])
            self.MAPs = tmp[-1][0]
        # non_unique = int(np.sum(counts[counts!=1]))
        # highest_count = np.max(counts)

        ### logL(data until step | particles), tempering modification
        self.logLL_particles = (self.current_T - self.former_T)*self.diff_tar_impLL + self.logLL_particles

        # return non_unique, highest_count

    def mh_correction(self, target_potential, proposal_kernel, title, first_step = False): ### v2.5: remove 'step'
        """
        improve the particle approximation by sampling several steps (here: n_steps)
         of the Metropolis-Hastings (MH) Markov Chain constructed using a given proposal kernel
         and the with the given target distribution as limiting distribution
        
        >> INPUT
        step                filtering step - tempering modification
        target_potential    potential function of the target distribution
                             (np.ndarray) -> np.ndarray
        proposal_kernel     function to sample proposals for the Metropolis-Hastings algorithm,
                            conditioned on the current particles, passed as a parameter
                             (np.ndarray) -> np.ndarray
            
        >> OUTPUT           the average acceptance rate over the correction steps
        """
        scaling_factor = 2. # adaptivity parameter: how much to increase/decrease the variance each time
        nMCMC_factor = 2. # adaptivity parameter for how to adjust number of MCMC correction steps adaptively
        
        # adaptive choice of number of MCMC applications 
        self.correction_steps = np.floor(nMCMC_factor/self.scaling**2)
        if self.correction_steps < self.nsteps_min:
            self.correction_steps = self.nsteps_min
        elif self.correction_steps > self.nsteps_max:
            self.correction_steps = self.nsteps_max

        print('Applying MCMC kernel {} times \n'.format(self.correction_steps))
            
        total_accepted = 0.
        
        for MCMCstep in range(int(self.correction_steps)):
            # Note: we are now using random walk MH for proposals. This part of the code might need changes for pCN (e.g., checking the prior support will not be needed anymore)
            
            ### Sample from the proposal kernel, conditioned on currect particles
            proposals = proposal_kernel(self)

            tmp = self.vlogPDF(list(np.concatenate([proposals.T, self.particles.T])))
            prop_pdf = np.array(tmp[:self.M]); part_pdf = np.array(tmp[-self.M:])
            prop_supp = np.isfinite(prop_pdf)

            ### current_potentials (after resampling) = logL_{t+1}(X)
            if np.sum(prop_supp) > 0:
                ### tempering modification
                if first_step:
                    proposal_diffpotentials = target_potential(proposals[..., prop_supp])
                    proposal_potentials = self.current_T * proposal_diffpotentials
                else:
                    prev, proposal_diffpotentials = target_potential(proposals[..., prop_supp])
                    proposal_potentials = self.current_T * proposal_diffpotentials + prev

                acceptance_ratio = np.zeros(self.M)

                ### Compute the acceptance ratio
                potential_ratio = proposal_potentials - self.logLL_particles[prop_supp] ### log(L_{t+1}(Xprop)/L_{t+1}(X))
                prior_ratio = prop_pdf[prop_supp] - part_pdf[prop_supp] ### log(prior(Xprop)/prior(X))

                ### accept_ratio of Xprop outside prior_supp stays zero; on supp:
                acceptance_ratio[prop_supp] = np.exp(potential_ratio + prior_ratio) ### (L_{t+1}*prior)(Xprop)/(L_{t+1}*prior)(X)

                ### Randomly accept the transitions based on the acceptance ratio
                tmp = np.random.uniform(size=self.M) ; accepted = tmp < acceptance_ratio; total_accepted += np.sum(accepted)
                self.particles[...,accepted] = proposals[...,accepted]

                ### update logLL, the difference potentials and priorPDF of mixed particles
                self.logLL_particles[prop_supp][accepted[prop_supp]] = proposal_potentials[accepted[prop_supp]]
                self.diff_tar_impLL[prop_supp][accepted[prop_supp]] = proposal_diffpotentials[accepted[prop_supp]] # tempering modification (addition)
                part_pdf[accepted] = prop_pdf[accepted]
                self.total_accepted = total_accepted; 
            else:
                print('all proposals outside prior support')

            self.lastMCMCdone += 1

            self.save(title)
            
        self.acceptance_rate = self.total_accepted / (self.correction_steps * self.M)
        
        # Adaptive tuning for next SMC iteration
        # adaptive choice of variance scaling
        if self.acceptance_rate > 0.3:
            print('Acceptance ratio: {} >0.3 \n'.format(self.acceptance_rate),
            " => variance scaling: {}->{}".format(self.scaling, self.scaling * scaling_factor))
            self.scaling *= scaling_factor
        elif self.acceptance_rate < 0.15:
            print('Acceptance ratio: {} <0.15 \n'.format(self.acceptance_rate),
            " => variance scaling: {}->{}".format(self.scaling, self.scaling / scaling_factor))
            self.scaling /= scaling_factor
        else:
            print('Acceptance ratio: {} ~ 0.25 \n'.format(self.acceptance_rate),
            " => keep variance scaling: {}".format(self.scaling))


    ### Tempering modification
    def compute_temperingexponent(self):
        """
        bisection algorithm for adaptive search of next temperature
        """

        if (1.-self.former_T)<= self.min_increase:
            self.current_T = 1.
            # self.temperingfinished = True
            # print('Tempering finsihed!')
        else:
            ### Bisection algorithm
            l = self.former_T ### left limit for bisection
            r = 1. ### right limit for bisection
            maxiter = 1000 ### max number of bisection iterations
            iter = 0 ### number of current iterations

            while iter<maxiter and (r-l)>self.min_increase:
                mid = (l+r)/2. # midpoint
                diffweights_tmp = (mid-self.former_T)*self.diff_tar_impLL # @Sabrina: please check sign of first summand!

                weights_tmp = np.exp(np.log(self.weights)+diffweights_tmp) ## update weights
                weights_tmp /= weights_tmp.sum() ### normalize
                ess_tmp = 1. / np.dot(weights_tmp, weights_tmp)
                # cess_tmp = self.M / np.sum(np.exp(2. * diffweights_tmp + np.log(self.weights)))

                if ess_tmp > self.M / self.ess_ratio: # mid can be larger: take half interval on the right
                    l = mid
                else: ### mid has to be smaller: take half interval on the left
                    r = mid

                iter = iter + 1

            self.current_T = mid
            if (1.-self.current_T) <= self.min_increase:
                self.current_T = 1.
                # self.temperingfinished = True
                # print('Tempering finished!')
        

    def smc_update(self, proposal_kernel, target_potential, num_data, title, first_step=False): ### v2.5: remove 'step'
        """
        updates current particle approximation using SMC with Metropolis-Hastings (MH) MCMC correction
        
        >> INPUT
        proposal_kernel         function to sample proposals for the MH algorithm, conditioned on the
                                current particles, passed as a parameter
                                 (np.ndarray) -> np.ndarray
        importance_potential    potential function of the importance distribution
                                 (np.ndarray) -> np.ndarray
        target_potential        potential function of the target distribution
                                 (np.ndarray) -> np.ndarray
        num_data                number of considered data point at the current SMC step - all are considered for tempering
            
        >> OUTPUT               ess: effective sample size after reweighting
        """
        ### in tempering, the difference potential is the same in all steps; this array is updated for new particles in mh_correction if needed
        self.compute_diffpotentials(target_potential)
        self.temperingfinished = False

        while not self.temperingfinished: ### iterate until the tempering exponent is 1
            print('------ TEMPERING-STEP {} --------------------------\n'.format(int(self.tempering_steps + 1)))

            if not self.reweightingfinished:
                self.compute_temperingexponent()
                self.list_exponents += [self.current_T]

                self.reweight(num_data)
                print(">> REWEIGHTED with exponent", self.current_T) ### keeping track of tempering exponent

                self.ess_temp = self.effective_sample_size(); self.list_essTemp += [self.ess_temp]

                ### keeping track of resampling
                ### we always resample
                if self.current_T != 1: 
                    self.resample()
                    print(">> RESAMPLED (ESS {}% , exponent < 1)\n".format(np.round(self.ess_temp * 100 / self.M, 2)))
                else:
                    print(">> NOT resampled (ESS {}% , exponent = 1)\n".format(np.round(self.ess_temp * 100 / self.M, 2)))
                self.reweightingfinished = True
                
                ### for debugging and algorithmic optimization purposes:
                ### keeping track of duplicates for...
                ### ... all parameters
                unique_particles, counts = np.unique(self.particles, axis=1, return_counts=True)
                non_unique = int(np.sum(counts[counts != 1]))
                highest_count = np.max(counts)
                print("--- duplicates: {}/{}  |  highest count: {}x".format(non_unique, self.M, highest_count))

                ### MH correction step
                self.mh_correction(target_potential, proposal_kernel, title, first_step=first_step)
                self.list_MCMC += [int(self.correction_steps)]
                self.list_scaling += [self.scaling]
                self.lastMCMCdone = 0
                print(">> MCMC complete!")
            
                ### for debugging and algorithmic optimization purposes:
                ### keeping track of duplicates for...
                ### ... all parameters
                unique_particles, counts = np.unique(self.particles, axis=1, return_counts=True)
                non_unique = int(np.sum(counts[counts != 1])); highest_count = np.max(counts)
                self.list_unique += [non_unique]
                self.list_highest += [highest_count]
                print("--- duplicates: {}/{}  |  highest count: {}x".format(non_unique, self.M, highest_count))

                self.former_T = self.current_T  ### update tempering exponent
                self.tempering_steps += 1
                self.reweightingfinished = False
                self.save(title); print("[saved intermediate result]")

                ### check if tempering is finished
                if self.current_T == 1:
                    self.temperingfinished = True
                    print('Tempering finished!')

            print('Tempering steps done: ',self.tempering_steps)

            self.list_essTemp = np.array(self.list_essTemp)
            self.list_exponents = np.array(self.list_exponents)
            self.list_unique = np.array(self.list_unique)
            self.list_MCMC = np.array(self.list_MCMC)
            self.list_highest = np.array(self.list_highest)
            self.list_scaling = np.array(self.list_scaling)

################################################

class SMC:
    def __init__(self, prior, get_loglikelihood, title, data, M, MC_steps = 1, num_data = 1, corr_steps=4, ess_ratio= 0.75, init=True):
        """
        CLASS OBJECT: SMC

        >> INPUT to initialize
        
        prior               prior object to sample from the prior and evaluate its logpdf (for RWMH steps)
        MC_steps            number of SMC steps, default = 1 for tempering
        get_loglikelihood   function calculating the log-likelihood of the parameters given model and data
        title               to save files
        data                tuple of ndarrays
        M                   number of particles
        num_data            list of number of data points, 1 for tempering
        corr_steps          number of applications of the MCMC kernel (per mutation step)
        ess_ratio           1/(ratio of M), at which particles should be resampled

        >> ATTRIBUTES
        aprx                corresponding ParticleAprx object
        title               (see INPUT)
        get_loglikelihood   (see INPUT)
        data                (see INPUT)
        num_data            (see INPUT)
        prior               (see INPUT)
        M                   (see INPUT)
        MC_steps            (see INPUT)
        corr_steps          (see INPUT)
        ess_ratio           (see INPUT)
        ESS                 effective sample size
        """
        if init:
            self.title = title
            self.get_loglikelihood = get_loglikelihood
            self.data = data; self.num_data = num_data

            self.prior = prior

            self.M = M; self.MC_steps = MC_steps
            self.corr_steps = corr_steps
            self.ess_ratio = 1/ess_ratio
            self.ESS = np.append(M,np.zeros(MC_steps))
            self.aprx = ParticleAprx(self.prior, self.M, self.corr_steps, self.ess_ratio)

    def smc_tempering(self, cal_pars=None):
        """
        performs the parameter calibration with SMC algorithm using tempering

        >> INPUT
        cal_pars        array of fixed parameter values

        >> no OUTPUT    (all information is stored in the SMC object in the end)
        """

        print("SMC for", self.title, "- Started:", time.ctime())
        
        MC_steps = 1

        ### initialize arrays for storing evidence/BIC evolution
        self.aprx.evidence = np.zeros(MC_steps+1) # log10(Z)
        self.aprx.BIC_evol = np.zeros(MC_steps)
        self.aprx.evidence[0] = self.aprx.Z
        with open(self.title + '_TMP.txt', 'w') as f:
            f.write('{}'.format(0))
            f.close()
        # plt_currHIST(self, 0, MAPs=True, save=True)
        with open(self.title+'_ESS.txt', 'a') as f:
            np.savetxt(f, [0,self.M,1,0,1,self.aprx.nsteps_min,1], fmt='%.5f', newline=" ")
            f.write("\n")
            f.close()

        with open(self.title + '_TMP.txt', 'r') as f:
            var_check = int(f.readlines()[0])
            f.close()

        start = time.time()

        ### Tempering 
        self.numdata = 1
        def proposal(aprx):
            """
                proposes particles using MCMC kernel (random-walk Metropolis)
                and using marginal variance scaled according to previous avg. acceptance rate

                >> INPUT    aprx : ParticleAprx object containing the current SMC approximation
                >> OUTPUT   array of proposed particles
            """
            # cov_mat = np.cov(aprx.particles,aweights=aprx.weights,ddof=0)
            # prop = stats.multivariate_normal(np.repeat(0.,aprx.d),cov_mat*aprx.scaling**2.).rvs(size=int(aprx.M))
            # return np.array(aprx.particles+prop.T)
            std = aprx.scaling * np.sqrt(aprx.variances)
            return np.array([aprx.particles[par] + stats.norm.rvs(loc=0.0, scale=std[par], size=aprx.M) for par in range(aprx.d)])
        
        self.aprx.smc_update(proposal, lambda x: self.get_loglikelihood(x), self.num_data,self.title, first_step= True)

        ### update and print current marginal mean/variance
        self.aprx.update_mean_var()
       
        ### save SMC approximation
        if (self.aprx.variances > 1e-5).all():
            self.aprx.save(self.title+'_tmp{}'.format(var_check))
        else:
            var_check += 1
            with open(self.title + '_TMP.txt', 'w') as f:
                f.write('{}'.format(var_check))
                f.close()
        self.aprx.save(self.title)

        with open(self.title + '_MAPs.txt', 'a') as f:
            for line in np.matrix(self.aprx.MAPs):
                np.savetxt(f, line, fmt='%.5f')
            f.close()

        with open(self.title + '_means.txt', 'a') as f:
            for line in np.matrix(self.aprx.means):
                np.savetxt(f, line, fmt='%.5f')
            f.close()

        with open(self.title + '_vars.txt', 'a') as f:
            for line in np.matrix(self.aprx.variances.T):
                np.savetxt(f, line, fmt='%.5f')
            f.close()
    print("Finished:", time.ctime())

