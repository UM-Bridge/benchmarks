# MLMC function for level l

from numpy import random, zeros, array
from numpy import sum as npsum
from time import time
from concurrent.futures import ThreadPoolExecutor, wait, as_completed, FIRST_EXCEPTION


def futures_output(futures):
    """
    Blocking function to retrieve parallel model outputs
    Raise exception and cancel all jobs when error is caught
    """
    outputs = []
    completed, running_or_failed = wait(futures, return_when=FIRST_EXCEPTION)
    for future in completed:
        try:
            future.result()
        except Exception as e:
            for job in futures:
                job.cancel()
            raise e
        else:
            outputs.append(future.result())

    for future in running_or_failed:
        try:
            future.result()
        except Exception as e:
            for job in futures:
                job.cancel()
            raise e
        else:
            outputs.append(job.result())
    return array(outputs)


def mlmc_fn(l, N, problems, coupled_problem=False, sampler=None, N1 = 1):
    """
    Inputs:
        l: level
        N: number of paths
        problems: list of problems
            problems[l-1]: application-specific coarse problem (for l>0)
            problems[l]: application-specific fine problem 
            Problems must have an evaluate method such that
            problems[l].evaluate(sample) returns output P_l.
            Optionally, user-defined problems.cost
        coupled_problem: if True,
             problems[l].evaluate(sample) returns both P_l and P_{l-1}.
        sampler: sampling function, by default standard Normal.
            input: N, l
            output: (samplef, samplec). The fine and coarse samples.
         N1: number of paths to generate concurrently.

    Outputs:
        (sums, cost) where sums is an array of outputs:
        sums[0] = sum(Pf-Pc)
        sums[1] = sum((Pf-Pc)**2)
        sums[2] = sum((Pf-Pc)**3)
        sums[3] = sum((Pf-Pc)**4)
        sums[4] = sum(Pf)
        sums[5] = sum(Pf**2)
        cost = user-defined computational cost. By default, time
    """

    if sampler is None:
        def sampler(N, l):
            sample = random.randn(N)
            return (sample, sample)

    sums = zeros(6)
    cpu_cost = 0.0
    problemf = problems[l]
    if l>0:
        problemc = problems[l-1]
   
    samplesc, samplesf = [], []
    
    for i in range(1, N+1, N1):
        N2 = min(N1, N - i + 1)

        samplef, samplec = sampler(N2,l)
        samplesf.append(samplef)
        samplesc.append(samplec)
    
    workers = 1 # Change to match number of nodes
    with ThreadPoolExecutor(max_workers=workers) as executor:
        """
        Pf, Pc are np arrays now
        """
        start = time()
        if coupled_problem:
            futures = {executor.submit(problems[l].evaluate, samplesf[i-1], i): i for i in range(1, N+1, N1)}
            Ps = futures_output(futures)
            Pf = Ps[:, 0]
            Pc = Ps[:, 1]
        else:
            futures = {executor.submit(problemf.evaluate, samplesf[i-1], i): i for i in range(1, N+1, N1)}
            Pf = futures_output(futures)
            if l == 0:
                Pc = 0.
            else:
                futures = {executor.submit(problemc.evaluate, samplesc[i-1], i): i for i in range(1, N+1, N1)}
                Pc = futures_output(futures)
            
        end = time()
        cpu_cost += (end - start) * workers # cost defined as total computational time
        sums += array([npsum(Pf - Pc),
                       npsum((Pf - Pc)**2),
                       npsum((Pf - Pc)**3),
                       npsum((Pf - Pc)**4),
                       npsum(Pf),
                       npsum(Pf**2)])

        
    problem_cost_defined = hasattr(problemf, 'cost')
    problem_cost_defined = problem_cost_defined and problemf.cost is not None

    if problem_cost_defined:
        cost = N*problemf.cost
        if l>0:
            cost += N*problemc.cost # user-defined problem-specific cost
    else:
        cost = cpu_cost

    return (sums, cost)
