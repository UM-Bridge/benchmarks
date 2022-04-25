import numpy as np
import h5py

# Import forward model class
from BeamModel import EulerBernoulli

# MUQ Includes
import muq.Modeling as mm 


# Load the problem definition and data ("GenerateObservations.py" must have been previously run to generate this)
f = h5py.File('ProblemDefinition.h5','r')

x = np.array( f['/ForwardModel/NodeLocations'] )
B = np.array( f['/Observations/ObservationMatrix'] )
obsData = np.array( f['/Observations/ObservationData'] )

length = f['/ForwardModel'].attrs['BeamLength']
radius = f['/ForwardModel'].attrs['BeamRadius']

loads = np.array( f['/ForwardModel/Loads'])

numObs = obsData.shape[0]
numPts = x.shape[1]
dim = 1


# Create a linear operator that maps the 3 lumped parameter values to the entire field
numIntervals = 3
endPts = np.linspace(0,1,numIntervals+1)
intervals = [(endPts[i],endPts[i+1]) for i in range(numIntervals)]

A = np.zeros((numPts,numIntervals))
for i in range(len(intervals)):
    A[(x[0,:]>=intervals[i][0]) & (x[0,:]<=intervals[i][1]), i] = 1.0

# Define the prior distribution 
logPriorMu = 10*np.ones(numIntervals)
logPriorCov = 4.0*np.eye(numIntervals)

logPrior = mm.Gaussian(logPriorMu, logPriorCov).AsDensity()

# Create some pieces that will be needed in the graph
mField      = mm.DenseLinearOperator(A)
expmVals    = mm.ExpOperator(numIntervals)
loadPiece   = mm.ConstantVector(loads)
obsOperator = mm.DenseLinearOperator(B)

# EulerBernoullie is a child of ModPiece with two inputs:
# 1. A vector of loads at each finite difference node
# 2. A vector containing the material property (exp(m(x))) at each finite difference node
beamModel = EulerBernoulli(numPts, length, radius)



# Define the likelihood function 
noiseVar = 1e-4
noiseCov = noiseVar*np.eye(obsData.shape[0])
likelihood = mm.Gaussian(obsData, noiseCov).AsDensity()


# Construct the posterior graph 
posteriorPiece = mm.DensityProduct(2)
mPiece = mm.IdentityOperator(numIntervals)

graph = mm.WorkGraph()

# Forward model nodes and edges
graph.AddNode(mPiece, "m_i")
graph.AddNode(expmVals, "exp(m_i)")
graph.AddNode(mField, "exp(m(x))")
graph.AddNode(loadPiece, "f")
graph.AddNode(obsOperator, "B")
graph.AddNode(beamModel, "u")

graph.AddEdge("m_i", 0, "exp(m_i)", 0)
graph.AddEdge("exp(m_i)", 0, "exp(m(x))", 0)
graph.AddEdge("exp(m(x))", 0, "u", 1)
graph.AddEdge("f", 0, "u", 0)
graph.AddEdge("u", 0, "B", 0)

# Other nodes and edges
graph.AddNode(likelihood, "Likelihood")
graph.AddNode(logPrior, "Prior")
graph.AddNode(posteriorPiece,"Posterior")

graph.AddEdge("B", 0, "Likelihood", 0)
graph.AddEdge("m_i", 0, "Prior", 0)
graph.AddEdge("Prior",0,"Posterior",0)
graph.AddEdge("Likelihood",0, "Posterior",1)


## Serve up the log posterior density on port 4243
logPost = graph.CreateModPiece("Posterior")

mm.serveModPiece(logPost, "0.0.0.0", 4243)