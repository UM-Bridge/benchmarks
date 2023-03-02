# %% Based on CUQIpy version v0.3.0
import numpy as np
import cuqi
import matplotlib.pyplot as plt
import scipy

# Check cuqi version
#assert cuqi.__version__ == "0.3.0"
print ("CUQIpy version: ", cuqi.__version__)

# %%Utilities
def mysavefig(filename):
    plt.savefig(filename,
                bbox_inches='tight', 
                dpi=300)

# %% Benchmark data generation


# %% CUQI imports
from cuqi.model import LinearModel
from cuqi.geometry import Image2D
from cuqi.array import CUQIarray
from cuqi.distribution import Gaussian, GMRF, Laplace_diff, Cauchy_diff
from cuqi.problem import BayesianProblem

# %%  Specification of image grid

N = 256

lower = -1.0
upper = 1.0
width = upper - lower

dx = width/N
dx5 = dx/2.0

xx = np.linspace(lower+dx5,upper-dx5,N)

X, Y = np.meshgrid(xx,xx)

# %%  Function to create filled disk with value v inside, center (x0.y0), radius r
def disk_fun(x0, y0, r, v):
    im = np.zeros((N,N))
    im[(X - x0)**2 + (Y-y0)**2 <= r**2] = v
    return im

# %% Outer disk
im1 = disk_fun(0.025, 0.05, 0.8, 0.5)
im1

plt.imshow(im1)
plt.colorbar()

# %%  Add 3x3 disks, 3 radii, 3 contrasts
r0 = 0.025
r1 = 0.05
r2 = 0.1
r3 = 0.2

va = 0.1
vb = 0.2
vc = 0.4
#d0c = disk_fun( 0.5, 0.1, r0,-vc)

d0a = disk_fun( 0.3, 0.0, r0,-va)
d0b = disk_fun( 0.1,-0.5, r0, vb)
d0c = disk_fun(-0.5,-0.3, r0,-vc)

d1a = disk_fun( 0.0,-0.1, r1,-va)
d1b = disk_fun( 0.0, 0.7, r1, vb)
d1c = disk_fun( 0.4, 0.6, r1, vc)

d2a = disk_fun(-0.2,-0.4, r2, va)
d2b = disk_fun(-0.3, 0.5, r2,-vb)
d2c = disk_fun( 0.6, 0.2, r2,-vc)

d3a = disk_fun( 0.1, 0.3, r3, va)
d3b = disk_fun( 0.4,-0.3, r3,-vb)
d3c = disk_fun(-0.5, 0.1, r3, vc)

imall = im1 + d1a + d1b + d1c + d2a + d2b + d2c + d3a + d3b + d3c
imall = im1 + d0a + d0b + d0c + d1a + d1b + d1c + d2a + d2b + d2c + d3a + d3b + d3c

plt.imshow(imall)
plt.colorbar()

# %%  Reshape image into vector for use with sparse matrix as forward model.
imall.shape = N**2

# %%  Choose how many angles, load system matrix and scale by pixel size.
nv = 30

Amat = scipy.io.loadmat('A' + str(N) + '_' + str(nv) + '.mat')["A"]

Amat = dx*Amat

# %% Create the linear model

A = LinearModel(Amat)

# %%  Create the visual_only Image2D geometries of image and sinogram spaces
dg = Image2D((N,N),visual_only=True)
rg = Image2D((nv,N),visual_only=True)

# %% Equip linear operator with geometries
A.domain_geometry = dg
A.range_geometry = rg

# %% Create the CUQI data structure from vectorized image and geometry
imC = CUQIarray(imall,geometry=dg)
imC.plot()
plt.colorbar()

# %%  Demonstrate creation of sinogram with CUQI 
sinoC = A(imC)
sinoC.plot()
plt.colorbar()

# %%  Backprojection possible as adjoint.
bp = A.adjoint(sinoC)
bp.plot()
plt.colorbar()

# %% Gaussian prior - have to provide geometry to make work

x = Gaussian(mean=np.zeros(A.domain_dim), 
                           cov=0.01,
                           geometry=A.domain_geometry)

# %%  Choose noise std and create data distribution
s = 0.01
y = Gaussian(A@x, s**2)

# %%  Generate data and equip with missing geometry
np.random.seed = 0

y_data = y(x=imC).sample()
y_data.geometry = rg

# %%  Display data, clean data, and noise
y_data.plot()
plt.colorbar()

# %%

sinoC.plot()
plt.colorbar()

#%%
(y_data-sinoC).plot()
plt.colorbar()




#%%  Run for Gaussian prior

BP = BayesianProblem(y, x).set_data(y=y_data)
print(BP)

samples = BP.sample_posterior(1000)

samples.plot_ci(exact=imC)




# %% Run for GMRF prior

x = GMRF(mean=np.zeros(A.domain_dim), 
                           prec=1/0.01,
                           physical_dim=2,
                           geometry=dg) 

y = Gaussian(A@x, s**2)

BP = BayesianProblem(y, x).set_data(y=y_data)
print(BP)

samples = BP.sample_posterior(1000)

samples.plot_ci(exact=imC)




# %% Run for Laplace_diff prior

x = Laplace_diff(np.zeros(A.domain_dim), 
                                    scale=0.1,
                                    physical_dim=2,
                                    geometry=dg)

y = Gaussian(A @ x, s**2)

BP = BayesianProblem(y, x).set_data(y=y_data)

print(BP)

samples = BP.sample_posterior(1000)

#%% 
samples.plot_mean(), plt.colorbar()
mysavefig('LaplaceDiff_mean.png')


#%%
samples.plot_std(), plt.colorbar()
mysavefig('LaplaceDiff_std.png')




# %%  Run with Cauchy_diff prior

x = Cauchy_diff(np.zeros(A.domain_dim), 
                scale=0.01,
                physical_dim=2,
                geometry=dg)

y = Gaussian(A @ x, s**2)

BP = BayesianProblem(y, x).set_data(y=y_data)
print(BP)

samples = BP.sample_posterior(500)

samples.plot_ci(exact=imC)



#%%

np.random.seed = 0

TP = cuqi.testproblem.Deconvolution1D(
    dim=128,
    PSF="gauss",
    PSF_param=5,
    phantom="square",
    noise_std=0.01,
    noise_type="gaussian",
)

# Save data code (uncomment to overwrite data file)
#np.savez("data_square.npz", data=TP.data, exact=TP.exactSolution)

# %% Plot data

TP.data.plot()
plt.title("Data")

# %% Plot exact solution

TP.exactSolution.plot()
plt.title("Exact solution")

# %% Plot posterior samples for Gaussian

TP.prior = cuqi.distribution.Gaussian(np.zeros(128), 0.01, name="x")
samples = TP.UQ()
plt.title("Gaussian with delta=0.01")

# %% Plot posterior samples for GMRF

TP.prior = cuqi.distribution.GMRF(np.zeros(128), 1/0.01, name="x")
samples = TP.UQ()
plt.title("GMRF with delta=0.01")

# %% plot posterior samples for LMRF

TP.prior = cuqi.distribution.Laplace_diff(np.zeros(128), 0.01, name="x")
samples = TP.UQ()
plt.title("LMRF with delta=0.01")

# %% Plot posterior samples for CMRF

TP.prior = cuqi.distribution.Cauchy_diff(np.zeros(128), 0.01, name="x")
samples = TP.UQ()
plt.title("CMRF with delta=0.01")


# %%