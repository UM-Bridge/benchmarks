from dolfin import *
import logging
import sys
import gc
import numpy as np
from petsc4py import PETSc


class CookiePDE:
    """
    Class representing an elliptic or parabolic partial differential equation with ``cookie'' diffusion coefficient.
    """

    def __init__(self, N, basis_degree):
        """
        Initialize the PDE problem.

        Parameters:
        N (int): Number of mesh divisions in each direction.
        """
        sh = logging.StreamHandler(sys.stdout)
        logging.basicConfig(level=logging.DEBUG, handlers=[sh])
        set_log_active(True)
        set_log_level(LogLevel.DEBUG)

        # Create FEniCS mesh and define function space
        # Define the corners of the rectangular domain
        x0, y0 = 0.0, 0.0  # Bottom-left corner
        x1, y1 = 1.0, 1.0  # Top-right corner

        # Define the number of elements in each direction
        nx, ny = N, N

        p0 = Point(x0,y0)
        p1 = Point(x1,y1)
        # Create a rectangular mesh with quadrilateral elements
        mesh = RectangleMesh.create([p0,p1],[nx,ny],CellType.Type.quadrilateral)

        # Can alternatively use a triangular mesh
        # mesh = UnitSquareMesh(N, N)

        # Define approximation space
        V = FunctionSpace(mesh, "Lagrange", basis_degree)
        self.mesh = mesh
        self.V = V

    def setupProblem(self, difftype, y, quad_degree=8, varcoeffs=None, advection=0):
        """
        Set up the variational problem for the PDE.

        Parameters:
        difftype (str): Type of diffusion field (only 'cookie' implemented).
        y (array): Array of coefficients for the diffusion field.
        quad_degree (int, optional): Quadrature degree for integration (default is 8).
        varcoeffs (array, optional): Array of variable coefficients (default is None).
        advection (integer, optional): Whether to include advection term (default is 0 (false)).
        """
        V = self.V

        # Define boundary condition
        u0 = Constant(0.0)

        # Homogeneous Dirichlet BC for unit square domain
        def boundary(x):
            return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS or x[1] < DOLFIN_EPS or x[1] > 1.0 - DOLFIN_EPS

        bc = DirichletBC(V, u0, boundary)

        # Define variational problem
        u = TrialFunction(V)
        v = TestFunction(V)

        # Forcing on subdomain F defined by f_indicator
        # F is an indicator function for the square subdomain [0.4,0.6]^2
        # NOTE: I think using SpatialCoordinate avoids interpolating the discontinuity into the function space.
        x = SpatialCoordinate(self.mesh)
        f_indicator = conditional(lt(abs(x[0] - 0.5),0.1), 1.0, 0.0)*conditional(lt(abs(x[1] - 0.5),0.1), 1.0, 0.0)
        # Define f using the indicator function
        f = 100 * f_indicator

        # Diffusion field defined via indicators on appropriate subdomains
        # See the definition of the benchmark problem (https://doi.org/10.48550/arXiv.2402.13768)
        if difftype == 'cookie':
            # Input varcoeffs allows the background diffusion field to be changed
            if varcoeffs is None:
                a0 = 1.0
            else:
                a0 = varcoeffs[0]

            # Cookie centres
            cxlist = [0.2, 0.5, 0.8, 0.2, 0.8, 0.2, 0.5, 0.8]
            cylist = [0.2, 0.2, 0.2, 0.5, 0.5, 0.8, 0.8, 0.8]

            # Cookie radius
            r = 0.13
            diff = a0
            for ii in range(8):
                cx_ii = cxlist[ii]
                cy_ii = cylist[ii]
                distance = sqrt(pow(x[0] - cx_ii, 2) + pow(x[1] - cy_ii, 2))
                indicator_ii = conditional(lt(distance, r), 1.0, 0.0)
                diff = diff + indicator_ii * Constant(1.0-a0 + y[ii])
        else:
            error("Unknown diffusion field")

        # Diffusion bilinear form
        dx_q = Measure('dx', domain=self.mesh, metadata={'quadrature_degree': quad_degree})
        a = inner(diff*grad(u), grad(v))*dx_q

        # Optional: additional advection term. This is not included in the benchmark
        if advection==1:
            w = as_vector((4*(x[1]-0.5)*(1-4*pow(x[0]-0.5,2)), -4*(x[0]-0.5)*(1-4*pow(x[1]-0.5,2))))
            a = a + inner(w, grad(u)) * v * dx_q

        # Define linear form on RHS
        L = f*v*dx_q

        # Bilinear form for mass matrix
        m = inner(u, v)*dx_q

        # Generic function in approx space V
        u = Function(V)

        # Write to object
        # self.diffproject = project(diff, V)
        # self.fproject = project(f, V)
        self.f_indicator = f_indicator
        self.a = a
        self.L = L
        self.bc = bc
        self.m = m
        self.u = u
        self.y = y

    def solve(self, directsolver, pctype, tol):
        """
        Solve the elliptic PDE problem using PETSc.

        Parameters:
        directsolver (int): Use direct LU solver if equal to 1, otherwise GM-RES iterative solver.
        pctype (str): Type of preconditioner to use with iterative solver ('ILU', 'JACOBI', or 'none').
        tol (float): Tolerance for the iterative solver (if used).

        Returns:
        numpy.ndarray: Solution vector.
        """
        # Assemble FEniCS matrices
        A = assemble(self.a)
        b = assemble(self.L)
        self.bc.apply(A)
        self.bc.apply(b)

        print("Solving for y=" + str(self.y))

        # Convert from fenics to petsc
        A_petsc = as_backend_type(A).mat()
        b_petsc = as_backend_type(b).vec()
        u_petsc = as_backend_type(self.u.vector()).vec()

        # Construct solver
        ksp = PETSc.KSP().create()
        pc = PETSc.PC().create()

        if directsolver == 1:
            # Direct solver is full LU preconditioner (and hence no iterations required)
            pc.setType(PETSc.PC.Type.LU)  # Set the preconditioner type
            ksp.setType(PETSc.KSP.Type.PREONLY)  # Choose a solver type
            pc.setOperators(A_petsc)  # Attach the matrix to the preconditioner
            ksp.setPC(pc)
        else:
            # Use GMRES to approximate solution of linear system
            ksp.setType(PETSc.KSP.Type.GMRES)  # Choose a solver type
            ksp.setTolerances(atol=tol)  # Set tolerance

            # Preconditioning
            if pctype == "ILU":
                pc.setType(PETSc.PC.Type.ILU)  # Set the preconditioner type
            elif pctype == "JACOBI":
                pc.setType(PETSc.PC.Type.JACOBI)  # Set the preconditioner type
            else:
                pass
            pc.setOperators(A_petsc)  # Attach the matrix to the 
            ksp.setPC(pc)
        # Set up solver operator
        ksp.setOperators(A_petsc)  # Set the operator
        ksp.setFromOptions()  # Set PETSc options for the solver


        # Solve the linear system
        ksp.solve(b_petsc, u_petsc)  # Solve A*x = b for x
        iternum = ksp.getIterationNumber()
        print("Solved in: "+str(iternum)+" iterations")

        # Write vector to FEniCS vector
        self.u.vector()[:] = u_petsc[:]

        # Destroy petsc objects
        A_petsc.destroy()
        b_petsc.destroy()
        u_petsc.destroy()
        ksp.destroy()
        pc.destroy()

        # Return dof as vector
        return self.u.vector().get_local()

    def projectref(self, n_ref):
        """
        Interpolate the solution to a finer mesh.

        Parameters:
        n_ref (int): Number of mesh divisions in each direction for the finer mesh.

        Returns:
        numpy.ndarray: Interpolated solution vector on the finer mesh.
        """
        # Define the reference mesh
        N = n_ref
        mesh = UnitSquareMesh(N, N)
        V = FunctionSpace(mesh, "Lagrange", 1)
        # Interpolate on mesh
        u_int = interpolate(self.u, V)
        # Return values at nodes
        return u_int.vector().get_local()

    def computebenchmarkqoi(self):
        """
        Compute the benchmark quantity of interest.

        Returns:
        float: Integral of the solution over the region F.
        """
        integral = assemble(self.f_indicator*self.u*dx)
        return integral

    def computenorm(self, x, normtype):
        """
        Compute the norm of the solution.

        Parameters:
        x (array): Solution vector.
        normtype (str): Type of norm to compute.

        Returns:
        float: Computed norm of the solution.
        """
        u = Function(self.V)
        u.vector().set_local(x[:])
        return norm(u, normtype)

    def writeSln(self, filename):
        """
        Write the solution to a Paraview file.

        Parameters:
        filename (str): Base name for the output file.
        """
        fileout = File(filename + ".pvd")
        fileout << self.u
        # fileout << self.diffproject
        # fileout << self.fproject


    def solveTime(self, tol, finalTime):
        """
        Solves the time-dependent problem using PETSc's TS solver.

        Args:
            tol (float): The relative tolerance for the solver.
            finalTime (float): The final time up to which the solution is computed.

        Returns:
            numpy.ndarray: The local part of the solution vector.
        """
        # Compute solution
        A = assemble(self.a)
        b = assemble(self.L)
        M = assemble(self.m)
        self.bc.apply(A)
        self.bc.apply(b)
        self.bc.apply(M)

        print("Solving for y=" + str(self.y))

        # Set up PETSc formulation
        A_petsc = as_backend_type(A).mat()
        M_petsc = as_backend_type(M).mat()
        b_petsc = as_backend_type(b).vec()
        u_petsc = as_backend_type(self.u.vector()).vec()

        f = u_petsc.copy()
        fm = u_petsc.copy()

        ts = PETSc.TS().create()

        def rhs_function(ts, t, u, F, A, b):
            A.mult(u, F)
            F.axpby(1.0, -1.0, b)
            return

        def rhs_function_specific(ts, t, u, F):
            rhs_function(ts, t, u, F, A_petsc, b_petsc)

        ts.setRHSFunction(rhs_function_specific, f)

        def mass_matrix_multiply(ts, t, u, udot, F, M):
            M.mult(udot, F)  # F = M * u

        def mass_matrix_multiply_specific(ts, t, u, udot, F):
            mass_matrix_multiply(ts, t, u, udot, F, M_petsc)

        ts.setIFunction(mass_matrix_multiply_specific, fm)

        vtkfile = File("output.pvd")  # Output file name (with .pvd extension)

        def monitor(ts, step, t, u):
            print(f"Time Step {step}: Time = {t}")
            self.u.vector()[:] = u[:]
            self.u.rename("u", "label")
            vtkfile << (self.u, t)  # Write function u to VTK file

        # Set the monitor function for the TS solver
        ts.setMonitor(monitor)

        # Create a PETSc options object
        opts = PETSc.Options()

        # Set the adaptive basis type (e.g., default is "basic")
        opts.setValue('-ts_adapt_type', 'basic')
        # Other adaptive basis options can be set here

        # Set solver options for adaptive time-stepping

        ts.setProblemType(ts.ProblemType.LINEAR)
        # ts.setEquationType(ts.EquationType.IMPLICIT)
        ts.setType(ts.Type.BEULER)

        ts.setSolution(u_petsc)
        ts.setMaxTime(float(finalTime))

        ts.setTolerances(rtol=tol)
        ts.setTimeStep(1e-8)
        ts.setFromOptions()

        # ts.setTimeStep(1e-8)
        ts.solve(u_petsc)

        PETScOptions.clear()
        ts.reset()
        ts.destroy()

        self.u.vector()[:] = u_petsc[:]

        # Destroy PETSC objects
        A_petsc.destroy()
        M_petsc.destroy()
        b_petsc.destroy()
        u_petsc.destroy()
        f.destroy()
        fm.destroy()

        gc.collect()
        # PETSc.Log.view()

        # Return solution vector at time T
        return self.u.vector().get_local()


    def solveTimeSimple(self, tol, finalTime):
        """
        Solves the time-dependent problem using a simple time-stepping approach.

        Args:
            tol (float): The relative tolerance for the solver.
            finalTime (float): The final time up to which the solution is computed.

        Returns:
            numpy.ndarray: The local part of the solution vector.
        """
        # PETSc.Options().setValue('-log_view', '')
        # PETSc.Options().setValue('-malloc_view', '')
        # PETSc.Log.begin()

        # Compute solution
        A = assemble(self.a)
        b = assemble(self.L)
        M = assemble(self.m)
        self.bc.apply(A)
        self.bc.apply(b)
        self.bc.apply(M)

        print("Solving for y=" + str(self.y))

        A_petsc = as_backend_type(A).mat()
        M_petsc = as_backend_type(M).mat()
        b_petsc = as_backend_type(b).vec()
        u_petsc = as_backend_type(self.u.vector()).vec()
        A = None
        b = None
        M = None
        gc.collect()

        # Output file name (with .pvd extension)
        vtkfile = File("output" + str(tol) + ".pvd")

        def monitor(ts, step, t, u, u2, enorm):
            print(f"Time = {t:.4g},\t dt = {ts:.4g},\t E = {enorm:.4g},\t acc = {acc:.4g}")
            # self.u.vector()[:] = u[:]
            # self.u.rename("u", "label")
            # vtkfile << (self.u, t)  # Write function u to VTK file

        ksp = PETSc.KSP().create()

        # Initialise
        t = 0.0
        dt = 1e-5
        dtold = 1e9

        lhs = A_petsc.copy()
        lhs_old = A_petsc.copy()
        rhs = u_petsc.copy()
        rhs_temp = u_petsc.copy()
        abf2 = u_petsc.copy()
        u_petsc_m1 = u_petsc.copy()
        e_petsc = u_petsc.copy()
        e_temp = u_petsc.copy()

        ii = 0
        # Timestepping loop.
        # Timestep dt is adapted using TR-AB2 pair and the Milne device to estimate local error.
        # Solution is returned at finalTime
        while t < finalTime:
            if t + dt > finalTime:
                dt = finalTime - t

            lhs.zeroEntries()
            lhs_old.zeroEntries()
            rhs.zeroEntries()

            # lhs_old = M_petsc + 0.5*dt*A_petsc
            lhs.axpy(0.5 * dt, A_petsc)
            lhs.axpy(1.0, M_petsc)

            # rhs = (M_petsc - 0.5*dt*A_petsc) * u_petsc + dt*b_petsc
            M_petsc.mult(u_petsc, rhs)
            A_petsc.mult(u_petsc, rhs_temp)
            rhs.axpy(-0.5 * dt, rhs_temp)
            rhs.axpy(dt, b_petsc)

            ksp.setOperators(lhs)  # Set the operator

            # abf2 = u_petsc + dt * (-1.5 * A_petsc * u_petsc + 0.5 * A_petsc * u_petsc_m1)
            A_petsc.mult(u_petsc, abf2)
            A_petsc.mult(u_petsc_m1, rhs_temp)
            abf2.axpby(0.5 * (dt / dtold), -(1 + 0.5 * dt / dtold), rhs_temp)
            abf2.axpby(1.0, dt, u_petsc)
            abf2.axpy(dt, b_petsc)

            u_petsc_m1.axpby(1.0, 0.0, u_petsc)

            # Solve linear system for timestep
            ksp.solve(rhs, u_petsc)

            # print(u_petsc[:])

            # Estimate local timestepping error
            e_petsc.zeroEntries()
            e_petsc.axpby(1.0, 0.0, abf2)
            e_petsc.axpy(-1.0, u_petsc)
            M_petsc.mult(e_petsc, e_temp)
            enorm = (dt / (3 * (1 + dtold / dt))) * np.sqrt(e_petsc.dot(e_temp))
            ii += 1
            t += dt

            # Compute timestep acceleration
            acc = float(np.power(tol / enorm, 2 / 3))
            if acc > 10:
                acc = 10  # Limit to exponential growth
            monitor(dt, ii, t, u_petsc, abf2, enorm)
            dtold = dt
            if ii != 1:
                dt = dt * acc  # Don't adapt first step

            ksp.reset()

        # Clean up
        ksp.destroy()

        self.u.vector()[:] = u_petsc[:]

        PETScOptions.clear()

        A_petsc.destroy()
        M_petsc.destroy()
        b_petsc.destroy()
        u_petsc.destroy()

        lhs.destroy()
        lhs_old.destroy()
        rhs.destroy()

        # Return approximation for finalTime T
        return self.u.vector().get_local()