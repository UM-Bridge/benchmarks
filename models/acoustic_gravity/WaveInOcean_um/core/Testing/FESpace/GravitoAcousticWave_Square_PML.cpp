
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;

// Physical pqrqmeters
Real bdny_density_value = 400.0;
Real density_value = 64.0 / 9.0;
Real length = 2;
Real width = 1;

// Time sampling
Real T_end = 40.0;
Real OutputDeltat = 0.2;

// Parameters value for in the pml

Real xi = 50.0;
Real Size_PML = 0.1;

Index Labels_PML(const RealVector &p)
{

    if (p[0] > length - Size_PML || p[0] < Size_PML)
        return 1;

    return 0;
}

int main(int argc, char *argv[])
{

    std::cout << "Running AcousticWave_Square........................" << std::endl;

    std::cout << "Simulation of a transversaly isotropic acoustic wave in a square." << std::endl;

    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;

    Scaling<2> scaling(length, width);
    Square<Scaling<2>> fespace(10, 10, 40, 20, scaling);

    auto surface_fespace = fespace.getBndyFESpace(2);
    auto bottom_fespace = fespace.getBndyFESpace(0);

    bottom_fespace.setAsTraceSpace();
    surface_fespace.setAsTraceSpace();
    surface_fespace.setLabels(Labels_PML);

    Index N = fespace.getNumDoFs();
    Index M = fespace.getNumDoFsDisc();

    // Create the subdomain where the pml are defined

    auto fespace_pml = fespace.getSubDomainFESpace([](RealVector xyz)
                                                   {if (xyz[0] > length - Size_PML || xyz[0] < Size_PML)
                            return true;
                    else return false; });

    fespace_pml.setAsSubSpace(); // Non optimal version -> auxiliary unknown stored everywhere

    std::cout << "The simulation has " << N << " degree of freedom" << std::endl;

    Field::Scalar one(1.0);
    Field::Scalar density(density_value);
    Field::Scalar bndy_density(bdny_density_value);
    std::map<Index, Real> xi_map;

    xi_map[0] = 0.0;
    xi_map[1] = bdny_density_value * xi;

    Field::Scalar<Field::Regularity::PiecewiseConstant> field_xi(xi_map);

    LAL::DiagonalMatrix massMatrix;
    LAL::DiagonalMatrix surface_massMatrix;
    LAL::DiagonalMatrix surface_massMatrix_dissip;
    LAL::DiagonalMatrix bottom_massMatrix;
    LAL::DiagonalMatrix MassTotal;

    CreateMass(fespace, density, massMatrix, 1.0);
    CreateMassBndy(surface_fespace, bndy_density, surface_massMatrix, 1.0);

    CreateMassBndy(surface_fespace, field_xi, surface_massMatrix_dissip, 1.0);
    CreateMassBndy(bottom_fespace, one, bottom_massMatrix, 1.0);

    MassTotal = massMatrix + surface_massMatrix;

    Index N_iter;
    Real Residue;
    Real rho;

    rho = LAL::PowerIteration([&MassTotal](LAL::Vector &v)
                              { LAL::Solve(MassTotal, v); },

                              [&fespace, &one](LAL::Vector &u, LAL::Vector &v)
                              { MltAddStiffness(fespace, one, 1.0, u, 0.0, v); },

                              N, N_iter, Residue);

    Real dt = 0.95 * (2.0 / sqrt(rho));

    std::cout << "Using " << N_iter << " iterations of a Power Iteration method, the estimated time step is" << std::endl;

    std::cout << "\t" << dt << std::endl;

    // Operators for the PML

    RealMatrix2x2 dissip_mat;
    RealMatrix2x2 coupling_mat;

    coupling_mat[0][0] = -xi;
    coupling_mat[0][1] = 0.0;
    coupling_mat[1][0] = 0.0;
    coupling_mat[1][1] = +xi;

    dissip_mat[0][0] = xi;
    dissip_mat[0][1] = 0.0;
    dissip_mat[1][0] = 0.0;
    dissip_mat[1][1] = 0.0;

    Field::Scalar dissip_coeff(xi);
    Field::Matrix<2, 2> dissip(dissip_mat);
    Field::Matrix<2, 2> coupling(coupling_mat);

    LAL::DiagonalMatrix DissipMatrix;
    LAL::DiagonalMatrix MassMatrix_Disc_Dissip;
    LAL::DiagonalMatrix MassMatrix_Disc;

    CreateMass(fespace_pml, dissip_coeff, DissipMatrix, 1.0);
    CreateMassDG<2>(fespace_pml, dissip, MassMatrix_Disc_Dissip, 1.0);
    CreateMassDG<2>(fespace, density, MassMatrix_Disc, 1.0); // non optimal -> computed everywhere

    LAL::DiagonalMatrix Matrix_p;

    LAL::DiagonalMatrix Matrix_m;

    Matrix_p = MassTotal + 0.5 * dt * DissipMatrix + 0.5 * dt * surface_massMatrix_dissip;

    Matrix_m = MassTotal - 0.5 * dt * DissipMatrix + 0.5 * dt * surface_massMatrix_dissip;

    LAL::DiagonalMatrix MassMatrix_Disc_p;

    LAL::DiagonalMatrix MassMatrix_Disc_m;

    MassMatrix_Disc_p = MassMatrix_Disc + 0.5 * dt * MassMatrix_Disc_Dissip;

    MassMatrix_Disc_m = MassMatrix_Disc - 0.5 * dt * MassMatrix_Disc_Dissip;

    // Compute the source

    LAL::Vector source;

    LAL::Allocate(source, N);

    for (Index i = 0; i < fespace.getNumDoFs(); ++i)
    {
        RealVector P;

        fespace.getDoFCoordinate(i, P);
        RealVector Center = {1.0, 0.0, 0.};
        source(i) = 10e6 * CInfty0<1>(P, Center, 0.1, 1);
    }

    LAL::VectorSequence solution(3);
    solution.Allocate(N);

    LAL::VectorSequence auxilliary_unknown(2);
    auxilliary_unknown.Allocate(2 * M);

    LAL::Vector tmp;
    Index OutputFreq = floor(OutputDeltat / dt);

    Index NOutput = 0;
    std::string FileName = argv[1];

    Mesh mesh = fespace.getMesh();
    std::ofstream stream(FileName + ".mesh");

    print_MESH(mesh, stream);

    stream.close();

    // Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;

    for (Index iStep = 0; iStep * dt < T_end; iStep++)
    {

        // Leap-frog scheme
        MltAddStiffness(fespace, one, -dt * dt, solution.getVector(1), 0.0, solution.getVector(0));

        TransposeMltAddGradient(fespace_pml, -dt * dt, auxilliary_unknown.getVector(1), 1.0, solution.getVector(0));

        solution.getVector(0) -= Matrix_m * solution.getVector(2);
        solution.getVector(0) += MassTotal * 2 * solution.getVector(1);
        solution.getVector(0) += surface_massMatrix * dt * dt * Ricker(iStep * dt, 0.20, 100) * source;

        LAL::Solve(Matrix_p, solution.getVector(0));

        // update the auxiliary variable

        tmp = dt * 0.5 * (solution.getVector(0) + solution.getVector(1));

        Core::MltAdd<decltype(fespace_pml), decltype(fespace_pml), true, true, false,
                     Core::QuadratureOpt::Default,
                     1, // DimU
                     2, // DimV
                     1, // DimI
                     2, // DimJ
                     Field::Matrix<2, 2>,
                     Core::DiffOp::Gradient,
                     Core::DiffOp::Identity>(fespace_pml, fespace_pml, coupling, 1.0, tmp, 0.0, auxilliary_unknown.getVector(0));

        auxilliary_unknown.getVector(0) += MassMatrix_Disc_m * auxilliary_unknown.getVector(1);

        LAL::Solve(MassMatrix_Disc_p, auxilliary_unknown.getVector(0));

        // Writing solution.
        if (iStep % OutputFreq == 0)
        {
            // Writing solution

            ParallelWriteVTK(FileName + "." + std::to_string(NOutput) + ".vtk", fespace, solution.getVector(0));

            ParallelWriteSolutionVIZIR(FileName, fespace, mesh, solution.getVector(0), NOutput);

            NOutput++;
        }

        // Swapping solutions.

        solution.Swap();
        auxilliary_unknown.Swap();
    }

    ParallelWriteVTKTerminate();
    ParallelWriteVIZIRTerminate();

    // Output of the time spent
    time(&end);

    std::cout << "... Leap-Frog iteration end. Elasped time :" << difftime(end, start) << " seconds " << std::endl;
}
