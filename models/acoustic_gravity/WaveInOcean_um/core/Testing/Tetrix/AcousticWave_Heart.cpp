#include <iostream>
#include <fstream>
#include <ctime>
#include <string>
#include <algorithm>

#include "OndoMathX.h"

using namespace OndoMathX;
  
 

int main(int argc, char *argv[])
{
    
    {
        vector<Real> Vec(15);
        std::srand(unsigned(std::time(nullptr)));
        std::generate(Vec.begin(), Vec.end(), std::rand);
        Real Buffer_Grad_U[42];

        time_t start, end;

        time(&start);
        for (Index i=0;i<20000000;++i)
        {
              TetriX::Tet_15_14::GradientInterpolation(Vec.data(),Buffer_Grad_U);
            
        }
        time(&end);
    
        std::cout << "Elapsed time Non Optimized Gradient Interpolation :" << difftime(end, start) << " seconds " << std::endl;
        

        time(&start);
        for (Index i=0;i<20000000;++i)
        {
            TetriX::Tet_15_14::GradientInterpolation_Opt(Vec.data(),Buffer_Grad_U);
        }
        time(&end);
    
        std::cout << "Elapsed time Optimized Gradient Interpolation :" << difftime(end, start) << " seconds " << std::endl;

    }
 


   
    std::cout << "Running AcousticWave_Heart........................" << std::endl;
    
    std::cout << "Simulation of an scalar acoustic wave in a heart geometry" << std::endl;
    std::cout << "A number of " << omp_get_max_threads() << " threads are used by OpenMP" << std::endl;
    
    Mesh mesh_with_surface(argv[1],FormatMesh);

    std::cout << "N. vertices of the mesh \t \t \t" << mesh_with_surface.getNumVertices() << std::endl; 
    std::cout << "N. faces of the mesh \t \t \t \t" <<  mesh_with_surface.getNumFaces() << std::endl; 
    std::cout << "N. edges of the mesh \t \t \t \t" <<  mesh_with_surface.getNumEdges() << std::endl; 

    std::cout << "N. Elements of the mesh \t \t \t" << mesh_with_surface.getNumElements() << std::endl;  
  
    Mesh mesh = mesh_with_surface.extract(mesh_with_surface.getElementList(3));
 

    std::vector<Index> coloring;
    Index numColors;
    
    getColoring(mesh,coloring,numColors);
    mesh.setLabels(coloring);
    ofstream output_mesh_stream(std::string(argv[2])+ "_coloring.msh");
    print_MSH(mesh,output_mesh_stream);

    TetriXP2Space FESpace(mesh);

    std::cout << "N. of tetra \t" << FESpace.getNumElements() << std::endl; 

    std::cout << "N. Dof of P2-enriched Finite element space \t" << FESpace.getNumDoFs() << std::endl; 

    std::cout << "N. Dof of P2 Finite element space \t \t" << FESpace.getP2NumDoFs() << std::endl; 

    std::cout << "N.colors \t \t" <<  FESpace.getNumColors() << std::endl; 

    for (Index i = 0; i < FESpace.getNumColors()  ; ++i)
        std::cout << FESpace.getNumElements(i) << " "; 
    std::cout << "\n";

    Index N = FESpace.getNumDoFs();

    LAL::Vector source;
    LAL::Allocate(source, N);

    for (Index i=0;i<N;++i)
    {
        RealVector xyz;
        
        FESpace.getDoFCoordinate(i,xyz);
        
        Real r = sqrt((xyz[0]-20)*(xyz[0]-20) +  (xyz[1]-70)*(xyz[1]-70) + (xyz[2]-70)*(xyz[2]-70));

        source[i] = Gaussian(r,0,0.004);

       // source[i] = xyz[0];

    }

    std::string FileName = argv[2];

    CreateOutputDirectory(FileName);
    WriteVTK(FileName + "_Source.vtk", FESpace, source);


    




/*
    for (Index n=0;n < FESpace.getNumElements();++n)
    {
         RealMatrix3x3 F;
         FESpace.getGradDef(n,0,F);

         cout << n << " : Jacobian " << ArrayAlgebra::Det(F) << "\n";
    }

*/
    std::ofstream file_mesh_tet(FileName+"_tet.mesh");
    print_MESH(mesh,file_mesh_tet);
    file_mesh_tet.close();

    WriteSolutionVIZIR(FileName+"_Source", FESpace, mesh, source);

    LAL::DiagonalMatrix massMatrix;

    LAL::Allocate(massMatrix,N,N);
 
    FESpace.AssembleMass(massMatrix,1.0);

    Index N_iter;
    Real Residue;

    std::vector<Real> Buffer_Sol_U(15,0.0);
    std::vector<Real> Buffer_Sol_V(15,0.0);
    std::vector<Real> Buffer_Grad_U(42,0.0);
    std::vector<Real> Buffer_Grad_V(42,0.0);


/*
    cout << "test with the function x" << endl;
    for (Index i = 0; i < 15; i++) 
    {
        Buffer_Sol_U[i] =  TetriX::Tet_15::Point[i][0];
    }

    TetriX::Tet_15_14::GradientInterpolation(Buffer_Sol_U,Buffer_Grad_U);
      
    for (Index i = 0; i < 42; i++)
        cout << Buffer_Grad_U[i] << endl;
 
    cout << "test with the function y" << endl;
    for (Index i = 0; i < 15; i++) 
    {
        Buffer_Sol_U[i] =  TetriX::Tet_15::Point[i][1];
    }

    TetriX::Tet_15_14::GradientInterpolation(Buffer_Sol_U,Buffer_Grad_U);
      
    for (Index i = 0; i < 42; i++)
        cout << Buffer_Grad_U[i] << endl;

    cout << "test with the function z" << endl;
    for (Index i = 0; i < 15; i++) 
    {
        Buffer_Sol_U[i] =  TetriX::Tet_15::Point[i][2];
    }

    TetriX::Tet_15_14::GradientInterpolation(Buffer_Sol_U,Buffer_Grad_U);
      
    for (Index i = 0; i < 42; i++)
        cout << Buffer_Grad_U[i] << endl;



    //Test de TransposeGradientInterpolation

    for (Index i = 0; i < 15; i++) Buffer_Sol_U[i] =  i;
    for (Index i = 0; i < 42; i++) Buffer_Grad_U[i] =  i/3.0;

    Real sp1 = 0;
    TetriX::Tet_15_14::GradientInterpolation(Buffer_Sol_U,Buffer_Grad_V);
    for (Index i = 0; i < 42; i++) sp1+= Buffer_Grad_V[i]*Buffer_Grad_U[i];

    Real sp2 = 0;
    TetriX::Tet_15_14::TransposeGradientInterpolation(Buffer_Grad_U,Buffer_Sol_V);
    for (Index i = 0; i < 15; i++) sp2+= Buffer_Sol_U[i]*Buffer_Sol_V[i];

    cout << "Test transposé " << sp1 << " " << sp2 << endl;

*/
 

    Real rho = LAL::PowerIteration([&massMatrix](LAL::Vector &v) {LAL::Solve(massMatrix,v);},
                              [&FESpace](LAL::Vector &u,LAL::Vector &v) {FESpace.MltStiffness(u,v,1.0);},
                              N,N_iter,Residue);
    

    cout << rho << " " << N_iter << endl; 



    Real dt =  0.95 * (2.0 / sqrt(rho));

    cout << "Time step " << dt << endl;


    LAL::VectorSequence solution(3);
    solution.Allocate(N);
    
    Real T_end = 200.0;
    Real OutputDeltat = 1;
    Index OutputFreq = floor(OutputDeltat/dt);
    Index NOutput = 0;
 

    //Set the time variable
    time_t start, end;
    time(&start);

    cout << "Leap-frog iterations starts..." << std::endl;
    for (Index iStep = 0; iStep*dt < T_end; iStep++)
    {
       //Leap-frog scheme
      
        FESpace.MltStiffness(solution.getVector(1),solution.getVector(0),-dt*dt);
        LAL::Solve(massMatrix,solution.getVector(0));

        solution.getVector(0) -= solution.getVector(2); 
        solution.getVector(0) += 2*solution.getVector(1);
        solution.getVector(0) += dt*dt*Ricker(iStep*dt,0.20,1)*source;       
  
        // Writing solution.
        if (iStep % OutputFreq == 0 && iStep>0)
        {
            cout << "Output " << NOutput << endl;
            // Writing solution
            WriteVTK(FileName+"."+std::to_string(NOutput) + ".vtk", FESpace,solution.getVector(0));
             
            WriteSolutionVIZIR(FileName, FESpace, mesh, solution.getVector(0), NOutput);

            NOutput++;
        }
  
        // Swapping solutions.
        solution.Swap();
    }


}
 
