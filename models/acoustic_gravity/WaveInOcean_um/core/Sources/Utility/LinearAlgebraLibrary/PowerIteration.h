#pragma once
 
 
namespace OndoMathX
{

namespace LAL
{
    inline Real PowerIteration(std::function<void(Vector&)> Solve,
                             std::function<void(Vector&,Vector&)> MltAdd,
                             Index N,
                             Index & IP_NIter, Real & IP_Residue,
                             Real Tol = 10e-7, Index NItermax = 200)
        {
            // Algorithm variables.
            Real Lambda0, Lambda1;
       
            // Iteration counters.
            IP_NIter = 0;
            
            // Creating relevant vectors.
            Vector FEVect1, FEVect2;
            
            // Allocating vectors.
            Allocate(FEVect1, N);
            Allocate(FEVect2, N);
            
            Real * dataV1 = getData(FEVect1);
            
            // Initializing a random like vector.
            for (Index i = 0; i < N; i++)
                dataV1[i] = i % 7 + 3*i % 3 + (i % 2)*1e-6 + i % 20;
            
            // Initializing maximum eigenvalue.
            Lambda0 = 0.0;
            Lambda1 = 0.0;
            
            do{
                // Swapping values of eigen value.
                Lambda0 = Lambda1;
                
                // Applying operator appearing on the right-hand side of the generalized eigen value problem.
                MltAdd(FEVect1, FEVect2);
            
                // Solving linear system.
                Solve(FEVect2);
                
                FEVect1 = FEVect2;
                
                // Extracting infinite norm.
                Lambda1 = GetInfiniteNorm(FEVect1);
                
                // Normalizing.
                Scale(1.0 / Lambda1, FEVect1);
                
                // Computing new tolerance.
                if (Lambda0 != 0.0) IP_Residue = fabs(Lambda1 - Lambda0) / fabs(Lambda0);
                else IP_Residue = std::numeric_limits<double>::max();
                
                // Updating number of iterations.
                IP_NIter++;
                
            } while (IP_Residue > Tol && IP_NIter <= NItermax);
            
        
            return Lambda1;
        }
        //---------------------------------------------------------------------------------//
        
 
} //LAL

} // OndoMathX
 
