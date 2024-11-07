#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  

int main(int argc, char *argv[])
{
    std::cout << "Testing CubeBndy........................" << std::endl;

    Cube fespace(3,2,1,1,2,3);
    
    for (Index b=0;b<6;++b)
    {
        auto bndy = fespace.getBndyFESpace(b);
        
        std::cout << "Elements of boundary "<<b<<" ................." << std::endl;
        
        for (Index i=0;i<bndy.getNumElements();++i)
        {
            
            Index iElt_cube;
            Index iloc_cube;
            
            bndy._loc_bdny2loc_cube(i,0,iElt_cube,iloc_cube);
            
            std::cout << i << ": " << iElt_cube << std::endl;
        }
        
        std::cout << "DoF of boundary "<<b<<" ......................" << std::endl;
        
        bndy.setAsFESpace();
        
        for (Index i=0;i<bndy.getNumDoFs();++i)
            std::cout << i << ": " << bndy._glob_bdny2glob_cube(i) << std::endl;
    }
    
 
    
    std::cout << "End test................................" << std::endl;
}
 
