#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;


int main(int argc, char *argv[])
{
    std::cout << "Testing CylinderBndy........................" << std::endl;

    Index FEorderXY = 4;
    Index FEorderZ = 2;
    Index NX_4 = 4;
    Index NZ = 3;

    Cylinder fespace(FEorderXY, FEorderZ, NX_4, NZ);
    
    for (Index b=0;b<3;++b)
    {
        auto bndy = fespace.getBndyFESpace(b);
        
        std::cout << "Elements of boundary "<<b<<" ................." << std::endl;
        
        for (Index i=0;i<bndy.getNumElements();++i)
        {
            
            Index iElt_cylinder;
            Index iloc_cylinder;
            
            bndy._loc_bdny2loc_cylinder(i,0,iElt_cylinder,iloc_cylinder);
            
            std::cout << i << ": " << iElt_cylinder << std::endl;
        }
        
        std::cout << "DoF of boundary "<<b<<" ......................" << std::endl;
        
        bndy.setAsFESpace();
        
        for (Index i=0;i<bndy.getNumDoFs();++i)
            std::cout << i << ": " << bndy._glob_bdny2glob_cylinder(i) << std::endl;
    }
    
 
    
    std::cout << "End test................................" << std::endl;
    
    
    std::cout << "Output of the cylinder boundary" << std::endl;
    
    
    std::string FileName = argv[1];
    std::ofstream ofs(FileName + "_circle.msh");
    
    Circle circle(4,4);
    Mesh circle_mesh = circle.getMesh();
    print_MSH(circle_mesh,ofs);
    ofs.close();
    

    
    for (Index b=0;b<3;++b)
    {
        std::string FileName = argv[1];
        std::ofstream ofs(FileName + "_bndy_" + std::to_string(b) + ".msh");
        
        auto bndy = fespace.getBndyFESpace(b);
        
        Mesh boundary_mesh = bndy.getMesh();
        
        print_MSH(boundary_mesh,ofs);
        
        ofs.close();
    }
    
   
}
