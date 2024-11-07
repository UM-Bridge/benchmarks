#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  
int main(int argc, char *argv[])
{
    assert(argc==3);
    
    Mesh surface(argv[1],FormatMsh);

    auto rasterization = CreateRasterization(surface,8e5, 5e5, 60, 60, 2.5e3, 2.5e3);
    Mesh rasterized_mesh = CreatMesh(rasterization);
    
    ofstream output_mesh_stream(argv[2]);
    print_MESH(rasterized_mesh,output_mesh_stream);
}
 
