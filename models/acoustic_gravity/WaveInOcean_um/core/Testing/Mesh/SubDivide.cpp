#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  
int main(int argc, char *argv[])
{
    assert(argc==3);
    
    Mesh mesh(argv[1],FormatMsh);

  

  




    ofstream output_mesh_stream(argv[2]);
    print_MESH(mesh,output_mesh_stream);
}
 
