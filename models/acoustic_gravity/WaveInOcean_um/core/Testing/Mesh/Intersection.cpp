#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

#include "OndoMathX.h"

using namespace OndoMathX;
  
int main(int argc, char *argv[])
{
    assert(argc==4);
    
    Mesh background_mesh(argv[1],FormatMsh);
    Mesh mesh_2D_for_intersection(argv[2],FormatMsh);
    
    std::vector<std::vector<Index> > IntersectedElem1;
    std::vector<std::vector<Index> > IntersectedElem2;
    
    Mesh intersected_mesh = meshIntersection(mesh_2D_for_intersection,IntersectedElem1,
                                             background_mesh,         IntersectedElem2);
    
    ofstream output_mesh_stream(argv[3]);
    print_MESH(intersected_mesh,output_mesh_stream);
}
 
