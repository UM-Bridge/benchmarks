#Exe#################################

add_executable (Rasterization Mesh/Rasterization.cpp)
add_executable (Intersection Mesh/Intersection.cpp)
add_executable (Coloring Mesh/Coloring.cpp)
add_executable (SubDivide Mesh/SubDivide.cpp)

add_executable (AcousticWave_2D_Mesh FESpace/AcousticWave_2D_Mesh.cpp)
add_executable (AcousticWave_Square FESpace/AcousticWave_Square.cpp)
add_executable (AcousticWave_Square_PML FESpace/AcousticWave_Square_PML.cpp)
add_executable (GravitoAcousticWave_Square_PML FESpace/GravitoAcousticWave_Square_PML.cpp)
add_executable (ParaxialWave_Square FESpace/ParaxialWave_Square.cpp)
add_executable (ParaxialWave45_Square FESpace/ParaxialWave45_Square.cpp)
add_executable (Transport_Square FESpace/Transport_Square.cpp)
add_executable (AcousticWave_Circle FESpace/AcousticWave_Circle.cpp)
add_executable (AcousticWave_Cube FESpace/AcousticWave_Cube.cpp)
add_executable (AcousticWave_Cylinder FESpace/AcousticWave_Cylinder.cpp)
add_executable (ElasticWave_Plate FESpace/ElasticWave_Plate.cpp)
add_executable (IncompressibleElasticWave_Square FESpace/IncompressibleElasticWave_Square.cpp)

add_executable (CubeBndy Other/CubeBndy.cpp)
add_executable (CylinderBndy Other/CylinderBndy.cpp)
add_executable (HardCodedGLE Other/HardCodedGLE.cpp)

add_executable (NLElasticWave_Cube FESpace/NonLinear/NLElasticWave_Cube.cpp)
add_executable (NLElasticWave_Bar FESpace/NonLinear/NLElasticWave_Bar.cpp)

add_executable (AcousticWave_Heart Tetrix/AcousticWave_Heart.cpp)

####################################

#Launcher#################################

configure_file (Mesh/MeshTest.sh.in MeshTest.sh)
configure_file (FESpace/FESpaceTest.sh.in FESpaceTest.sh)
configure_file (Other/OtherTest.sh.in OtherTest.sh)
configure_file (Tetrix/Tetrix.sh.in Tetrix.sh)

####################################
