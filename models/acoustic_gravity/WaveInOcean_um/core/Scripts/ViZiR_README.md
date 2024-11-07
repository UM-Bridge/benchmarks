# ViZiR Readme

In order to generate files needed by ViZiR to deal with transient simulations, we can use the `generate_vizir_movie.py` script.

It requires two input arguments and assumes that you are inside the `Scripts` folder: 
- `--solution` which corresponds to the prefix that all solutions share. For instance it can be `AcousticWave_Circle`
- `--result_dir` the relative path to the result directory which contains the results. It assumes that there is only mesh in this folder and that the path is relative to the `Scripts` folder.

In order to make sure that there is only one mesh, use the `CreateOutputDirectory()` function in the `main.cpp` of your computation. Refer to the function's documentation to see what is expected.

For instance: 
```c++ 
std::string FileName = argv[1];
CreateOutputDirectory(FileName);
/*
    ...
*/
WriteSolutionVIZIR(FileName, fespace, mesh, solution.getVector(0), NOutput);
```

Here is an example of usage for the script: 
```bash
cd Scripts 
python3 generate_vizir_movie.py --solution AcousticWave_Circle --result_dir ../Results/AcousticWaveCircle
```
This should generate binary solutions inside `../Results/AcousticWaveCircle` that correspond to the surfacic elements of the mesh, as well as a binary version of the mesh with surface elements. 

The script will also generate a `vizir.movie` file and `AcousticWave_Circle.sols` file inside the results folder, as well as an empty `VizirOut` directory. 

The `.sols` file can be used to visualize the transient solution: 
```bash
cd ../Results/AcousticWaveCircle
vizir4 -in AcousticWave_Circle.surfOK.meshb -sols AcousticWave_Circle.sols
```
Inside ViZiR, using the `Utilities/Non Stationary Solutions Menu` we can see the dynamics of the solution. 
Once the visualisation has been set up (ie showing the solution on the mesh, setting up the camera angle, color palette, ...), press `w` inside ViZiR to write the `vizir.state` file. This will save the current view, which is the one that will be used to generate the screenshots of the simulation.

Finally, we can generate screenshots of each timestep using the `.movie`  and the `.state` files:
```bash
vizir4 -movie -state vizir.state   
```
This will generate the screenshots inside `Results/AcousticWaveCircle/VizirOut`.

Then you can generate a movie using the images, for instances with `ffmpeg` which is a nice open source tool.
It can be installed simply through Homebrew:
```bash
brew install ffmpeg
```

Hereafter is an example of how to use it to generate a movie:
```bash
cd VizirOut
ffmpeg -y -framerate 5 -i AcousticWave_Circle.%05d.jpg ./AcousticWaveCircle.mp4
```