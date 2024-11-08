#Exe#################################

add_executable (WaveInOcean_um WaveInOcean_um.cpp)

####################################

#Launcher#################################

configure_file (WaveInOcean_um_launcher.sh.in WaveInOcean_um_launcher.sh)

####################################
