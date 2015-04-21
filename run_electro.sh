#! /bin/sh

# a small shell script to clean up the files in the folder and run the 
# immersed boundary code. Be sure to run this in its own folder. It will remove all 
# data and .png files in order to "start fresh" 


echo 'removing folders and .png files'

rm *.png
rm -r forces_valveless_test markers_valveless_test particles_valveless_test velu_valveless_test velv_valveless_test vorticity_valveless_test ave_vel_test

echo 'compiling code (killdevil and Kure)'

icc electro_heart_tube_peri_ciona.c -I/nas02/apps/fftw-3.2.2/include -L/nas02/apps/fftw-3.2.2/lib -lfftw3

#echo 'running code now'

#bsub -o out.%J ./a.out
