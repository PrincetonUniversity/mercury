echo "Compiling mercury7 package"

gfortran -w -O3 -o element7 element7.for 
gfortran -w  -O3 -o mercury7 frag.f90
#gfortran -w  -O3 -o mercury7_noclo frag.noclo.f90

