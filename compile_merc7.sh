echo "Compiling mercury7 package"

gfortran -w -O3 -o element7 element7.for 
gfortran -w -O3 -o mercury7 frag.f90

