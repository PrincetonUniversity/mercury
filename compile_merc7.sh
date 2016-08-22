echo "Compiling mercury7 package"
echo "First for Condor"
condor_compile gfortran -w  -O3 -o mercury7 frag.f90

echo "Now not for Condor"
gfortran -w  -O3 -g -o mercury7_nocondor frag.f90

#echo "Now element7"
#gfortran -w -O3 -o element7 element7.for 

#gfortran -w  -O3 -o mercury7_noclo frag.noclo.f90

