echo "Compiling mercury7 package"
echo "First for Condor"
condor_compile gfortran -w  -O3 -o mercury7_multcollisioncounter frag.check_number_of_multiple_collisions.f90

echo "Now not for Condor"
gfortran -w  -O3 -g -o mercury7_multcollisioncounter_nocondor frag.check_number_of_multiple_collisions.f90

#echo "Now element7"
#gfortran -w -O3 -o element7 element7.for 

#gfortran -w  -O3 -o mercury7_noclo frag.noclo.f90

#gfortran -w   -g -o mercury7_nocondor frag.f90
#gfortran -w   -g -o mercury7_lotsofnames frag.lotsofnewnames.f90
