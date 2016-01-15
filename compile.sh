echo "Compiling mercury6 package"

gfortran element6.for -o element6
gfortran close6.for -o close6
gfortran mercury6_2.for -o mercury6
#gfortran trying_other.for -o tryingother
#gfortran mercury6_2.noce.for -o mercury6_noce

#echo "Create sample data files? (yes/no)"
#read create_samples

#if [ $create_samples = 'yes' ] ; then
#  for i in *.sample; do
#    cp $i `echo $i | sed 's/.sample//g'`
#  done
#fi