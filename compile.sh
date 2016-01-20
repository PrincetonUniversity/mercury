echo "Compiling mercury6 package"

gfortran -w -o element6 element6.for 
gfortran -w -o close6 close6.for 
gfortran -w -o mercury6 mercury6_2.for 
#gfortran trying_other.for -o tryingother
#gfortran mercury6_2.noce.for -o mercury6_noce

#echo "Create sample data files? (yes/no)"
#read create_samples

#if [ $create_samples = 'yes' ] ; then
#  for i in *.sample; do
#    cp $i `echo $i | sed 's/.sample//g'`
#  done
#fi