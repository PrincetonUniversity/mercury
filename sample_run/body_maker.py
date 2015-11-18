import numpy as np
import numpy.random as rand
import sys
import scipy.interpolate as intp


rand.seed(1844)
e_0 = 0.01
i_0 = 0.5*np.pi/180.
mass = 0.0167 * 0.000003003

if len(sys.argv) != 2:
    print "I need one argument to this code!  The number of protoplanets"
    quit()

number = int(sys.argv[1])

if number != 150:
    print "You aren't giving me an input of 150"
    print "If you're sure about this, and want to proceed, change the code"
    print "Now quitting..."
    quit()

def sigma_1(a):
    if a < 0.3 or a > 2.0:
        return 0.
    if a < 0.7:
        return 34.1494 * a - 10.2448
    else:
        return 8.*np.power(a,-1.5)

a_distribution = np.linspace(.30000000000001,2.,4251)
sigma_distribution = []
sigmas = []
cumulative = 0.
for i in range(len(a_distribution)):
    cumulative+=sigma_1(a_distribution[i])
    sigmas.append(sigma_1(a_distribution[i]))
    sigma_distribution.append(cumulative)

print sigmas[0]
print sigmas[1]
print a_distribution[np.argmax(sigmas)]
    
sigma_distribution = np.subtract(sigma_distribution,sigma_distribution[0])
sigma_distribution = np.divide(sigma_distribution,sigma_distribution[-1])

a_func_sigma = intp.interp1d(sigma_distribution,a_distribution)
a_to_start_out = a_func_sigma(np.linspace(0.,1.,number))

#s = 'Initial_aes/initial_a_sim1_' + str(number) + '.txt'
#np.savetxt(s,np.transpose(a_to_start_out))

#Now, for the rest of the oribtal elements:

eccentricity = rand.random(number)
eccentricity = np.multiply(eccentricity,e_0)
inclination = rand.random(number)
inclination = np.multiply(inclination,i_0)
capital_omega = rand.random(number)
capital_omega = np.multiply(capital_omega,2.*np.pi)
omega = rand.random(number)
omega = np.multiply(omega,2.*np.pi)
mean_anomaly = rand.random(number)
mean_anomaly = np.multiply(mean_anomaly,2.*np.pi)

particle_number = 2

f = open('big.in','w')

f.write(')O+_06 Big-body initial data\n')
f.write(')\n')
f.write(')----------------------------------------\n')
f.write(' style (Cartesian, Asteroidal, Cometary) = Asteroidal\n')
f.write(' epoch (in days) = 0.\n')
f.write(')-----------------------------------------\n')

for i in range(number):
    f.write(str(particle_number) + '    m=' + '%E' % mass + ' r=5.0D0' + ' d=3.0' + '\n')
    f.write(' %E  %E  %E\n' % (a_to_start_out[i], eccentricity[i], inclination[i]))
    f.write(' %E  %E  %E\n' % (omega[i], capital_omega[i], mean_anomaly[i]))
    f.write(' 0. 0. 0.\n')
    particle_number += 1

#Jupiter_inclination = -2.1198273483
#Saturn_minus_Jupiter_inclination = 2.48446 - 1.30530

#Now, add in Jupiter and Saturn
f.write('Jupiter' + '    m=' + '9.5440e-4' + ' r=3.0D0' + ' d=1.33' + '\n')
f.write(' %E  %E  %E\n' % (5.2053529587525267, 0.043367030904308053, 0.406143799218))
f.write(' %E  %E  %E\n' % (303.968060227, 202.09530903, 359.970599313))
f.write(' 0. 0. 0.\n')

f.write('Saturn' + '    m=' + '2.8575e-4' + ' r=3.0D0' + ' d=0.70' + '\n')
f.write(' %E  %E  %E\n' % (9.53707032, 0.05415060, 1.20425544538))
f.write(' %E  %E  %E\n' % (136.436954896, 177.62498702, 359.520222278))
f.write(' 0. 0. 0.\n')
f.close()
