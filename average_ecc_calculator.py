import numpy as np
import sys
import os
def tail(f, n, offset=0): # Taken from http://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail
  stdin,stdout = os.popen2("tail -n "+n+offset+" "+f)
  stdin.close()
  lines = stdout.readlines(); stdout.close()
  return lines[:,-offset]

path = "./"


with open("final_bodies.txt","r") as f:
    names_of_final_bodies = f.readlines()

print "Number of final bodies: " + len(names_of_final_bodies)

final_bodies_info = []

for i in range(len(names_of_final_bodies)):
    with open(names_of_final_bodies[i] + ".aei","r"):
        last_line = tail(f, 1)
    last_line_split = last_line.split()

    if (float(last_line_split[0]) - 3e5) > 0.0001:
        raise RuntimeError("Times don't match up, this must not be a final body")

    e = float(last_line_split[2])
    x = float(last_line_split[5])
    y = float(last_line_split[6])
    z = float(last_line_split[7])
    final_bodies_info.append( {'e':float(e), 'xyz': (float(x),float(y),float(z))})

total_e = 0.
total_num = 0

for i in range(len(final_bodies_info)):
    dist_squared = final_bodies_info[i]['xyz'][0]**2 +  final_bodies_info[i]['xyz'][1]**2 +  final_bodies_info[i]['xyz'][2]**2 
    if np.sqrt(dist_squared) > 0.008877:
        total_num += 1
        total_e += final_bodies_info[i]['e']

print "Got " + str(total_num) + " bodies outside Roche radius"
print "Average ecc. is: " + str(total_e/float(total_num) )
