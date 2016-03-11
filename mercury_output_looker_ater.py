#Created 03102016 by Joshua Wallace
# This works closely with mercury_outputreader.py in examining the output of
# the mercury7 stuff.
# This one deals more with post-processing, than collecting the data
#
#AGBTG

import numpy as np
import collections

import mercury_outputreader as outputreader


def in_order(the_list):
    for i in range(len(the_list)-1):
        if the_list[i+1] < the_list[i]:
            return False
    else:
        return True


def mutual_hill_radii_checker(planets,central_object_mass=1.):
    """This will calculate the mutual hill radii of the final planets 
    (which are handed as an aei_info object) and compare to 
    their separations.  The central object mass is an optional parameter
    with default value of 1 in solar masses.
    Returned is a list of delta a/mutual hill radii.
    """
    if not isinstance(planets, outputreader.aei_singletime):#If not an aei_info object
        raise TypeError("I wasn't given an aei_singletime instance!")

#    if not (len(planets) > 1):
#        raise TypeError("The thing passed to me doesn't have more than one object in it.")
    if not isinstance(planets.a,collections.Sequence):
        raise TypeError("The thing passed to me doesn't have more than one object in it.")


    #Check that objects are in order
    if not in_order(planets.a):
        print "Planets weren't in order"
        indexes = range(len(planets.a))
        indexes.sort(key=planets.a.__getitem__)

        aes = map(planets.a.__getitem__, indexes)
        masses = map(planets.mass.__getitem__, indexes)

        print aes
        print masses

    else:
        aes = planets.a
        masses = planets.mass

    output = []

    for i in range(len(planets.a)-1):
        delta_a = aes[i+1] - aes[i]
        mutual_hill_radius = 0.5 * (aes[i+1] + aes[i]) * np.power( (masses[i+1] + masses[i])/ (3.0*central_object_mass) , 1./3.)
        output.append(delta_a/mutual_hill_radius)

    return output
