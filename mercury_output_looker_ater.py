#Created 03102016 by Joshua Wallace
# This works closely with mercury_outputreader.py in examining the output of
# the mercury7 stuff.
# This one deals more with post-processing, than collecting the data
#
#AGBTG

import numpy as np
import matplotlib.pyplot as pp
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

    print "  received " + str(len(planets.a)) + " planets"

    #Check that objects are in order
    if not in_order(planets.a):
        print "Planets weren't in order"
        indexes = range(len(planets.a))
        indexes.sort(key=planets.a.__getitem__)

        aes = map(planets.a.__getitem__, indexes)
        masses = map(planets.mass.__getitem__, indexes)
        ees    = map(planets.e.__getitem__, indexes)

        print aes
        print masses
        print ees

    else:
        aes = planets.a
        masses = planets.mass

    output = []

    for i in range(len(planets.a)-1):
        print i
        delta_a = aes[i+1] - aes[i]
        mutual_hill_radius = 0.5 * (aes[i+1] + aes[i]) * np.power( (masses[i+1] + masses[i])/ (3.0*central_object_mass) , 1./3.)
        output.append(delta_a/mutual_hill_radius)

    return output


def plot_mutual_hill_radii(planets,central_object_mass=1.):
    """This will calculate the mutual hill radii of the final planets 
    (which are handed as an aei_info object) and compare to 
    their separations.  The central object mass is an optional parameter
    with default value of 1 in solar masses.  This information is then 
    plotted and the function returns a figure
    """

    if not isinstance(planets,outputreader.aei_singletime):
        raise TypeError("plot_mutual_hill_radii was not given an aei_singletime instance!")

    fig = pp.figure()

    stuff = mutual_hill_radii_checker(planets)

    a_list = []
    for i in range(len(planets.a)-1): #This will allow us to plot separation as a function of position
        a_list.append( (planets.a[i] + planets.a[i+1])/2.0)

    pp.scatter(a_list,stuff)
    pp.xlabel("Semi-major axis (AU), defined to be between the planets")
    pp.ylabel("Separation / mutual Hill radii")
    pp.ylim(bottom=0)
    pp.xlim(left=0)

    return fig



def plot_a_func_time(aes, times, which_are_final_bodies=None,year_unit='kyr'):
    """This function takes a list of each individual body's semimajor axes,
    which itself can be stored as a list, as well as a list of lists of the 
    times corresponding to semimajor axes passed with the first argument
    and then plots semi-major axis as a function of time for the objects 
    passed to me.
    which_are_final_bodies: pass the index of the final bodies if you want 
    those lines plotted thicker
    """  
    year_unit_dict  = {"Myr":1.e6,"kyr":1.e3}

    fig = pp.figure()

    for i in range(len(aes)):
        pp.plot(np.divide(times[i],year_unit_dict[year_unit]),aes[i],color='blue',linewidth=0.5)

    if which_are_final_bodies != None:
        for i in range(len(which_are_final_bodies)):
            pp.plot(np.divide(times[which_are_final_bodies[i]],year_unit_dict[year_unit]),aes[which_are_final_bodies[i]],color='red')#,color='blue',linewidth=1.5)

    pp.xscale(u'log')

    pp.xlabel("Time ("+year_unit+")")
    pp.ylabel("Semimajor axis (AU)")


    return fig
