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



def plot_a_func_time(aes, times, which_are_final_bodies=None,year_unit='kyr',title=None):
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
            #print "   final body plotting as red: " + str(i)
            pp.plot(np.divide(times[which_are_final_bodies[i]],year_unit_dict[year_unit]),aes[which_are_final_bodies[i]],color='red')#,color='blue',linewidth=1.5)

    pp.xscale(u'log')
    
    if title != None:
        pp.title(title)

    pp.xlabel("Time ("+year_unit+")")
    pp.ylabel("Semimajor axis (AU)")


    return fig


def adhoc_number_insideoutside_roche(time_list,aei_functime_list,num_objects_list):
    """This function does an ad hoc calculation of the number
    of bodies inside and outside the Roche radius by simply comparing
    the semi-major axes against the Roche radius.  In the case of 
    zero eccentricity, this would provide the correct answers, but with 
    eccentricity an actual calculation would need to be made based on 
    the actual position of the bodies.
    Input: list of times, list of aeis function of time, list of total number
       Note: the output of aei_func_time will
        work here
    Output: figure object
    """


    if not all(isinstance(item,outputreader.aei_singletime) for item in aei_functime_list):
        raise TypeError("I was not passed a list entirely full of aei_singletime objects!")

    if len(time_list) != len(aei_functime_list):
        raise TypeError("The time list and aei list are not the same length!")

    if len(time_list) != len(num_objects_list):
        raise TypeError("The time list and aei list are not the same length!")

    num_inside = []
    num_outside = []

    for i in range(len(aei_functime_list)):
        n_in = 0
        n_out = 0
        for value in aei_functime_list[i].a:
            if value < 0.00888: #Hardcoded value for roche radius
                n_in += 1
            else:
                n_out += 1


        num_inside.append(n_in)
        num_outside.append(n_out)

    temp = np.add(num_inside,num_outside)
    for i in range(len(temp)):
        if temp[i] != num_objects_list[i]:
            raise TypeError("The sum of bodies inside and outside do not equal the total number of bodies! " + str(i) )

    fig = pp.figure()

    time_list = np.multiply(time_list,1e-3)
    unit = 'kyr'

    pp.step(time_list,num_objects_list,where='post',lw=2.5,label='all bodies')
    pp.step(time_list,num_outside,where='post',lw=1.8,label='outside Roche')
    pp.step(time_list,num_inside,where='post',lw=1.8,label='inside Roche')

    pp.legend(loc='best')

    pp.xlabel("Time (" + unit + ")")
    pp.xscale('log')
    pp.ylabel("Number of bodies")

    return fig
