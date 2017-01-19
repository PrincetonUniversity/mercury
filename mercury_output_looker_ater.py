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
    Returned is a tuple: first element, list of delta a/mutual hill radii.
       Second element, list of "mutual semimajor axis", the average of the 
       semimajoraxes of the two bodies.
       Third element, list of delta a'/mutual hill radii, where delta a' is the
       separation between apoastron of inner object and periastron of outer object.
       Fourth element, the average of the apoastron of inner object and 
       periastron of outer object
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

        print "aes, masses, and ees:"
        print aes
        print masses
        print ees

    else:
        aes = planets.a
        masses = planets.mass

        print "aes, masses, and ees:"
        print aes
        print masses
        print ees

    output1 = []
    output2 = []
    output3 = []
    output4 = []

    for i in range(len(planets.a)-1):
        print i
        delta_a = aes[i+1] - aes[i]
        delta_aprime = aes[i+1]*(1-ees[i+1]) - aes[i]*(1+ees[i])
        mutual_hill_radius = 0.5 * (aes[i+1] + aes[i]) * np.power( (masses[i+1] + masses[i])/ (3.0*central_object_mass) , 1./3.)
        output1.append(delta_a/mutual_hill_radius)
        output2.append( (aes[i+1] + aes[i])/2.0 )
        output3.append(delta_aprime/mutual_hill_radius)
        output4.append( ( aes[i+1]*(1-ees[i+1]) + aes[i]*(1+ees[i]) )/2.0 )

    return (output1, output2, output3, output4)


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

    (stuff1, stuff2) = mutual_hill_radii_checker(planets)

    pp.plot( (0.008877,0.008877),(0,50),linestyle='--',color='black')

    pp.scatter(stuff2,stuff1)
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
    pp.step(time_list,num_inside,where='post',lw=1.8,label='inside Roche',color='red')
    pp.step(time_list,num_outside,where='post',lw=1.8,label='outside Roche',color='green')

    np.savetxt("inout_adhoc_roche_output.txt",(time_list,num_objects_list,num_inside,num_outside) )

    pp.legend(loc='best')

    pp.xlabel("Time (" + unit + ")")
    pp.xscale('log')
    pp.ylabel("Number of bodies")

    return fig

def plot_collision_scatterplot_divide_by_radius(collision_info,radialbins,title=None):
    """This function will take the output of 
    outputreader.collision_info_extractor through the input
    collision_info and a range of radial bins (in AU) through 
    radial bins (also a list of titles for each of the plots)

    and then output a list of figures to plot"""

    if  len(radialbins) < 2:
        raise RuntimeError("radialbins is not long enough to even be a bin range")

    collisions = collision_info[0]

    collisions_binned_by_radius = []
    
    if not (all(collisions[i].radius != None for i in collisions) ):
        raise RuntimeError("At least one of the collisions doesn't have a radius!")

    for i in range(len(radialbins)-1):
        temp = []
        lower_r = radialbins[i]
        upper_r = radialbins[i+1]
        for j in range(len(collisions)):
            if collisions[j].radius >= lower_r and collisions[j].radius < upper_r:
                temp.append(collisions[j])

        collisions_binned_by_radius.append(temp)


    to_return = []
    if title==None:
        for i in range(len(collisions_binned_by_radius)):
            title.append(plot_collision_scatterplot(filename="literallydoesntmatter",
                                                    collision_info = (collisions_binned_by_radius[i],
                                                                      collision_info[1],collision_info[2])))

    else:
        if len(title) != len(collisions_binned_by_radius):
            raise RuntimeError("don't have enough titles for all the plots!")
        for i in range(len(collisions_binned_by_radius)):
            title.append(plot_collision_scatterplot(filename="literallydoesntmatter",
                                                    collision_info = (collisions_binned_by_radius[i],
                                                                      collision_info[1],collision_info[2]),
                                                    title=title[i]))

    return to_return


def collision_type_func_time(input_file="info.out"):
    """Takes a specified info.out file and returns a tuple,
    first element is times of collisions, second is list of the
    collision types (collision_info instances), third is number of bodies
    after collision
    """

    collision, central, eject = outputreader.collision_info_extractor(input_file)

    time_list = []
    collision_classification_list = []
    numbers_list = []
    for i in range(len(collision)):
        time_list.append(collision[i].time)
        collision_classification_list.append(collision[i].classification)
        numbers_list.append(len(collision[i].remnant_names))
    

    return (time_list,collision_classification_list,numbers_list)


def plot_collision_type_func_time(input_file="info.out",x_limits=None):
    """Given an info.out file, this figures out the collision types as 
    a function of time and then plots this
    and returns a pp.figure to plot.
    """

    fig = pp.figure()
    time, collisions, numbers = collision_type_func_time(input_file=input_file)


    time_merger = [time[i] for i in range(len(time)) if collisions[i] in [outputreader.collision_type.SIMPLE_MERGER,outputreader.collision_type.EFFECTIVE_MERGER,outputreader.collision_type.GRAZE_MERGER] ]
    #merger = [item for item in collisions if item in [outputreader.collision_type.SIMPLE_MERGER,outputreader.collision_type.EFFECTIVE_MERGER,outputreader.collision_type.GRAZE_MERGER] ]


    time_hitandrun = [time[i] for i in range(len(time)) if collisions[i] in [outputreader.collision_type.HIT_AND_RUN] ]
    #hitandrun = [item for item in collisions if item in [outputreader.collision_type.HIT_AND_RUN] ]

    time_fragment_increasenum = [time[i] for i in range(len(time)) if (collisions[i] in [outputreader.collision_type.NONGRAZING_FRAG,outputreader.collision_type.GRAZING_FRAG]  and numbers[i] > 2 )  ]
    time_fragment_numstaysame = [time[i] for i in range(len(time)) if (collisions[i] in [outputreader.collision_type.NONGRAZING_FRAG,outputreader.collision_type.GRAZING_FRAG]  and numbers[i] == 2 )  ]
    #fragment =  [item for item in collisions if item in [outputreader.collision_type.NONGRAZING_FRAG,outputreader.collision_type.GRAZING_FRAG] ]


    pp.scatter(time_merger,[3]*len(time_merger),color='red',alpha=0.4,s=10,label="merge")
    pp.scatter(time_hitandrun,[2]*len(time_hitandrun),color='blue',alpha=0.4,s=10,label="h&r")
    pp.scatter(time_fragment_increasenum,[1.2]*len(time_fragment_increasenum),color='cyan',alpha=0.4,s=10,label="frag increase")
    pp.scatter(time_fragment_numstaysame,[1]*len(time_fragment_numstaysame),color='green',alpha=0.4,s=10,label="frag stay same")

    pp.legend(loc='upper left')

    pp.xlabel("Time (yr)")
    pp.ylabel("Type of collision")
    if x_limits != None:
        pp.xlim(x_limits)

    return fig
