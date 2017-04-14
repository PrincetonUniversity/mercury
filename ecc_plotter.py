# Created by JJW Mar. 13 2016
# Originally designed to plot average ecc. as a function of time
# for just the final bodies


import numpy as np
import matplotlib.pyplot as pp
import matplotlib.markers as markers
import json
import os
import sys
import glob
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/joshuajw/PlanetProject/mercury'))
if not path in sys.path:
    sys.path.insert(1, path)
import mercury_outputreader as mercury
import mercury_output_looker_ater as looker_ater


def close_encounter_getter(list_of_bodies,threshold=0):
    """ This will extract all the close encounters from the bodies in list_of_bodies, 
    then return those that are closer than the given threshold (in AU)
    """
    
    close_encounters = []

    for body in list_of_bodies:
        with open(body + ".clo",'r') as f:
            lines = f.readlines()
            for line in lines:
                if line[0] != "#":
                    splitted = line.split()
                    close_encounters.append( [float(splitted[0]), float(splitted[2])] )

    encounters, times = zip(*sorted(zip([item[1] for item in close_encounters], [item[0] for item in close_encounters])))
    return times[:18]


def plot_average_ecc_func_time(time_values_to_use,ecc_values_to_use,num_body_func_time,second_time=None,ce_to_plot=None):
    """This will plot the average ecc. as a function of time, with the 
    number of bodies also as a function to time, for comparison

    time_values_to_use -- a list of the time values to be plotted
    ecc_values_to_use  -- a list of lists, each list entry is a list of the eccentricities at that time, from which to calculate the average ecc. at a given time
    num_body_func_time: a list of the number of bodies as a function of time
    """

    avg_ecc = []
    max_ecc = []
    for i in range(len(ecc_values_to_use)):
        max_ecc.append(max(ecc_values_to_use[i]) )
        total = np.sum(ecc_values_to_use[i])
        avg_ecc.append(total/float(len(ecc_values_to_use[i])) )


    fig = pp.figure()

    ax = fig.add_subplot(111)
    ax2 = ax.twinx()

    line1 = ax.plot(time_values_to_use,avg_ecc,label='<e>',lw=0.5,color='red')
    line2 = ax.plot(time_values_to_use,max_ecc,label='max e',lw=0.5,color='cyan')

    if second_time == None:
        second_time = time_values_to_use
    line3 = ax2.step(second_time,num_body_func_time,where='post',lw=2.2,label='numbodies',color='blue')
    ax.set_xscale('log')
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('ecc. of final planets')
    ax2.set_ylabel('Num bodies outside Roche')

    lines_for_legend = line1 + line2 + line3
    fig.legend(lines_for_legend, [l.get_label() for l in lines_for_legend], loc='best')

    if ce_to_plot != None:
        ax.scatter( ce_to_plot, [0.002]*len(ce_to_plot), marker=markers.CARETUP,color='red')

    ax.set_ylim(bottom=0)
    return fig


def plot_ecc_vector_func_time(time, aei_functime):
    """ This plots the mass-weighted eccentricity vector as a function of time"""

    ecc_vec = []

    for i in range(len(aei_functime)):
        sum_massweighted_eccvector_x = 0.
        sum_massweighted_eccvector_y = 0.
        sum_mass = 0.
        for j in range(len(aei_functime[i].mass)):
            sum_mass += aei_functime[i].mass[j] 
            sum_massweighted_eccvector_x += aei_functime[i].mass[j] * aei_functime[i].e[j] * np.cos(np.deg2rad(aei_functime[i].pomega[j]))
            sum_massweighted_eccvector_y += aei_functime[i].mass[j] * aei_functime[i].e[j] * np.sin(np.deg2rad(aei_functime[i].pomega[j]))

        vector_magnitude = np.sqrt(sum_massweighted_eccvector_x * sum_massweighted_eccvector_x + 
                                   sum_massweighted_eccvector_y * sum_massweighted_eccvector_y)
        ecc_vec.append(vector_magnitude/sum_mass)
    
    fig = pp.figure()
    ax  = fig.add_subplot(111)

    ax.plot(time, ecc_vec)
    ax.set_xlabel("Time (years)")
    ax.set_ylabel("Ecc vector magnitude")
    ax.set_xscale('log')

    return fig
    
