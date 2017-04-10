# Created by JJW Mar. 13 2016
# Originally designed to plot average ecc. as a function of time
# for just the final bodies


import numpy as np
import matplotlib.pyplot as pp
import json
import os
import sys
import glob
path = os.path.abspath(os.path.join(os.path.dirname(__file__), '/u/joshuajw/PlanetProject/mercury'))
if not path in sys.path:
    sys.path.insert(1, path)
import mercury_outputreader as mercury
import mercury_output_looker_ater as looker_ater

def plot_average_ecc_func_time(time_values_to_use,ecc_values_to_use,num_body_func_time,second_time=None):
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

    line1 = ax.plot(time_values_to_use,avg_ecc,label='<e>',lw=2.2,color='red')
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

    return fig
