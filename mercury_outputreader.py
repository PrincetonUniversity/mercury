#Created 02172016 by Joshua Wallace
#This defines various functions that will be useful for reading the output
# of the mercury7 code I am using.
#
#AGBTG

import glob
import numpy as np
import matplotlib.pyplot as pp
#import matplotlib 
#print matplotlib.__version__
import os
import warnings
from matplotlib.ticker import MaxNLocator


long_string_of_dashes = "---------------------------------------------------------------------"

class aei_info:
    def __init__(self,time_,a_,e_,i_,mass_):
        try:
            assert( len(time_) == len(a_))
            assert( len(time_) == len(e_))
            assert( len(time_) == len(i_))
            assert( len(time_) == len(mass_))
        except TypeError:
            assert( isinstance(time_,float))
            assert( isinstance(a_,float))
            assert( isinstance(e_,float))
            assert( isinstance(i_,float))
            assert( isinstance(mass_,float))
            time_ = [time_]
            a_    = [a_]
            e_    = [e_]
            i_    = [i_]
            mass_ = [mass_]

        self.time = time_
        self.a    = a_
        self.e    = e_
        self.i    = i_
        self.mass = mass_


class aei_singletime:
    def __init__(self,a_,e_,i_,mass_):
        #assert( not isinstance(a_,list))
        #assert( not isinstance(e_,list))
        #assert( not isinstance(i_,list))
        #assert( not isinstance(mass_,list))
        self.a     = a_
        self.e     = e_
        self.i     = i_
        self.mass  = mass_


class collision_type:
    SIMPLE_MERGER = 1 #impact velocity less than mutual escape velocity
    EFFECTIVE_MERGER = 2 #if non-grazing regime and fragment mass less than minimum
    NONGRAZING_FRAG = 3   #if non-grazing and fragments
    GRAZE_MERGER = 4  #if grazing regime and impact velocity less than v2gm
    HIT_AND_RUN = 5   #if grazing regime and largest fragment larger than target mass
    GRAZING_FRAG = 6  #if grazing regime and fragmentstime

class collision_information:
    def __init__(self,target_name_,projectile_name_,time_,classification_,vimpact_vescape_ratio_,vgrazemerge_vescape_ratio_,B_Rtarg_ratio_,masslargestremnant_msum_ratio_,masslargestremnant_mtarget_ratio_,mfrag_mfragmin_ratio_):
        self.target_name = target_name_
        self.projectile_name = projectile_name_
        self.time = time_
        self.classification = classification_
        self.vimpact_vescape_ratio = vimpact_vescape_ratio_
        self.vgrazemerge_vescape_ratio = vgrazemerge_vescape_ratio_
        self.B_Rtarg_ratio = B_Rtarg_ratio_
        self.masslargestremnant_msum_ratio = masslargestremnant_msum_ratio_
        self.masslargestremnant_mtarget_ratio = masslargestremnant_mtarget_ratio_
        self.mfrag_mfragmin_ratio = mfrag_mfragmin_ratio_


class central_collision_information:
    def __init__(self,name_,time_):
        self.name = name_
        self.time = time_

class ejection_information:
    def __init__(self,name_,time_):
        self.name = name_
        self.time = time_

#class clo_info:
#    def __init__(self,time_,a_,e_,i_,mass_):
#        self.time = time_
#        self.a    = a_
#        self.e    = e_
#        self.i    = i_
#        self.mass = mass_

class numberofbodies_functime:
    def __init__(self,time_,number_):
        self.time = time_
        self.numberbodies = number_


def aei_aggregator(path_="./"):
    """This reads in all the *.aei files in the given path
    and returns a tuple, the first element of which is the names of all the bodies
    and the second element of which is a list of aei_info instances"""
    if not isinstance(path_, str):
        raise TypeError("The given path is not a string!")
    filelist, body_names = get_files("aei",path=path_)
    
    files_with_no_info = 0
    aei_list = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for filename in filelist:
            temp = aei_file_reader(filename)
            if temp == False:
                files_with_no_info += 1
            else:
                aei_list.append(temp)

    if len(filelist) != len(aei_list):
        if len(filelist) != ( len(aei_list) + files_with_no_info):
            raise TypeError("Lengths are not the same between the two things to return! Why?")
        else:
            print "There were " + str(files_with_no_info) + " aei files with no info"

    return (body_names,aei_list)


def aei_file_reader(filename):
    """This reads in a specified *.aei file and returns an
    instance of the aei_info class, corresponding to the object"""
    temp = np.loadtxt(filename,unpack=True)
    try:
        toreturn = aei_info(temp[0].tolist(),temp[1].tolist(),temp[2].tolist(),temp[3].tolist(),temp[4].tolist())
    except IndexError:
        toreturn = False #In this case, there is no information on the object in the file, and we return false
    return toreturn


def aei_func_time(aei_info_list):
    """This function takes a list of aei_info objects and returns,
    as a function of time, the aei's of all objects at the corresponding time.
    Returnts a tuple: list of time, list of aei_info objects, and list of number of objects"""
    if not isinstance(aei_info_list,list):
        raise TypeError("Umm, this isn't even a list.  I don't know what to do with this")
    if not isinstance(aei_info_list[0],aei_info):
        raise TypeError("This list doesn't contain aei_info instances, at least not the first element")

    i_with_maxlength = None
    maxlength = 0
    for i in range(len(aei_info_list)):
        temp = len(aei_info_list[i].time)
        if temp > maxlength:
            i_with_maxlength = i
            maxlength = temp

    print "Index of body with max length: " + str(i_with_maxlength)
            
    number_of_bodies_withmaxlength = 0
    for i in range(len(aei_info_list)):
        if len(aei_info_list[i].time) == maxlength:
             number_of_bodies_withmaxlength += 1
             if not np.array_equal(aei_info_list[i].time , aei_info_list[i_with_maxlength].time):
                 raise TypeError("There are two same-lengthed time arrays with differing values")

    print "Number of bodies with same time array as max length: " + str(number_of_bodies_withmaxlength)

    time_list = aei_info_list[i_with_maxlength].time
    aei_list_full = []

    for i in range(len(time_list)):
        time_lookingat = time_list[i]
        a_just_thistime = []
        e_just_thistime = []
        i_just_thistime = []
        mass_just_thistime = []
        for j in range(len(aei_info_list)):
            try:
                indextouse = aei_info_list[j].time.index(time_lookingat)
                a_just_thistime.append(aei_info_list[j].a[indextouse])
                e_just_thistime.append(aei_info_list[j].e[indextouse])
                i_just_thistime.append(aei_info_list[j].i[indextouse])
                mass_just_thistime.append(aei_info_list[j].mass[indextouse]) 
            except ValueError:
                continue

        aei_list_full.append(aei_singletime(a_just_thistime, e_just_thistime, i_just_thistime, mass_just_thistime) )

    number_of_objects = [ len(aei_list_full[k].a) for k in range(len(aei_list_full)) ]

    return (time_list,aei_list_full,number_of_objects)
        

#def clo_file_reader(filename):
#    """This reads in a specified *.aei file and returns an
#    instance of the aei_info class, corresponding to the object"""
#    temp = np.loadtxt(filename,unpack=True)
#    toreturn = aei_info(temp[0],temp[1],temp[2],temp[3],temp[4])
#    return toreturn

#Problem I see: Object column won't be a number in general, need to save as string


def get_files(extension,path="./"):
    """This returns a tuple, the first element of which is a 
    list of all the files in the current directory
    that end with the given string 'extension'.  No need to include wildcard
    or period. It also returns, as the second element of the tuple,
    the names of all the bodies (basically, the part of the filename
    that precedes the extension."""
    if not isinstance(path, str):
        raise TypeError("The given path is not a string!")
    temp = glob.glob(path + "*." + extension)
    if len(temp) == 0:
        raise TypeError("Found no files with path " + path + " and extension ." + extension)
    body_names = []
    for item in temp:
        body_names.append(os.path.splitext(os.path.split(item)[-1])[0])
    return(temp,body_names)
    


def numberofbodies_functime_reader(filename):
    """This reads in a specified *.stdout file and returns a list of
    instances of the numberofbodies_functime class"""
    f = open(filename, 'r')
    toreturn = []
    for line in f:
        if line[0:21] != " Number of big bodies":
            continue
        temp = line.split()
        if len(temp) !=10:
            raise TypeError("The number of bodies line doesn't have the right number of words/numbers")
        toreturn.append( numberofbodies_functime(float(temp[6]),int(temp[9])) )

    f.close()
    return toreturn


def collision_info_extractor(filename):
    """This reads in a specified *.info file and returns a tuple, 
    the first element of which is a list of collision_information instances,
    the second element of which is a list of central_collision_information 
    instances, and the third of which is a list of ejection_information 
    instances"""
    collision_info = []
    central_collision_info = []
    ejection_info = []
    f = open(filename, 'r')

    info_found = False
    

    while True:
        M2_info = False
        try:
            line = f.next()
            if line[0:66] == long_string_of_dashes[0:66]:
                info_found = True
                #Pull out the necessary information
                info = f.next().split()
                mass_ratio = float(info[-1])
                info = f.next().split()
                B_Rtarg_ratio = float(info[-1])
                info = f.next().split()
                vimpact_vescape_ratio = float(info[-1])
                info = f.next().split()
                vgrazemerge_vescape_ratio = float(info[-1])
                f.next()
                info = f.next().split()
                masslargestremnant_msum_ratio = float(info[-1])
                info = f.next().split()
                masslargestremnant_mtarget_ratio = float(info[-1])
                info = f.next().split()
                mfrag_mfragmin_ratio = float(info[-1])

                #Some will have this next line, others will not
                info = f.next().split()
                try:
                    if info[0] == "M2":
                        M2_info = True
                        M2_mproj_ratio = float(info[-1])
                        f.next()
                except IndexError:
                    pass
                info = f.next().split()
                projectile_name = info[0]
                target_name = info[-4]
                time = float(info[-2])
                classification = -1
                #Now to classify the collision
                if info[1] == "simply" and info[2] == "merged":
                    classification = collision_type.SIMPLE_MERGER
                elif info[1] == 'effectively' and info[2] == 'merged':
                    classification = collision_type.EFFECTIVE_MERGER
                elif info[1] == 'head-on' and info[2] == 'smashed':
                    classification = collision_type.NONGRAZING_FRAG
                elif info[1] == "grazed" and info[3] == "merged":
                    classification = collision_type.GRAZE_MERGER
                elif info[1] == "hit" and info[3] == "run":
                    classification = collision_type.HIT_AND_RUN
                elif info[1] == 'grazing' and info[2] == 'smashed':
                    classification = collision_type.GRAZING_FRAG
                else:
                    raise TypeError("Don't recognize collision type" + info[1] + " " + info[2] + " " + info[3])
                collision_info.append( collision_information(target_name,projectile_name,time,classification,vimpact_vescape_ratio,vgrazemerge_vescape_ratio,B_Rtarg_ratio,masslargestremnant_msum_ratio,masslargestremnant_mtarget_ratio,mfrag_mfragmin_ratio) )
            elif "collided with the central body" in line:
                info_found = True
                temp = line.split()
                central_collision_info.append( central_collision_information(temp[0],float(temp[-2]) ) )
            elif "ejected at" in line:
                info_found = True
                temp = line.split()
                ejection_info.append( ejection_information(temp[0],float(temp[-2]) ) )
        except StopIteration:
            break

    if info_found == True:
        return (collision_info, central_collision_info, ejection_info)
    else:
        raise TypeError("Did not find any information in file " + filename)

##########################
##########################
###Plotting functions

##########################
##########################
##########################
def mass_to_pointsize_converter(mass,scale=100):
    """Calculates a rough radius from the mass
    and scales it according to a given parameter, for the size of 
    circle to be plotted.  Default scale radius is sqrt(20) for an earth mass"""

    two_thirds = 2./3.
    temp = np.multiply(scale, np.power(np.multiply(mass,333000.),two_thirds)) #Numerical factor in parentheses converts to earth mass from solar
    #The above takes the given scale and then multiplies it by the radius (mass to the one-third power) squared (which gives the two-thirds power)
    return temp


def plot_aei_multiple(time_values_to_use,times_list,aeis_list,parameter_1,parameter_2,round_time=True,number_of_digits_to_round_to=1,year_unit="Myr",ylimits=None,xlimits=None):
    """This will plot a bunch of aei info similar to what John Chambers
    does.  time_values_to_use shows the time values that you want to plot
    at (the closest ones will be chosen)
    time_list is a list of the times outputted from aei_func_time,
    aeis_list is a list of the aei info outputted from aei_func_time,
    parameter_1 is the vertical axis parameter, parameter_2 is the 
    horizontal axis parameter, round_time asks whether you want to round
    the time values, and number_of_digits_to_round_to is the number 
    of digits to round to.
    """

    param_name_dict = {'e':"Eccentricity", 'a':"Semi-Major Axis", 'i':"Inclination", 'mass':"Mass"}
    param_unit_dict = {'e':"", 'a':" (AU)", 'i':" (degrees)", 'mass':" (Mass Units)"}
    param_limit_dict= {'e':(0,1),'a':(.2,2.6),'i':(0,90),'mass':(0,2e-5)}
    year_unit_dict  = {"Myr":1.e6}

    if not ( (parameter_1 in param_name_dict) and (parameter_2 in param_name_dict) ):
        raise TypeError("I can't recognize at least one of the two parameters given to me, " + parameter_1 + "  " + parameter_2)

    if year_unit not in year_unit_dict:
        raise TypeError("I do not recognize this (case-sensitive) year unit name: " + year_unit)

    if len(time_values_to_use) > 6:
        raise TypeError("I can't plot more than six times at a time!")

    figuresizex = 8.0
    figuresizey = 8.0
    lowerx = .082
    lowery = .093
    upperx = .97
    uppery = .97
    xwidth = (upperx-lowerx)/2.
    ywidth = (uppery-lowery)/3.

    fig = pp.figure(figsize=(figuresizex,figuresizey))
    if len(time_values_to_use) >= 5:
        ax1 = fig.add_axes([lowerx,lowery+2*ywidth,xwidth,ywidth])
        ax2 = fig.add_axes([lowerx+xwidth,lowery+2*ywidth,xwidth,ywidth])
        ax3 = fig.add_axes([lowerx,lowery+ywidth,xwidth,ywidth])
        ax4 = fig.add_axes([lowerx+xwidth,lowery+ywidth,xwidth,ywidth])
        ax5 = fig.add_axes([lowerx,lowery,xwidth,ywidth])
        ax6 = fig.add_axes([lowerx+xwidth,lowery,xwidth,ywidth])
        axlist = [ax1,ax2,ax3,ax4,ax5,ax6]
        ax_ylist = [ax3,ax5]
        for ax in (ax1,ax3,ax5):
            ax.set_ylabel(param_name_dict[parameter_1] + param_unit_dict[parameter_1],size=16)
        for ax in (ax2,ax4,ax6):
            ax.set_yticklabels([])
        for ax in (ax5,ax6):
            ax.set_xlabel(param_name_dict[parameter_2] + param_unit_dict[parameter_2],size=16)
        for ax in (ax1,ax2,ax3,ax4):
            ax.set_xticklabels([])
    else:
        ax1 = fig.add_axes([lowerx,lowery+ywidth,xwidth,ywidth])
        ax2 = fig.add_axes([lowerx+xwidth,lowery+ywidth,xwidth,ywidth])
        ax3 = fig.add_axes([lowerx,lowery,xwidth,ywidth])
        ax4 = fig.add_axes([lowerx+xwidth,lowery,xwidth,ywidth])
        axlist = [ax1,ax2,ax3,ax4]
        ax_ylist = [ax3]
        for ax in (ax1,ax3):
            ax.set_ylabel(param_name_dict[parameter_1] + param_unit_dict[parameter_1],size=16)
        for ax in (ax3,ax4):
            ax.set_xlabel(param_name_dict[parameter_2] + param_unit_dict[parameter_2],size=16)
        for ax in (ax2,ax4):
            ax.set_xlabel(param_name_dict[parameter_2] + param_unit_dict[parameter_2],size=16)
        for ax in (ax1,ax2):
            ax.set_xticklabels([])

    for ax in (axlist):  #This sets the x and y limits
        if xlimits==None:
            ax.set_xlim(param_limit_dict[parameter_2])
        else:
            ax.set_xlim(xlimits)
        if ylimits==None:
            ax.set_ylim(param_limit_dict[parameter_1])
        else:
            ax.set_ylim(ylimits)

    fig.canvas.draw() #So that later I can update the y-axis labels


    #Now, to hide the top ylabel value in some of the axes, the ones that are getting covered up

    for ax in ax_ylist:
        labels = [tick.get_text() for tick in ax.get_yticklabels()]
        ax.set_yticklabels(labels[:-1]) 

    for i in range(len(time_values_to_use)): #Plot the data
        argument = min(range(len(times_list)), key=lambda j: abs(times_list[j]-time_values_to_use[i]))
        print "Argument associated with time " + str(time_values_to_use[i]) + " is " + str(argument)
        print "    Time is " + str(times_list[argument]) + "  instead of " + str(time_values_to_use[i])
        if time_values_to_use[i] > times_list[-1] and (time_values_to_use[i] - times_list[-1]) > 1.01*(times_list[-1] - times_list[-2]): #If the wanted time exceeds the greatest time available by an amount greater than time printout step size
            raise TypeError("The desired time is greater than the times that data are available for")

        if 'mass' not in (parameter_1,parameter_2): #If mass is not one of the axes
            point_sizes = mass_to_pointsize_converter(aeis_list[argument].mass)
            axlist[i].scatter( getattr(aeis_list[argument],parameter_2), getattr(aeis_list[argument],parameter_1), s=point_sizes,facecolors='none')
        else:
            axlist[i].scatter( getattr(aeis_list[argument],parameter_2), getattr(aeis_list[argument],parameter_1))

    #Now, put time labels on each plot
    for i in range(len(time_values_to_use)):
        if abs(time_values_to_use[i] - 0.0) <= 1e-9:
            text = "Start"
        else: #round_time=False,number_of_digits_to_round_to=1,
            if round_time == True:
                time = time_values_to_use[i]
                order_of_magnitude = (np.log10(int(time_values_to_use[i])))
                taken_to_the_bottom = time / 10**(int(order_of_magnitude)+1)
                time = round(taken_to_the_bottom,number_of_digits_to_round_to)
                time = time * (10**(int(order_of_magnitude)+1)/year_unit_dict[year_unit])
                if time % 1 <= 1.e-3: #If the error introduced by multiplying this back up is small
                    time = round(time)
            else:
                time = time_values_to_use[i]/year_unit_dict[year_unit]

            if time >= 1.0:
                text = str(int(time)) + " " +year_unit
            else:#If time is a decimal
                text = str(time) + " " + year_unit
            if not time.is_integer():
                print "###### Possible warning!  Possible warning!  ######"
                print "The following cannot be represented as an integer with the given parameters:"

                print "     " + str(time_values_to_use[i]) + " will be represented as " + text

        axlist[i].text(0.03, 0.97,text, horizontalalignment='left',verticalalignment='top',transform=axlist[i].transAxes,size=16)

                
            


    return fig

def plot_number_func_time(filename="stdout.out"):
    """This will plot the number of big bodies as a function of time,
    as read from a specified file name.  The default file name is stdout.out,
    since it is to the stdout that the number of bodies is written 
    in the code."""

    fig = pp.figure()
    num_func_time = numberofbodies_functime_reader(filename)
    t = []
    num = []
    for i in range(len(num_func_time)):
        t.append(num_func_time[i].time)
        num.append(num_func_time[i].numberbodies)

    if t[-1] > 1e9:
        unit = "Gyr"
        t = np.divide(t,1e9)
    elif t[-1] > 1e6:
        unit = "Myr"
        t = np.divide(t,1e6)
    elif t[-1] > 1e3:
        unit = "Kyr"
        t = np.divide(t,1e3)
    else:
        unit = "years"

    pp.step(t,num,where='post',lw=2)
    pp.xlabel("Time (" + unit + ")")
    pp.ylabel("Number of bodies")

    return fig


def plot_all_aeis_here(times=(0.,3e6,10e6,30e6,60e6,300e6),a_limits=None):
    names, aei_functime = aei_aggregator()    
    times_aei_output, aeis, numbers = aei_func_time(aei_functime)

    fig = plot_aei_multiple(times,times_aei_output,aeis,'e','a',number_of_digits_to_round_to=2,xlimits=a_limits)
    fig.savefig("e_vs_a.pdf")

    fig = plot_aei_multiple(times,times_aei_output,aeis,'i','a',number_of_digits_to_round_to=2,ylimits=(0,45),xlimits=a_limits)
    fig.savefig("i_vs_a.pdf")

    fig = plot_aei_multiple(times,times_aei_output,aeis,'mass','a',number_of_digits_to_round_to=2,xlimits=a_limits)
    fig.savefig("mass_vs_a.pdf")

    fig = plot_aei_multiple(times,times_aei_output,aeis,'e','i',number_of_digits_to_round_to=2,xlimits=(0,45))
    fig.savefig("e_vs_i.pdf")



if __name__ == '__main__':
    """temp, tempcentral, tempejection = collision_info_extractor("info.out")
    print temp

    temp, tempcentral, tempejection = collision_info_extractor("info.2.out")
    print len(tempcentral)
    print len(tempejection)

    print "     "
    print "     "
    listoffiles = get_files("for")
    print listoffiles[0]
    print listoffiles[1]


    #################

    temp = get_files("aei")
    print temp"""

    names, aei_functime = aei_aggregator()
    
    times, aeis, numbers = aei_func_time(aei_functime)

    #print numbers

    print "     "
    print "     "

    fig = plot_aei_multiple((0.,3e6,10e6,30e6,60e6,299e6),times,aeis,'e','a',number_of_digits_to_round_to=2)
    fig.savefig("tempfig.pdf")