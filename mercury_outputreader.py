#Created 02172016 by Joshua Wallace
#This defines various functions that will be useful for reading the output
# of the mercury7 code I am using.
#
#AGBTG

import glob
import numpy as np
import math
import matplotlib.pyplot as pp
#import matplotlib 
#print matplotlib.__version__
import os
import warnings
from matplotlib.ticker import MaxNLocator


long_string_of_dashes = "---------------------------------------------------------------------"

class aei_info:
    def __init__(self,name_,time_,a_,e_,i_,mass_):
        try:
            assert( len(time_) == len(a_))
            assert( len(time_) == len(e_))
            assert( len(time_) == len(i_))
            assert( len(time_) == len(mass_))
        except TypeError:
            assert( isinstance(time_,float))
            assert( isinstance(name_,str))
            assert( isinstance(a_,float))
            assert( isinstance(e_,float))
            assert( isinstance(i_,float))
            assert( isinstance(mass_,float))
            time_ = [time_]
            a_    = [a_]
            e_    = [e_]
            i_    = [i_]
            mass_ = [mass_]

        self.name = name_
        self.time = time_
        self.a    = a_
        self.e    = e_
        self.i    = i_
        self.mass = mass_


class aei_singletime:
    def __init__(self,name_,a_,e_,i_,mass_):
        #assert( not isinstance(a_,list))
        #assert( not isinstance(e_,list))
        #assert( not isinstance(i_,list))
        #assert( not isinstance(mass_,list))
        self.name  = name_ #Names of the objects
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
    GRAZING_FRAG = 6  #if grazing regime and fragments

"""class collision_name_mass_coupled:
    def __init__(self,name_,mass_):
        self.name = name_
        self.mass = mass_"""

class collision_information:
    def __init__(self,target_name_,projectile_name_,time_,classification_,mt_over_mp_,vimpact_vescape_ratio_,vgrazemerge_vescape_ratio_,B_Rtarg_ratio_,masslargestremnant_msum_ratio_,masslargestremnant_mtarget_ratio_,mfrag_mfragmin_ratio_,number_of_fragments_,radius_,modified_escape_velocity_ratio_,names_,masses_):
        self.target_name = target_name_
        self.projectile_name = projectile_name_
        self.time = time_
        self.classification = classification_
        self.mt_over_mp = mt_over_mp_
        self.vimpact_vescape_ratio = vimpact_vescape_ratio_
        self.vgrazemerge_vescape_ratio = vgrazemerge_vescape_ratio_
        self.B_Rtarg_ratio = B_Rtarg_ratio_
        self.masslargestremnant_msum_ratio = masslargestremnant_msum_ratio_
        self.masslargestremnant_mtarget_ratio = masslargestremnant_mtarget_ratio_
        self.mfrag_mfragmin_ratio = mfrag_mfragmin_ratio_
        self.number_of_fragments = number_of_fragments_
        self.radius = radius_
        self.modified_escape_velocity_ratio = modified_escape_velocity_ratio_
        self.remnants = []
        if len(names_) != len(masses_):
            raise RuntimeError("names_ and masses_ not same length")
        self.remnant_names = names_
        self.remnant_masses = masses_
        


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

class numberofbodies_outsideroche_functime:
    def __init__(self,time_,number_inside_roche_,number_outside_roche_):
        self.time = time_
        self.number_inside_roche = number_inside_roche_
        self.number_outside_roche = number_outside_roche_


def aei_aggregator(path_="./",just_original_bodies=False):
    """This reads in all the *.aei files in the given path
    and returns a tuple, the first element of which is the names of all the bodies
    and the second element of which is a list of aei_info instances.
    If you only want to aggregate the original bodies, and none of the fragments,
    set just_original_bodies=True"""
    if not isinstance(path_, str):
        raise TypeError("The given path is not a string!")
    filelist, body_names = get_files("aei",path=path_)
    if just_original_bodies == True:
        #filelist = [item in filelist if 'F' not in item]
        filelist[:], body_names[:] = zip(*((x, y) for (x, y) in zip(filelist, body_names) if 'F' not in x))
        
    
    files_with_no_info = 0
    aei_list = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for j in range(len(filelist)):
            temp = aei_file_reader(filelist[j],body_names[j])
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


def aei_file_reader(filename,bodyname="blank"):
    """This reads in a specified *.aei file and returns an
    instance of the aei_info class, corresponding to the object"""
    temp = np.loadtxt(filename,unpack=True)
    try:
        toreturn = aei_info(bodyname,temp[0].tolist(),temp[1].tolist(),temp[2].tolist(),temp[3].tolist(),temp[4].tolist())
    except IndexError:
        toreturn = False #In this case, there is no information on the object in the file, and we return false
    return toreturn


def aei_func_time(aei_info_list):
    """This function takes a list of aei_info objects and returns,
    as a function of time, the aei's of all objects at the corresponding time.
    Returnts a tuple: list of time, list of aei_singletime objects, and list of number of objects"""
    if not isinstance(aei_info_list,list):
        raise TypeError("Umm, this isn't even a list.  I don't know what to do with this")
    if not isinstance(aei_info_list[0],aei_info):
        raise TypeError("This list doesn't contain aei_info instances, at least not the first element")


    i_with_maxlength = final_body_determiner(aei_info_list)


    time_list = aei_info_list[i_with_maxlength].time
    aei_list_full = []

    for i in range(len(time_list)):
        time_lookingat = time_list[i]
        name_just_thistime=[]
        a_just_thistime = []
        e_just_thistime = []
        i_just_thistime = []
        mass_just_thistime = []
        for j in range(len(aei_info_list)):
            try:
                indextouse = aei_info_list[j].time.index(time_lookingat)
                name_just_thistime.append(aei_info_list[j].name)
                a_just_thistime.append(aei_info_list[j].a[indextouse])
                e_just_thistime.append(aei_info_list[j].e[indextouse])
                i_just_thistime.append(aei_info_list[j].i[indextouse])
                mass_just_thistime.append(aei_info_list[j].mass[indextouse]) 
            except ValueError:
                continue

        aei_list_full.append(aei_singletime(name_just_thistime,a_just_thistime, e_just_thistime, i_just_thistime, mass_just_thistime) )

    number_of_objects = [ len(aei_list_full[k].a) for k in range(len(aei_list_full)) ]

    return (time_list,aei_list_full,number_of_objects)
        

#def clo_file_reader(filename):
#    """This reads in a specified *.aei file and returns an
#    instance of the aei_info class, corresponding to the object"""
#    temp = np.loadtxt(filename,unpack=True)
#    toreturn = aei_info(temp[0],temp[1],temp[2],temp[3],temp[4])
#    return toreturn

#Problem I see: Object column won't be a number in general, need to save as string


def final_body_determiner(aei_info_list,I_want_list_of_final_body_indices=False):
    """This will look at a list of aei_info instances, determine which
    objects lasted the whole integration time and also determine
    which objects are present at the end.
    Replaces some of the code that was in aei_func_time so it can now be used
    more generally.
    I_want_list_of_final_body_indices indicates whether one would like
    this function to also return a list of indices of bodies at the end
    Returns a tuple, first element is the index of the array with the max length
    and the second element is a list of the 
    """

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

    if I_want_list_of_final_body_indices == False:
        return i_with_maxlength

    elif I_want_list_of_final_body_indices == True:

        final_body_indices = []
        temp = aei_info_list[i_with_maxlength].time[-1]
        for i in range(len(aei_info_list)):
            if aei_info_list[i].time[-1] == temp: #If has as last time the time of a body that lasted the whole time
                final_body_indices.append(i)

        print "Number of final bodies: " + str(len(final_body_indices))


        return (i_with_maxlength,final_body_indices)



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
    toreturn_split_over_roche = []
    for line in f:
        if line[0:21] == " Number of big bodies":
            temp = line.split()
            if len(temp) != 10:
                print temp
                raise TypeError("The number of bodies line doesn't have the right number of words/numbers " + str(len(temp)) )
            toreturn.append( numberofbodies_functime(float(temp[6]),int(temp[9])) )
        elif line[0:26] == "  Num. bodies inside Roche":
            temp = line.split()
            if len(temp) != 14:
                raise TypeError("The number of bodies in/outside Roche line doesn't have the right number of words/numbers")
            toreturn_split_over_roche.append( numberofbodies_outsideroche_functime(float(temp[6]),float(temp[8]),float(temp[13]) ) )
    f.close()
    return (toreturn,toreturn_split_over_roche)


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
    no_remnant_mergers = False

    reported_no_radius = False
    reported_no_vescmodified = False

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
                info = f.next().split()     
                if info: #If info is not an empty list
                    radius = float(info[-1])
                    f.next()
                else:
                    radius = None #Set it to None because it wasn't there to be recorded
                if info: #If info is not an empty list
                    modified_escape_velocity_ratio = float(info[-1])
                    f.next()
                else:
                    modified_escape_velocity_ratio = None
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
                number_of_fragments = None
                #Now to classify the collision
                if info[1] == "simply" and info[2] == "merged":
                    classification = collision_type.SIMPLE_MERGER
                    number_of_fragments = 0
                elif info[1] == 'effectively' and info[2] == 'merged':
                    classification = collision_type.EFFECTIVE_MERGER
                    number_of_fragments = 0
                elif info[1] == 'head-on' and info[2] == 'smashed':
                    classification = collision_type.NONGRAZING_FRAG
                elif info[1] == "grazed" and info[3] == "merged":
                    classification = collision_type.GRAZE_MERGER
                    number_of_fragments = 0
                elif info[1] == "hit" and info[3] == "run":
                    classification = collision_type.HIT_AND_RUN
                elif info[1] == 'grazing' and info[2] == 'smashed':
                    classification = collision_type.GRAZING_FRAG
                else:
                    raise RuntimeError("Don't recognize collision type" + info[1] + " " + info[2] + " " + info[3])

                names = []
                masses = []
                no_second_remnant = None
            #    if number_of_fragments == None: #If we couldn't assume a fragment number from above
                info = f.next().split()
                if info:
                    if info[0] == "Remnant:":
                        names.append(info[1])
                        masses.append(float(info[3]))
                        info = f.next().split()
                        if info:
                            if info[0] == "Remnant:":
                                no_second_remnant = False
                                if float(info[3]) < 1e-20: #If the mass of the second remnant is zero
                                    second_mass_matters = False #Don't count this towards total number of final bodies
                                else:
                                    second_mass_matters = True
                                    names.append(info[1])
                                    masses.append(float(info[3]))
                            else:
                                no_second_remnant = True
                        else:
                            no_second_remnant = True
                    else:
                        second_mass_matters = False #Don't count this towards total number of final bodies
                        no_second_remnant = True
                        if no_remnant_mergers == False:
                            no_remnant_mergers = True
                            #raise RuntimeError("Expected to see Remnant information, but information not recorded!")
                            print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                            print "************************************************************************"
                            print "        Backwards version, no Remnants available for mergers"
                            print "************************************************************************"
                            print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                else:
                    second_mass_matters = False #Don't count this towards total number of final bodies
                    no_second_remnant = True
                    if no_remnant_mergers == False:
                        no_remnant_mergers = True
                            #raise RuntimeError("Expected to see Remnant information, but information not recorded!")
                        print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                        print "************************************************************************"
                        print "        Backwards version, no Remnants available for mergers"
                        print "************************************************************************"
                        print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                f_num = 0
                if no_second_remnant == False:
                    info=f.next().split()
                if info: #If list is not empty
                    while(info[0] == "Fragment:"):
                        names.append(info[1])
                        masses.append(float(info[3]))
                        f_num += 1
                        info=f.next().split()
                        if not info: #If list is empty
                            break
                if second_mass_matters == True:
                    number_of_fragments = f_num
                else:
                    number_of_fragments = f_num - 1

                if number_of_fragments < -1: #If somehow the number is very negative
                    print target_name
                    print projectile_name
                    print time
                    print classification
                    raise TypeError("The number of fragments is going onto record as being negative, " + str(number_of_fragments))

                if number_of_fragments == -1: #This is probably because two minimum frag mass bodies collided and merged just because they can't get smaller
                    number_of_fragments = 0



                collision_info.append( collision_information(target_name,projectile_name,time,classification,mass_ratio,vimpact_vescape_ratio,vgrazemerge_vescape_ratio,B_Rtarg_ratio,masslargestremnant_msum_ratio,masslargestremnant_mtarget_ratio,mfrag_mfragmin_ratio,number_of_fragments,radius,modified_escape_velocity_ratio,names,masses) )

                #### Raise some backwards compatibility notices
                if (radius == None and reported_no_radius == False):
                    print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                    print "************************************************************************"
                    print "        Backwards version, no radius recorded"
                    print "************************************************************************"
                    print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                    reported_no_radius = True
                if (radius == None and reported_no_vescmodified == False):
                    print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                    print "************************************************************************"
                    print "        Backwards version, no radius recorded"
                    print "************************************************************************"
                    print ""; print ""; print ""; print ""; print ""; print ""; print ""; print ""
                    reported_no_vescmodified = True
                
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


def mass_func_time_by_collisions(names_of_bodies,initial_body_mass,final_time,filename="info.out"):
    """This will get mass as a function of time for bodies specified in names_of_bodies
    based on the output of collision_info_extractor.  Initial_body_mass is the initial mass
    of the specific bodies, final_time is the time at the end of the integration,
    and filename is the file that contains the collision information
    (usually info.out)
    """
    collisions, central_collisions, ejections = collision_info_extractor(filename)

    t_lists = [ [0] for i in range(len(names_of_bodies)) ]
    m_lists = [ [initial_body_mass] for i in range(len(names_of_bodies)) ]

    for i in range(len(collisions)):
        for j in range(len(collisions[i].remnant_names)):
            if collisions[i].remnant_names[j] in names_of_bodies:
                index = names_of_bodies.index(collisions[i].remnant_names[j])
                m_lists[index].append(collisions[i].remnant_masses[j])
                t_lists[index].append(collisions[i].time)

    # Now to add on the final time/final mass so as to complete the eventual plot           
    for i in range(len(t_lists)):
        t_lists[i].append(final_time)
        m_lists[i].append(m_lists[i][-1])

    return (m_lists,t_lists)









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


def plot_aei_multiple(time_values_to_use,times_list,aeis_list,parameter_1,parameter_2,round_time=True,number_of_digits_to_round_to=3,max_number_of_decimals=3,year_unit="Myr",ylimits=None,xlimits=None):
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
    year_unit_dict  = {"Myr":1.e6,"kyr":1.e3}

    if not ( (parameter_1 in param_name_dict) and (parameter_2 in param_name_dict) ):
        raise TypeError("I can't recognize at least one of the two parameters given to me, " + parameter_1 + "  " + parameter_2)

    if year_unit not in year_unit_dict:
        raise TypeError("I do not recognize this (case-sensitive) year unit name: " + year_unit)

    if len(time_values_to_use) > 6:
        raise TypeError("I can't plot more than six times at a time!")

    figuresizex = 8.0
    figuresizey = 8.0
    lowerx = .091
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

                time = time/year_unit_dict[year_unit]


                time_temp = time * 10**max_number_of_decimals
                time_temp = round(time_temp)
                time_temp = time_temp  * 10**(-max_number_of_decimals) #Now, the number of decimal places is properly truncated, and the final value rounded.

                #print time_temp
                if abs(time_temp)<=1e-7:
                    time_temp_2 = 0.0
                else:
                    time_temp_2 = round(time_temp, -int(math.floor(np.log10(time_temp))) + (number_of_digits_to_round_to - 1))

                    if (number_of_digits_to_round_to - int(np.log10(time_temp_2)) + 1) == (max_number_of_decimals - 1) and max_number_of_decimals >= 1: #May need to adjust some things
                        super_temp = time * 10**(max_number_of_decimals - 2)
                        super_temp = int(  (super_temp - int(super_temp))*10)
                        super_temp_2 = time_temp * 10**(max_number_of_decimals - 2)
                        super_temp_2 = int(  (super_temp_2 - int(super_temp_2))*10)
                        if super_temp == 4 and super_temp_2 == 5:
                            print "Premature rounding affected things that will show up in time values, now taking care of"

                            time_temp_2 = time_temp_2 - 10**(-(max_number_of_decimals - 1))
                if time_temp_2 == 0.0:
                    if number_of_digits_to_round_to == 1:
                        time_temp_2 = int(time_temp_2)
                else:
                    if np.log10(time_temp_2) >= number_of_digits_to_round_to:
                        time_temp_2 = int(time_temp_2)

                time = time_temp_2

                
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

def plot_collision_scatterplot(filename="info.out",whichones=None,title="",collision_info=None):              
    """This function will take a info.out file and plot up a scatter plot of all the
    collisions in the v/vesc and r/R_target plane.  It returns this figure."""

    markers_touse = ("o","^","s","D")
    colors_touse = ( (146./255.,0,0),(0,109./255.,219./255.),(36./255.,255./255.,36./255.),(219./255.,209./255.,0))

    if collision_info == None:
        collisions, central_collisions, ejections = collision_info_extractor(filename)
    else:
        collisions = collision_info[0]
        central_collisions = collision_info[1]
        ejections = collision_info[2]
    print "Number of central collisions:   " + str(len(central_collisions))
    print "Number of ejections:            " + str(len(ejections))
    print " "
    print "Number of collisions:           " + str(len(collisions))

    merger = []
    effective_merger = []
    graze_merger = []
    grow   = []
    erode  = []
    hitandrun_nofrag=[]
    hitandrun_frag  =[]

    if whichones == None:
        for i in range(len(collisions)):
#            if collisions[i].classification in (collision_type.SIMPLE_MERGER, collision_type.EFFECTIVE_MERGER, collision_type.GRAZE_MERGER):
#                merger.append(collisions[i])
            if collisions[i].classification == collision_type.SIMPLE_MERGER:
                merger.append(collisions[i])
            elif collisions[i].classification == collision_type.EFFECTIVE_MERGER:
                effective_merger.append(collisions[i])
            elif collisions[i].classification ==  collision_type.GRAZE_MERGER:
                graze_merger.append(collisions[i])
            elif collisions[i].classification == collision_type.HIT_AND_RUN:
                if collisions[i].number_of_fragments == 0:
                    hitandrun_nofrag.append(collisions[i])
                else:
                    hitandrun_frag.append(collisions[i])
            elif collisions[i].masslargestremnant_mtarget_ratio >= 1.0:
                grow.append(collisions[i])
            elif collisions[i].masslargestremnant_mtarget_ratio < 1.0:
                erode.append(collisions[i])
            elif np.isnan(collisions[i].masslargestremnant_mtarget_ratio):
                continue
            else:
                print collisions[i].classification
                print collisions[i].target_name
                print collisions[i].projectile_name
                print collisions[i].time
                raise TypeError("I don't know what to do with this body, as far as its collision is concerned, body " + str(i))

    else:
        for i in range(len(collisions)):
            if collisions[i].target_name in whichones or collisions[i].projectile_name in whichones:
#                if collisions[i].classification in (collision_type.SIMPLE_MERGER, collision_type.EFFECTIVE_MERGER, collision_type.GRAZE_MERGER):
#                    merger.append(collisions[i])
                if collisions[i].classification == collision_type.SIMPLE_MERGER:
                    merger.append(collisions[i])
                elif collisions[i].classification == collision_type.EFFECTIVE_MERGER:
                    effective_merger.append(collisions[i])
                elif collisions[i].classification ==  collision_type.GRAZE_MERGER:
                    graze_merger.append(collisions[i])
                elif collisions[i].classification == collision_type.HIT_AND_RUN:
                    if collisions[i].number_of_fragments == 0:
                        hitandrun_nofrag.append(collisions[i])
                    else:
                        hitandrun_frag.append(collisions[i])
                elif collisions[i].masslargestremnant_mtarget_ratio >= 1.0:
                    grow.append(collisions[i])
                elif collisions[i].masslargestremnant_mtarget_ratio < 1.0:
                    erode.append(collisions[i])
                elif np.isnan(collisions[i].masslargestremnant_mtarget_ratio):
                    continue
                else:
                    raise TypeError("I don't know what to do with this body, as far as its collision is concerned, body " + str(i))

    fig = pp.figure()

    print "percentage that are perfect mergers: " + str( float(len(merger))/float(len(collisions)) * 100.0) + " %"

    pp.scatter([ item.B_Rtarg_ratio for item in merger], [ item.vimpact_vescape_ratio for item in merger],marker=markers_touse[0],color=colors_touse[0],label="merger")
    pp.scatter([ item.B_Rtarg_ratio for item in effective_merger], [ item.vimpact_vescape_ratio for item in effective_merger],marker=markers_touse[0],color="orange",label="effective merger")
    pp.scatter([ item.B_Rtarg_ratio for item in graze_merger], [ item.vimpact_vescape_ratio for item in graze_merger],marker=markers_touse[0],color="cyan",label="grazing merger")
    pp.scatter([ item.B_Rtarg_ratio for item in grow],  [item.vimpact_vescape_ratio for item in grow],marker=markers_touse[1],color=colors_touse[1],label="frag & grow")
    pp.scatter([ item.B_Rtarg_ratio for item in erode],  [item.vimpact_vescape_ratio for item in erode],marker=markers_touse[2],color=colors_touse[2],label="frag & erode")
    pp.scatter([ item.B_Rtarg_ratio for item in hitandrun_nofrag], [item.vimpact_vescape_ratio for item in hitandrun_nofrag],marker=markers_touse[3],color="magenta",label="hit & run no frag")
    pp.scatter([ item.B_Rtarg_ratio for item in hitandrun_frag], [item.vimpact_vescape_ratio for item in hitandrun_frag],marker=markers_touse[3],color=colors_touse[3],label="hit & run frag")
    pp.axhline(y=1,color='black',linestyle='--')
    pp.xlabel("b/R$_{\mathrm{target} }$",size=16)
    pp.ylabel("$v/v_{esc}$",size=16)
    pp.ylim(bottom=0)
    pp.xlim(left=0)
    pp.legend(loc="best")
    pp.title(title)

    return fig


def plot_collision_scatterplot_simplifiedcollisionclassification(filename="info.out",whichones=None,title="",collision_info=None,outside_Roche=False):              
    """This function will take a info.out file and plot up a scatter plot of all the
    collisions in the v/vesc and r/R_target plane.  It returns this figure.
    This is the same as plot_collision_scatterplot, except it takes a simplified collision classification scheme 
    like in the figures for the paper.  It also can differentiate between collisions that took place inside and outside 
    the Roche radius."""

    markers_touse = ("o","s","D","^")
    markersize = 4
    colors_touse = ('#f4320c','#0d75f8','gray','#fcc006')#,'#fac205','goldenrod')

    if collision_info == None:
        collisions, central_collisions, ejections = collision_info_extractor(filename)
    else:
        collisions = collision_info[0]
        central_collisions = collision_info[1]
        ejections = collision_info[2]
    print "Number of central collisions:   " + str(len(central_collisions))
    print "Number of ejections:            " + str(len(ejections))
    print " "
    print "Number of collisions:           " + str(len(collisions))

    if outside_Roche == False:
        radius_to_care_about = 0.
    else:
        radius_to_care_about = 0.0088781786


    merger = []
    effective_merger = []
    grow   = []
    erode  = []
    hitandrun = []

    if whichones == None:
        for i in range(len(collisions)):
#            if collisions[i].classification in (collision_type.SIMPLE_MERGER, collision_type.EFFECTIVE_MERGER, collision_type.GRAZE_MERGER):
#                merger.append(collisions[i])
            if collisions[i].classification == collision_type.SIMPLE_MERGER:
                merger.append(collisions[i])
            elif collisions[i].classification == collision_type.EFFECTIVE_MERGER:
                effective_merger.append(collisions[i])
            elif collisions[i].classification ==  collision_type.GRAZE_MERGER:
                merger.append(collisions[i])
            elif collisions[i].classification == collision_type.HIT_AND_RUN:
                hitandrun.append(collisions[i])
            elif collisions[i].masslargestremnant_mtarget_ratio >= 1.0:
                grow.append(collisions[i])
            elif collisions[i].masslargestremnant_mtarget_ratio < 1.0:
                erode.append(collisions[i])
            elif np.isnan(collisions[i].masslargestremnant_mtarget_ratio):
                print "Warning!  Warning!  Had a NaN mass ratio"
                continue
            else:
                print collisions[i].classification
                print collisions[i].target_name
                print collisions[i].projectile_name
                print collisions[i].time
                raise TypeError("I don't know what to do with this body, as far as its collision is concerned, body " + str(i))

    else:
        for i in range(len(collisions)):
            if (collisions[i].target_name in whichones or collisions[i].projectile_name in whichones) and (collisions[i].radius > radius_to_care_about): ###Just the final bodies
                if collisions[i].classification == collision_type.SIMPLE_MERGER:
                    merger.append(collisions[i])
                elif collisions[i].classification == collision_type.EFFECTIVE_MERGER:
                    effective_merger.append(collisions[i])
                elif collisions[i].classification ==  collision_type.GRAZE_MERGER:
                    merger.append(collisions[i])
                elif collisions[i].classification == collision_type.HIT_AND_RUN:
                    hitandrun.append(collisions[i])
            #if collisions[i].number_of_fragments == 0:
            #    hitandrun_nofrag.append(collisions[i])
            #else:
            #    hitandrun_frag.append(collisions[i])
                elif collisions[i].masslargestremnant_mtarget_ratio >= 1.0:
                    grow.append(collisions[i])
                elif collisions[i].masslargestremnant_mtarget_ratio < 1.0:
                    erode.append(collisions[i])
                elif np.isnan(collisions[i].masslargestremnant_mtarget_ratio):
                    print "Warning!  Warning!  Had a NaN mass ratio"
                    continue
                else:
                    print collisions[i].classification
                    print collisions[i].target_name
                    print collisions[i].projectile_name
                    print collisions[i].time
                    raise RuntimeError("I don't know what to do with this body, as far as its collision is concerned, body " + str(i))





    fig = pp.figure()

    print "percentage that are perfect mergers: " + str( float(len(merger))/float(len(collisions)) * 100.0) + " %"

    pp.scatter([ item.B_Rtarg_ratio for item in erode],  [item.vimpact_vescape_ratio for item in erode],marker=markers_touse[2],color=colors_touse[2],s=markersize,label="erode")
    pp.scatter([ item.B_Rtarg_ratio for item in hitandrun], [item.vimpact_vescape_ratio for item in hitandrun],marker=markers_touse[3],color=colors_touse[3],s=markersize,label="hit & run")
    pp.scatter([ item.B_Rtarg_ratio for item in merger], [ item.vimpact_vescape_ratio for item in merger],marker=markers_touse[0],color=colors_touse[0],s=(markersize+8),label="merger")
    pp.scatter([ item.B_Rtarg_ratio for item in merger], [ item.vimpact_vescape_ratio for item in merger],marker=markers_touse[0],color=colors_touse[0],s=(markersize+8),facecolor='white',label="eff. merger")
    pp.scatter([ item.B_Rtarg_ratio for item in grow],  [item.vimpact_vescape_ratio for item in grow],marker=markers_touse[1],color=colors_touse[1],s=markersize,label="grow")

    pp.axhline(y=1,color='black',linestyle='--')
    pp.xlabel("b/R$_{\mathrm{target} }$",size=16)
    pp.ylabel("$v/v_{esc}$",size=16)
    pp.ylim(bottom=0)
    pp.xlim(left=0)
    pp.legend(loc="best")
    pp.title(title)

    return fig





def plot_number_func_time(filename="stdout.out",xscale='log'):
    """This will plot the number of big bodies as a function of time,
    as read from a specified file name.  The default file name is stdout.out,
    since it is to the stdout that the number of bodies is written 
    in the code."""
    
    fig = pp.figure()
    (num_func_time,split_over_roche) = numberofbodies_functime_reader(filename)
    t = []
    num = []
    for i in range(len(num_func_time)):
        t.append(num_func_time[i].time)
        num.append(num_func_time[i].numberbodies)

    t_roche = []
    num_inside = []
    num_outside = []
    for k in range(len(split_over_roche)):
        t_roche.append(split_over_roche[k].time)
        num_inside.append(split_over_roche[k].number_inside_roche)
        num_outside.append(split_over_roche[k].number_outside_roche)

    #print t
    if t[-1] > 1e9:
        unit = "Gyr"
        divide_by_number = 1e9
    elif t[-1] > 1e6:
        unit = "Myr"
        divide_by_number = 1e6
    elif t[-1] > 1e3:
        unit = "Kyr"
        divide_by_number = 1e3
    else:
        unit = "years"
        divide_by_number = 1.

    pp.step(np.divide(t,divide_by_number),num,where='post',lw=2.2,label='all bodies')

    del t
    del num

    if split_over_roche: #If there is information for number of bodies inside and outside roche
        t_roche = []
        num_inside = []
        num_outside = []
        for k in range(len(split_over_roche)):
            t_roche.append(split_over_roche[k].time/1000.0)
            num_inside.append(split_over_roche[k].number_inside_roche)
            num_outside.append(split_over_roche[k].number_outside_roche)

        pp.step(np.divide(t_roche,divide_by_number),num_outside,where='post',lw=1.5,label='outside Roche')
        pp.step(np.divide(t_roche,divide_by_number),num_inside,where='post',lw=1.5,label='inside Roche')

    pp.legend(loc='best')

    pp.xlabel("Time (" + unit + ")")
    if xscale == 'log':
        pp.xscale('log')
    elif xscale == 'linear':
        pp.xscale('linear')
    else:
        print "Unrecognized xscale argument!  Defaulting to log."
        pp.xscale('log')
    pp.ylabel("Number of bodies")

    return fig


def plot_all_aeis_here(times=(0.,3e6,10e6,30e6,60e6,300e6),a_limits=None,e_limits=None,names_and_aeifunctime=None,just_original_bodies=False,year_unit='kyr'):
    if names_and_aeifunctime == None:
        names, aei_functime = aei_aggregator()
    else:
        if len(names_and_aeifunctime) != 2:
            raise TypeError("names_and_aeifunctime was not length 2!")
        names = names_and_aeifunctime[0]
        aei_functime = names_and_aeifunctime[1]

    if e_limits==None:
        e_limits=(0,0.05)
    times_aei_output, aeis, numbers = aei_func_time(aei_functime)

    fig = plot_aei_multiple(times,times_aei_output,aeis,'e','a',number_of_digits_to_round_to=3,max_number_of_decimals=3,xlimits=a_limits,ylimits=e_limits,year_unit=year_unit)
    fig.savefig("e_vs_a.pdf")

    #fig = plot_aei_multiple(times,times_aei_output,aeis,'i','a',number_of_digits_to_round_to=2,ylimits=(0,45),year_unit=year_unit,xlimits=a_limits)
    #fig.savefig("i_vs_a.pdf")

    #fig = plot_aei_multiple(times,times_aei_output,aeis,'mass','a',number_of_digits_to_round_to=2,year_unit=year_unit,xlimits=a_limits)
    #fig.savefig("mass_vs_a.pdf")

    #fig = plot_aei_multiple(times,times_aei_output,aeis,'e','i',number_of_digits_to_round_to=2,year_unit=year_unit,xlimits=(0,45),ylimits=e_limits)
    #fig.savefig("e_vs_i.pdf")



"""if __name__ == '__main__':
    ""temp, tempcentral, tempejection = collision_info_extractor("info.out")
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
    print temp""

    names, aei_functime = aei_aggregator()
    
    times, aeis, numbers = aei_func_time(aei_functime)

    #print numbers

    print "     "
    print "     "

    fig = plot_aei_multiple((0.,3e6,10e6,30e6,60e6,299e6),times,aeis,'e','a',number_of_digits_to_round_to=2,year_unit=year_unit)
    fig.savefig("tempfig.pdf")
"""




def plot_mass_func_time(names_of_bodies,initial_body_mass=5.015010e-08,final_time=3e5,file_path="./info.out"):
    """
    This plots the mass of selected bodies as a function of time.  It returns a figure object, which can then
    be made into an image a saved.
    """
    fig = pp.figure()
    earthmass_converter = 3.0025e-6

    masses, times = mass_func_time_by_collisions(names_of_bodies,initial_body_mass,final_time,filename=file_path)

    divide_by = 1.0
    for i in range(len(masses)):
        pp.step(np.divide(times[i],divide_by),np.divide(masses[i],earthmass_converter),where='post',color='blue')

    pp.xlabel("Time (years)")
    pp.ylabel('Mass (M$_\oplus$)',size=18)
    pp.xscale(u'log')
    pp.xlim(left=1e-1,right=final_time)
    top_val = 0.6
    pp.ylim(top=top_val)

    return fig




##########################################################
##########################################################
##########################################################
##########################################################
####   Added later

def plot_collisions_mtovermp_func_time(collisions,whichonestofocus=None):
    """ This function will take a list of collision_information instances
    (optional list of names of bodies to highlight in the plot) 
    then plots mt over mp on the y-axis, logarithm time
    on the x-axis, and highlights the collisions of the objects that
    made it to the end,
    and returns a figure that can then be saved.
    """

    if not isinstance(collisions[0],collision_information):
        raise TypeError("You did not give this function a list of collision_information instances!")

    fig = pp.figure()

    point_size = 12
    times = []
    mt_over_mp = []
    for i in range(len(collisions)):
        times.append(collisions[i].time)
        mt_over_mp.append(collisions[i].mt_over_mp)

    pp.scatter(times,mt_over_mp,s=point_size,facecolor='none') #all of them

    #Now, just the final bodies
    if whichonestofocus != None:
        times = []
        mt_over_mp = []
        for i in range(len(collisions)):
            if (collisions[i].target_name in whichonestofocus) or (collisions[i].projectile_name in whichonestofocus):
                times.append(collisions[i].time)
                mt_over_mp.append(collisions[i].mt_over_mp)

        pp.scatter(times,mt_over_mp,color='blue',s=point_size)
    
    pp.xscale(u'log')
    pp.xlabel("Time (years)")
    pp.ylabel("Ratio of collider masses")
    pp.ylim(-.01,2.05)

    return fig


def calc_average_mtovermp(collisions,whichonestofocus):
    """This function will calculate the average mt_over_mp for 
    all collisions in collisions, as well as calculate the average
    mt_over_mp for the bodies listed in whichonestofocus
    """

    mt_over_mp_list = []
    target_name_list = []
    projectile_name_list = []
    for i in range(len(collisions)):
        mt_over_mp_list.append(collisions[i].mt_over_mp)
        target_name_list.append(collisions[i].target_name)
        projectile_name_list.append(collisions[i].projectile_name)


    if not mt_over_mp_list:
        mean = 0
        median = 0
        stddev = 0
    else:
        mean = np.mean(mt_over_mp_list)
        median = np.median(mt_over_mp_list)
        stddev = np.std(mt_over_mp_list)


    focused_list = []

    for i in range(len(whichonestofocus)):
        temp = whichonestofocus[i]
        mt_over_mp_thisone = []
        for j in range(len(collisions)):
            if (temp == target_name_list[j]) or (temp == projectile_name_list[j]):
                mt_over_mp_thisone.append(collisions[j].mt_over_mp)


        if not mt_over_mp_thisone:
            focused_list.append ( [0,0,0] )
        else:
            focused_list.append( [np.mean(mt_over_mp_thisone),np.median(mt_over_mp_thisone),np.std(mt_over_mp_thisone)] )





    return ([mean,median,stddev],focused_list)
