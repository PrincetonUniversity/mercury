#Created 02172016 by Joshua Wallace
#This defines various functions that will be useful for reading the output
# of the mercury7 code I am using.
#
#AGBTG

import glob
import numpy as np
import os


long_string_of_dashes = "---------------------------------------------------------------------"

class aei_info:
    def __init__(self,time_,a_,e_,i_,mass_):
        assert( len(time_) == len(a_))
        assert( len(time_) == len(e_))
        assert( len(time_) == len(i_))
        assert( len(time_) == len(mass_))
        self.time = time_
        self.a    = a_
        self.e    = e_
        self.i    = i_
        self.mass = mass_


class aei_singletime:
    def __init__(self,a_,e_,i_,mass_):
        assert( not isinstance(a_,list))
        assert( not isinstance(e_,list))
        assert( not isinstance(i_,list))
        assert( not isinstance(mass_,list))
        self.a     = a_
        self.e     = e_
        self.i     = i_
        self.mass = mass_


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
    def __init(self,time_,number_):
        self.time = time_
        self.numberbodies = number_


def aei_aggregator(path_="./"):
    """This reads in all the *.aei files in the given path
    and returns a tuple, the first element of which is the names of all the bodies
    and the second element of which is a list of aei_info instances"""
    if not isinstance(path_, str):
        raise TypeError("The given path is not a string!")
    filelist, body_names = get_files("aei",path=path_)

    aei_list = []
    for filename in filelist:
        aei_list.append(aei_file_reader(filename))

    if len(filelist) != len(aei_list):
        raise TypeError("Lengths are not the same between the two things to return! Why?")

    return (body_names,aei_list)


def aei_file_reader(filename):
    """This reads in a specified *.aei file and returns an
    instance of the aei_info class, corresponding to the object"""
    temp = np.loadtxt(filename,unpack=True)
    toreturn = aei_info(temp[0].tolist(),temp[1].tolist(),temp[2].tolist(),temp[3].tolist(),temp[4].tolist())
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
        aei_list_just_thistime = []
        for j in range(len(aei_info_list)):
            try:
                indextouse = aei_info_list[j].time.index(time_lookingat)
                aei_list_just_thistime.append(aei_singletime(aei_info_list[j].a[indextouse], aei_info_list[j].e[indextouse], aei_info_list[j].i[indextouse], aei_info_list[j].mass[indextouse]) )
            except ValueError:
                continue

        aei_list_full.append(aei_list_just_thistime)

    number_of_objects = [ len(aei_list_full[k]) for k in range(len(aei_list_full)) ]

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
        if line[0] != "Number of big bodies":
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


if __name__ == '__main__':
    temp, tempcentral, tempejection = collision_info_extractor("info.out")
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
    print temp

    names, aei_functime = aei_aggregator()
    
    times, aeis, numbers = aei_func_time(aei_functime)

    #print numbers
