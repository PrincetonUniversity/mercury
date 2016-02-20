#Created 02172016 by Joshua Wallace
#This defines various functions that will be useful for reading the output
# of the mercury7 code I am using.
#
#AGBTG

import numpy as np

long_string_of_dashes = "---------------------------------------------------------------------"

class aei_info:
    def __init__(self,time_,a_,e_,i_,mass_):
        self.time = time_
        self.a    = a_
        self.e    = e_
        self.i    = i_
        self.mass = mass_


class collision_type:
    SIMPLE_MERGER = 1 #impact velocity less than mutual escape velocity
    EFFECTIVE_MERGER = 2 #if non-grazing regime and fragment mass less than minimum
    NONGRAZING_FRAG = 3   #if non-grazing and fragments
    GRAZE_MERGER = 4  #if grazing regime and impact velocity less than v2gm
    HIT_AND_RUN = 5   #if grazing regime and largest fragment larger than target mass
    GRAZING_FRAG = 6  #if grazing regime and fragmentstime

class collision_object:
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



def aei_file_reader(filename):
    """This reads in a specified *.aei file and returns an
    instance of the aei_info class, corresponding to the object"""
    temp = np.loadtxt(filename,unpack=True)
    toreturn = aei_info(temp[0],temp[1],temp[2],temp[3],temp[4])
    return toreturn

#def clo_file_reader(filename):
#    """This reads in a specified *.aei file and returns an
#    instance of the aei_info class, corresponding to the object"""
#    temp = np.loadtxt(filename,unpack=True)
#    toreturn = aei_info(temp[0],temp[1],temp[2],temp[3],temp[4])
#    return toreturn

#Problem I see: Object column won't be a number in general, need to save as string


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
    """This reads in a specified *.info file and returns a list of
    collision_information instances"""
    toreturn = []
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
                toreturn.append( collision_object(target_name,projectile_name,time,classification,vimpact_vescape_ratio,vgrazemerge_vescape_ratio,B_Rtarg_ratio,masslargestremnant_msum_ratio,masslargestremnant_mtarget_ratio,mfrag_mfragmin_ratio) )

        except StopIteration:
            break

    if info_found == True:
        return toreturn
    else:
        raise TypeError("Did not find any information in file " + filename)


if __name__ == '__main__':
    temp = collision_info_extractor("info.out")
    print temp
