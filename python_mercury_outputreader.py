#Created 02172016 by Joshua Wallace
#This defines various functions that will be useful for reading the output
# of the mercury7 code I am using.
#
#AGBTG

import numpy as np

class aei_info:
    def __init__(self,time_,a_,e_,i_,mass_):
        self.time = time_
        self.a    = a_
        self.e    = e_
        self.i    = i_
        self.mass = mass_

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

    while True:
        try:
            line = f.next()
            if line[0] != '-------------------------------------------------------------------':
                continue
            else:
                #Pull out the necessary information
        except StopIteration:
            break


