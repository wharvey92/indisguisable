import numpy
import math
import operator

import binascii
# Methods common to both the transmitter and receiver.
def hamming(s1,s2):
    l = min(len(s1), len(s2))
    hd = sum(map(operator.xor,s1[:l],s2[:l]))
    err = 1.0*hd/float(l)
    return hd, err    

def get_headerlen():
    return 18

def get_hextlen():
    return 208
