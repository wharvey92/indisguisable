import numpy
import math
import operator

import binascii
# Methods common to both the transmitter and receiver.
def local_carrier(fc, n, samplerate, name=None):
    '''
    Generate a cosine waveform at fc/samplerate for n samples.
    '''
    if name is None:
        args = numpy.arange(0, n) * fc * 2 * math.pi / samplerate
        return numpy.cos(args)
    elif name == "demodquad":        # used in quadrature DEMODULATION
        args = numpy.arange(0, n) * fc * 2 * math.pi / samplerate
        return numpy.exp(1j*args)



def bits2samples(bits, one, spb):
        '''
        Expand bits to samples by replacing each bit with self.spb samples.
        Each bit is mapped to corresponding sample using onebit2onesample.
        '''
        xsamples = [0.0]*len(bits)*spb
        for i in range(len(xsamples)):
            curr_bit = bits[int(i/spb)]
            if curr_bit == 1:
               xsamples[i] = one
            else:
               xsamples[i] = 0
        return xsamples


def get_barker():
    barker =  [1,0,1,1,0,1,1,1,0,0,0,1,1,1,1,1,0,0,1,1,0,1,0,1] # Barker sequence used in 802.11b
    #barker = [1,1,1,1,1,0,1,1,1,1,0,0,1,1,1,0,1,0,1,1,0,0,0,0,1,0,1,1,1,0,0,0,1,1,0,1,1,0,1,0,0,1,0,0,0,1,0,0,1,1,0,0,1,0,1,0,1,0,0,0,0,0,0]
    return barker


def set_preamble(silence):
    return numpy.append([0]*silence, get_barker())

def get_header_len():
    return 20

def get_rep():
    return 9

def lpfilter(samples_in, omega_cut):
  '''
  A low-pass filter of frequency omega_cut.
  '''
  # set the filter unit sample response
  L = 50
  hf = numpy.sin(omega_cut*numpy.arange(-L,0))/(math.pi*numpy.arange(-L,0))
  hf = numpy.append(hf, omega_cut)
  hf = numpy.append(hf, numpy.sin(omega_cut*numpy.arange(1,L+1))/(math.pi*numpy.arange(1,L+1)))
  hf[L] = omega_cut/math.pi
  # convolve unit sample response with input samples
  samples_out=[0]*len(samples_in)
  for i in range(len(samples_out) - len(hf)):
    samples_out[i] = numpy.sum(samples_in[i:i+len(hf)]*hf[::-1])

  return numpy.array(samples_out)



