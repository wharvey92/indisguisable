import sys
import math
import numpy
import matplotlib.pyplot as p
import scipy.cluster.vq
import common_txrx as common
from graphs import *
from numpy import linalg as LA

import hamming_db

class Receiver:
    def __init__(self, carrier_freq, samplerate, spb):
        '''
        The physical-layer receive function, which processes the
        received samples by detecting the preamble and then
        demodulating the samples from the start of the preamble's
        Barker sequence. Returns the sequence of received bits (after
        demapping), soft information (the squared Euclidean distance
        between the symbol corresponding to each bit and the nominal
        (mean) received value for that bit, and the signal-to-noise
        ratio estimated from the preamble's Barker sequence.
        '''
        self.fc = carrier_freq
        self.samplerate = samplerate
        self.spb = spb 
        print 'Receiver: '
    def detect_threshold(self, demod_samples): 
        # Now, we have a bunch of values that, for on-off keying, are
        # either near amplitude 0 or near a positive amplitude
        # (corresp. to bit "1").  Because we don't know the balance of
        # zeroes and ones in the input, we use 2-means clustering to
        # determine the "1" and "0" clusters.  In practice, some
        # systems use a random scrambler to XOR the input to balance
        # the zeroes and ones. We have decided to avoid that degree of
        # complexity in audiocom (for the time being, anyway).
        one = max(scipy.cluster.vq.kmeans(demod_samples, 2)[0])
        zero = min(scipy.cluster.vq.kmeans(demod_samples, 2)[0])
        thresh = (one + zero)/ 2.0
        return one, zero, thresh
 
    def detect_preamble(self, demod_samples, thresh, one):

        # Find the sample corresp. to the first reliable bit "1"; this step 
        # is crucial to a proper and correct synchronization w/ the xmitter.
        offset = self.detect_one(demod_samples, thresh, one)
        if offset < 0:
            print '*** ERROR: Could not detect any ones (so no preamble). ***'
            print '\tIncrease volume / turn on mic?'
            print '\tOr is there some other synchronization bug? ***'
            sys.exit(1)

        barker = common.get_barker()
        presamples = common.bits2samples(barker, 1.0, self.spb)
        mod_presamples = presamples * common.local_carrier(self.fc, len(presamples), self.samplerate)
        m = self.demodulate(mod_presamples)
        pre_offset = self.correlate(m, demod_samples[offset:offset+3*len(barker)*self.spb])
        return offset + pre_offset
        
    def demap_and_check(self, demod_samples, barker_start):
        spb = self.spb
        subsamples = self.subsample(demod_samples[barker_start:], int(spb/4), int(3*spb/4))

        (bits,si,snr) = self.demap(subsamples)
        start = self.barker_check(bits)

        if start >= 0:
            print '\tI think the Barker sequence starts at bit', start + barker_start/spb, \
                '(sample', start*spb+barker_start,')'
            recd_bits = numpy.array(bits[start:], dtype=int)
            if snr > 0.0:
                print '\tSNR from preamble: %.1f dB' % (10.0*math.log(snr, 10))
            else:
                print '\tWARNING: Couldn\'t estimate SNR...'
            blen = len(common.get_barker())
            print '\tRecd preamble: ', recd_bits[:blen]
            return recd_bits[blen:]
        else:
            print '*** ERROR: Could not detect preamble. ***'
            print '\tIncrease volume / turn on mic / reduce noise?'
            print '\tOr is there some other synchronization bug? ***'
            sys.exit(1)            
         

    def demodulate(self, samples):
        return self.avgfilter(numpy.abs(samples), (self.samplerate/self.fc)/2)

####################
    def subsample(self, dsamples, start, end):
        subsamp = numpy.array([])
        for i in range(0, len(dsamples), self.spb):
            if i+start < len(dsamples):
                subsamp = numpy.append(subsamp, numpy.mean(dsamples[i+start:i+end]))
        return subsamp

    def detect_one(self, demod_samples, thresh, one):
        spb = self.spb
        for offset in xrange(len(demod_samples)):
            if numpy.mean(demod_samples[offset+spb/4:offset+3*spb/4]) > thresh + (one-thresh) / 2:
                return offset
        return -1

    def avgfilter(self, samples_in, window):
        '''
        A simple averaging filter whose every output is the
        mean of the last "window" samples.
        y[n] = (x[n] + x[n+1] + ... + x[n+window-1])/window, 
        where x is the input sequence of samples and y is the output
        '''
        y = [0]*len(samples_in)
        for i in range(len(y)):
            y[i] = numpy.sum(samples_in[i:i+window])/window
        return numpy.array(y)

    def correlate(self, x, y):
        '''
        Find the starting point in y of the best occurrence of x.
        Assume len(x) <= len(y). Returns -1 if len(x) > len(y).
        We assume that x shows up in y from the xmitter only once in the 
        stream, but if it happens to repeat, we'll return the "best" match.
        The match is done by correlating y with x, i.e., taking the dot product
        of the corresponding samples.
        '''
        if len(x) > len(y) or len(x) == 0:
            # actually an error, but we'll let the caller deal w/ it
            return 0 
        correl = []
        for i in range(len(y)-len(x)+1):
            correl.append(numpy.dot(y[i:i+len(x)],x)/LA.norm(y[i:i+len(x)]))


        return numpy.argmax(correl)

    def barker_check(self, bits):
        barker = common.get_barker()
        barkerlen = len(barker)
        for i in xrange(len(bits)-barkerlen):
            if (barker == bits[i:i+len(barker)]).all():
                return i
        # else, couldn't find preamble :(
        return -1

    def demap(self, samples):
        bits = numpy.array([])
        softinfo = []
        snr, thresh = self.snr_subsamp(samples)
        print '\t0/1 threshold: %.3f' % thresh
        for s in samples:
            softinfo.append((s-thresh)**2)
            if s > thresh:
                bits = numpy.append(bits, 1)
            else:
                bits = numpy.append(bits, 0)
        return bits, softinfo, snr

    def snr_subsamp(self, samples):

        barker = common.get_barker()
        samp = [ numpy.array([]), numpy.array([]) ]
        for i in xrange(len(barker)):
            samp[barker[i]] = numpy.append(samp[barker[i]], samples[i])
        noise = [0.0, 0.0]
        signal = [0.0, 0.0]
        for i in xrange(2):
            noise[i] = numpy.var(samp[i])
            signal[i] = numpy.mean(samp[i])
        #print 'On 0:', signal[0], noise[0]
        #print 'On 1:', signal[1], noise[1]
        noise = (len(samp[1])*noise[1] + len(samp[0])*noise[0]) / (len(samp[1]) + len(samp[0]))
        thresh = (signal[1] + signal[0])/2
        return (signal[1] - signal[0])**2 / noise, thresh        

    def bin2int(self, bits):
        return int(''.join(map(str, list(bits))),2)

    def decode(self, recd_bits):
        header_code_length = common_txrx.get_header_len()*3
        print('\tFrame Header - Channel Coding result')
        frame_header = self.hamming_decoding(recd_bits[0:header_code_length], 0)
        frame_length = self.bin2int(frame_header[:16])
        index = self.bin2int(frame_header[16:])
        print('\tFrame Payload - Channel Coding result')
        coded_data = recd_bits[header_code_length:frame_length]
        decoded_data = self.hamming_decoding(coded_data, index)
        return decoded_data

    def hamming_decoding(self, coded_bits, index):
        n, k, H = hamming_db.parity_lookup(index)

        decoded_bits = []
        error_count = 0
        for i in range(len(coded_bits)/n):
            block = coded_bits[i*n:(i+1)*n];
            syndrome = numpy.dot(H, block)%2
            if (sum(syndrome)):
                error_count = error_count + 1
                distance = ((H.T + syndrome)%2).sum(1)
                err_ind = ((numpy.nonzero(distance==0))[0])[0]
                #import pdb; pdb.set_trace()
                block[err_ind] = (block[err_ind] + 1)%2
            decoded_bits.extend(block[0:k])
        print '\t\tChannel coding rate: ', k*1.0/(n*1.0)
        print '\t\tBlocks with error count: ', error_count, '/', len(coded_bits)/n
        return decoded_bits 