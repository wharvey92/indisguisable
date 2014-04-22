import math
import common_txrx as common
import numpy

import hamming_db
 
class Transmitter:
    def __init__(self, carrier_freq, samplerate, one, spb, silence, cc_len):
        self.fc = carrier_freq  # in cycles per sec, i.e., Hz
        self.samplerate = samplerate
        self.one = one
        self.spb = spb
        self.silence = silence
        self.cc_len = cc_len
        print 'Transmitter: '
    def add_preamble(self, databits):
        '''
        Prepend the array of source bits with preamble bits
        '''
        databits_with_preamble = []
        preamble = common.get_barker()
        preamble_repeated = []
        for x in preamble:
            preamble_repeated = numpy.append(preamble_repeated, [x]*common.get_rep())
        preamble_with_silence = numpy.append([0]*self.silence, preamble_repeated)
        databits_with_preamble = numpy.append(preamble_with_silence, databits)
        print '\tSent Preamble: ', numpy.array(common.get_barker())
        return databits_with_preamble


    def bits_to_samples(self, databits_with_preamble):
        samples = common.bits2samples(databits_with_preamble, self.one, self.spb)
        return samples
        

    def modulate(self, samples):
        '''
        Multiply samples by a local sinusoid carrier of the same length.
        Return the multiplied result.
        '''
        print '\tNumber of samples being sent:', len(samples)
        return samples * common.local_carrier(self.fc, len(samples), self.samplerate)

    def encode(self, databits):
        index, coded_bits = self.hamming_encoding(databits)
        frame_length = len(coded_bits) + common.get_header_len()*3
        header = numpy.array(self.int2bin(frame_length,16) + self.int2bin(index,4))
        temp, coded_header = self.hamming_encoding(header, True)
        return coded_header + coded_bits

    def hamming_encoding(self, databits, is_header=False):
        if is_header:
            n, k, index, G = hamming_db.gen_lookup(3)
        else:
            n, k, index, G = hamming_db.gen_lookup(self.cc_len)

        coded_bits = []

        if len(databits) % k:
            databits = numpy.append(databits, [0]*(k - len(databits)%k))
        #import pdb; pdb.set_trace()
        for i in range(len(databits)/k):
            coded_bits.extend(numpy.dot(databits[i*k:(i+1)*k],G)%2)

        return index, coded_bits

    def int2bin(self, num, length):
        num_bin = [int(x) for x in list(bin(num)[2:])]
        return [0]*(length - len(num_bin)) + num_bin