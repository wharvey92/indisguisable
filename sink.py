# audiocom library: Source and sink functions
import common_srcsink as common 
import Image
from graphs import *
import binascii
import random
import sys

import numpy


class Sink:
    def __init__(self, compression, encrypt):
        self.hlen = common.get_headerlen()
        self.hextlen = common.get_hextlen()
        self.compression = compression
        self.encrypt = encrypt
        print 'Sink:'

    def process(self, recd_bits, pubkey):
        header_bits = recd_bits[:self.hlen]
        srctype, length = self.read_header(header_bits)
        self.srctype = srctype
        print '\tLength from header: ', length
        startdata = self.hlen # data begins at bit offset 27 (= 11 + 16)

        if self.compression:
            header_extension = recd_bits[self.hlen:self.hlen+self.hextlen]
            #print 'Header extenstion: ', len(header_extension)
            stat = self.read_header_ext(header_extension)

            startdata = self.hlen+self.hextlen
            rcd_coded_data = recd_bits[startdata:startdata+length]

            rcd_data = self.huffman_decode(rcd_coded_data, stat)
        else:
            rcd_data = recd_bits[startdata:startdata+length]
            
        if self.encrypt:
            after_decrypt = []
            for databit, keybit in zip(list(rcd_data), pubkey):
                after_decrypt.append(databit^keybit)
            rcd_data = after_decrypt
            
        print '\tRecd', len(rcd_data), 'source bits'
        if self.srctype == "image":
            self.image_from_bits(rcd_data, 'rcd-image.png')
        elif self.srctype == "text":
            print '\tText recd:', self.bits2text(rcd_data)

        return numpy.array(rcd_data, dtype=int)

    def bits2text(self, bits):
        text = []
        intbits = numpy.array([], dtype=numpy.uint8)
        for i in xrange(len(bits)/8):
            intbits = numpy.append(intbits, self.bits2int(bits[i*8:(i+1)*8]))
        for c in intbits:
            text.append(chr(c))
        return  "".join([t for t in text])

    def image_from_bits(self, bits,filename):
        i = Image.new('L',(32,32))
        d = []
        for ix in xrange(32*32):
           b = bits[ix*8 : (ix+1)*8]
           while len(b) < 8:
              b.append(0) # pad with zeros if we got too few bits
           d.append(self.bin_to_int(b))
        i.putdata(d)
        i.save(filename)
        return d

    def colored_image_from_bits(self, bits, filename):
        i = Image.new('RGB', (32,32))
        d = []
        tmp = []
        for ix in xrange(32*32*3):
            b = bits[ix*8 : (ix+1)*8]
            while len(b) < 8:
                b.append(0)
            tmp.append(bin_to_int(b))
            if len(tmp) == 3:
                d.append(tuple(tmp))
                tmp = []
        i.putdata(d)
        i.save(filename)
        return d

    def read_header(self, header_bits): 
        print '\tRecd header: ', numpy.array(header_bits)
        payload_length = self.bits2int(header_bits[2:])
        if header_bits[0] == 0:
           if header_bits[1] == 0:
              srctype = "monotone"
           else:
              srctype = "image"
        else:
           if header_bits[1] == 0:
              srctype = "text"
           else:
              srctype = "NULL"
             
        print '\tSource type: ', srctype
        return srctype, payload_length

    def read_header_ext(self, h_ext):
        stat = {}
        for i in range(len(h_ext)/10):
            stat[i] = self.bits2int(h_ext[i*13:(i+1)*13])
        return stat

    def huffman_decode(self, codedbits, stat):
        block_len = 4
        srcbits = []
        stat_nz = {}
        codebook = {}
        offsprings = {}
        symbook = {}

        # initialize look-up table
        for i in range(2**block_len):
            codebook[i] = []
            offsprings[i] = [i]

        # Generate Huffman tree
        for i in stat:
            if stat[i]:
                stat_nz[i] = stat[i]

        if (len(stat_nz) == 1):
            for sym in stat_nz.keys():
                codebook[sym] = [0]

        while(len(stat_nz)>1):
            # find smallest node and remove
            least_sym = min(stat_nz, key=stat_nz.get)
            least_sym_val = stat_nz[least_sym]
            del stat_nz[least_sym]
            for i in offsprings[least_sym]:
                codebook[i].append(0)
            # find 2nd smallest node and remove
            next_sym = min(stat_nz, key=stat_nz.get)
            # merge two nodes
            stat_nz[next_sym] = stat_nz[next_sym] + least_sym_val
            for i in offsprings[next_sym]:
                codebook[i].append(1)
            offsprings[next_sym].extend(offsprings[least_sym])

        # reverse codeword and generate reverse dictionary
        for sym_key in codebook:
            codebook[sym_key].reverse()

        for sym_key in codebook:
            if codebook[sym_key]:
                # adding 1 in front is a hacky way to save codeword in integers.
                # e.g. preventing [0] to be read as [0,0,0]
                codebook[sym_key] = self.bits2int([1]+codebook[sym_key])
                symbook[codebook[sym_key]] = sym_key

        # convert to srcbits
        index = 0
        while (index < len(codedbits)):
            cw_len = 0
            FLAG = False
            while (not FLAG):
                cw_len = cw_len + 1;
                if (index+cw_len > len(codedbits)):
                    print '\tCannot decode Huffman coded bits'
                    print '\t(Bit sequence cannot be found in codebook)'
                    sys.exit(1)
                codeword = list(codedbits[index:index+cw_len])
                codeword_int = self.bits2int([1]+codeword)
                FLAG = symbook.has_key(codeword_int)
            symbol = symbook[codeword_int]
            srcbits.extend(self.int2bits(symbol, 4))
            index = index + cw_len

        return srcbits
####################################

    def bits2int(self, bits):
        out = 0
        for ix in xrange(len(bits)):
            out += bits[ix] * (2**(len(bits) - 1 - ix))
        return int(out)

    def int2bits(self, number, length):
        temp = [int(x) for x in list(bin(number)[2:])]
        bits = [0]*(length-len(temp)) + temp
        return bits 

    def bin_to_int(self, binary_list):
        out = 0
        for ix in xrange(len(binary_list)):
           out += binary_list[ix] * (2**(len(binary_list) - 1 - ix))
        return out
