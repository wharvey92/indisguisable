# audiocom library: Source and sink functions
import common_srcsink as common
import Image
from graphs import *
import binascii
import random
import numpy


class Source:
    def __init__(self, monotone, filename=None, compress=False, encrypt=False):
        '''
        Create a source of the given srctype, with an optional hdr (if
        hdr is True), with numbits payload bits, from an image file if
        specified or according to the srctype.
        self.srctype = {zeroes, carrier, unitstep, random, or imgfile}.
        Return the list of bits (optional header + payload)
        '''
        self.monotone = monotone
        self.fname = filename
        self.compress = compress
        self.encrypt = encrypt
        print 'Source: '

    def process(self):
            if self.fname is not None:
                if self.fname.endswith('.png') or self.fname.endswith('.PNG'):
                    self.srctype = "image"
                    self.srcbits = self.bits_from_image(self.fname)
                else:           # assume it's text
                    self.srctype = "text"
                    self.srcbits = self.text2bits(self.fname)
            else:               
                self.srctype = "monotone"
                self.srcbits = [1]*self.monotone

            if self.encrypt:
                # Generate public key
                pubkey = []
                for i in xrange(len(self.srcbits)):
                    key = random.uniform(0, 1)
                    if key > 0.5:
                        pubkey.append(1)
                    else:
                        pubkey.append(0)

                after_encrypt = []
                for srcbit, keybit in zip(list(self.srcbits), pubkey):
                    after_encrypt.append(srcbit^keybit)
            else:
                after_encrypt = self.srcbits
                pubkey = []

            if self.compress:
                self.stat, self.payload = self.huffman_encode(after_encrypt);    
            else:
                self.stat = {}
                self.payload = after_encrypt

            self.header, self.header_ext = self.get_header(len(self.payload), self.srctype, self.stat)
            self.bits = self.header + self.header_ext + self.payload
            print '\tSource type: ', self.srctype
            print '\tPayload Length: ', len(self.payload)
            print '\tHeader: ', numpy.array(self.header)
            return self.srcbits, self.payload, self.bits, pubkey

    def text2bits(self, filename):
        f = open(filename, 'r')
        bits = []
        line = f.read()
        for c in line:
            bits.extend(self.int2bits(ord(c), 8))
        return bits

    def bits_from_image(self, filename):
        im = Image.open(filename)
        im = im.convert('L')
        d = list(im.getdata())
        return [int(i) for i in "".join(self.bin_nbit(j, 8) for j in d)]

    def bits_from_colored_image(self, filename):
        i = Image.open(filename)
        i = i.convert('RGB')
        d = list(i.getdata())
        return [int(k) for k in "".join("".join(self.bin_nbit(i, 8) for i in j) for j in d)]

    def huffman_encode(self, srcbits):
        block_len = 4
        stat = {}
        offsprings = {}
        codebook = {}
        stat_nz = {}
        symbols = []
        codedbits = []

        # initialize statistics look-up table
        for i in range(2**block_len):
            stat[i] = 0
            codebook[i] = []
            offsprings[i] = [i]

        # Get statistics
        for i in range(len(srcbits)/4):
            sym_temp = self.bits2int(srcbits[i*4:(i+1)*4])
            symbols.append(sym_temp)
            stat[sym_temp] = stat[sym_temp]+1
        
        # Generate Huffman tree
        for i in stat:
            if stat[i]:
                stat_nz[i] = stat[i]
        print stat
        if (len(stat_nz) == 1):
            for sym in stat_nz.keys():
                codebook[sym] = [0]

        else:
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

            # reverse codeword
            for sym_key in codebook:
                codebook[sym_key].reverse()

        # convert to codeword
        for sym in symbols:
            codedbits.extend(codebook[sym])

        print "\tSource bit length: ", len(srcbits)
        print "\tSource-coded bit length: ", len(codedbits)
        print "\tCompression rate: ", (len(codedbits)*1.0)/len(srcbits)

        return stat, codedbits

    def get_header(self, payload_length, srctype, stat): 
        h = list(bin(payload_length)[2:])
        h = [int(x) for x in h]
        if srctype == "monotone":
            header_bit1 = 0
            header_bit2 = 0
        elif srctype == "image":
            header_bit1 = 0
            header_bit2 = 1
        elif srctype == "text":
            header_bit1 = 1
            header_bit2 = 0 
        else:
            header_bit1 = 1
            header_bit2 = 1

        header = [header_bit1] + [header_bit2] + [0]*(16 - len(h)) + h # get it up to 16 bitsi

        # extend header for compression
        header_extend = []
        for sym_key in stat:
            temp = self.int2bits(stat[sym_key], 13)
            header_extend.extend(temp)

        print '\tHeader length: ', len(header) + len(header_extend)

        return header, header_extend
###############
    def bin_nbit(self, number, n):
        out = bin(number)[2:]
        if len(out) > n:
           raise Exception, "number too big to be stored in %d bits: %d" % (n, number)
        while len(out) < n:
           out = '0' + out
        return out

    def bits2int(self, bits):
        out = 0
        for ix in xrange(len(bits)):
            out += bits[ix] * (2**(len(bits) - 1 - ix))
        return int(out)

    def int2bits(self, x, width): 
        return tuple((0,1)[x>>j & 1] for j in xrange(width-1,-1,-1)) 
