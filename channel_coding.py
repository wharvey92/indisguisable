'''
Channel coding module for transmitter and receiver
'''

import numpy

import hamming_db   # Matrices for Hamming codes

''' Transmitter side ---------------------------------------------------
'''
def get_frame(databits, cc_len):
    '''
    Perform channel coding to <databits> with Hamming code of n=<cc_len>
    Add header and also apply Hamming code of n=3
    Return the resulted frame (channel-coded header and databits)
    '''
    
    return frame
    
def get_header(payload, index):
    '''
    Construct and return header for channel coding information.
    Do not confuse this with header from source module.
    Communication system use layers and attach nested headers from each layers 
    '''

    return header
    
def encode(databits, cc_len):
    '''
    Perfrom channel coding to <databits> with Hamming code of n=<cc_len>
    Pad zeros to the databits if needed.
    Return the index of our used code and the channel-coded databits
    '''
    
    return index, coded_bits

''' Receiver side ---------------------------------------------------
'''    
def get_databits(recd_bits):
    '''
    Return channel-decoded data bits
    Parse the header and perform channel decoding.
    Note that header is also channel-coded    
    '''

    return databits

def decode(coded_bits, index):
    '''
    Decode <coded_bits> with Hamming code which corresponds to <index>
    Return decoded bits
    '''

    return decoded_bits
