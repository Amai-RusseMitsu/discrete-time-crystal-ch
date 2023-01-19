import numpy as np
from qiskit import *
import itertools

def convert_site(site,N):
    # this function perform the site index conversion as follows:
    # 01234, then convert, e.g. 0 to 4, 1 to 3, etc.
    if N%2 == 0:
        return int(N/2-site+1)
    elif N%2 == 1:
        M = (N-1)/2
        return int(M + (M-site))

def spin_combinations(j,N):
    # j is the site index which counts from the l.h.s. in the physical dimension
    # we convert j to ibm_site which counts from the r.h.s.
    combi = list(itertools.product(['0', '1'], repeat=N-1))
    bits = []
    for i in range(2**(N-1)):
        bits.append(''.join(combi[i]))
    spin_up = []
    ibm_site = convert_site(j,N)
    for i in range(2**(N-1)):
        spin_up.append(bits[i][:ibm_site] + '0' +  bits[i][ibm_site:]) # spin up is '0'. previously '1'
    spin_down = []
    for i in range(2**(N-1)):
        spin_down.append(bits[i][:ibm_site] + '1' +  bits[i][ibm_site:]) # spin down is '1'. previously '0'
    return spin_up,spin_down

def expectation_value(N,site,qobj,shot_num,circuit):
    spin_up,spin_down = spin_combinations(site,N)
    spin_up_amp = []
    spin_down_amp = []
    for i in range(2**(N-1)):
        spin_up_amp.append(qobj.get_counts(circuit).get(spin_up[i]))
        spin_down_amp.append(qobj.get_counts(circuit).get(spin_down[i]))
    spin_down_amp = [i/shot_num for i in spin_down_amp if i]
    spin_up_amp = [i/shot_num for i in spin_up_amp if i]
    return sum(spin_up_amp) - sum(spin_down_amp)
