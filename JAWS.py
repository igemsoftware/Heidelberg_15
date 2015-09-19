#!/usr/bin/env python

import sys
import os.path
import random
import RNA
import math
import argparse

def terminal_hairpin(sequence,aptamer,length,shift):
    RNAbet = ['a','c','g','u']
    DNAbet = ['a','c','g','t']
    final_seq = sequence
    for k in range(0,shift):
            final_seq += random.choice(DNAbet)
    final_seq += aptamer
    for k in range(0,shift):
            final_seq += random.choice(DNAbet)
    
    return final_seq

def terminal_hairpin_stem(sequence,aptamer,overhang,shift):
    seq = terminal_hairpin(sequence,aptamer,overhang,shift)
    RNA.pf_fold(seq)
    bp_probability = [[RNA.get_pr(i,j) for i in range(0,len(seq)+1)] for j in range(0,len(seq)+1)]
    Iexists = 0
    IIexists = 0
    pI = 1
    pII = 1
    length = shift + overhang
    seq_length = len(sequence)
    full_length = len(seq)
    for k in range(seq_length-length, seq_length+shift-1):
        if (bp_probability[k+1][full_length-(k-(seq_length-length))] > 0.1) and (bp_probability[k+1][full_length-(k-(seq_length-length))] < 0.99):
            Iexists += 1
            pI *= bp_probability[k][full_length-(k-(seq_length-length))]
    for k in range(seq_length - overhang, seq_length+shift-1):
        if (bp_probability[k+1][full_length-(k-(seq_length-overhang))] > 0.1) and (bp_probability[k+1][full_length-(k-(seq_length-overhang))] < 0.99):
            IIexists += 1
            pII *= bp_probability[k][full_length-(k-(seq_length-overhang))]
            
    #print length
    if (Iexists >= length + shift -2) and (IIexists >= length - 2):
        #if abs(pI-pII) < 0.2:
            print pI
            print pII
            print Iexists
            print IIexists
            return seq
        #else:
        #    return ""
    else:
        return ""
    
def inserted_hairpin(seq_five, seq_three, aptamer, shift):
    RNAbet = ['a','c','g','u']
    DNAbet = ['a','c','g','t']
    final_seq = seq_five
    for k in range(0,shift):
            final_seq += random.choice(DNAbet)
    final_seq += aptamer
    for k in range(0,shift):
            final_seq += random.choice(DNAbet)
    final_seq += seq_three
    return final_seq
    
def inserted_hairpin_stem(seq_five, seq_three, aptamer, shift):
    seq = terminal_hairpin(seq_five,aptamer,0,shift)
    seq_ins = seq + seq_three
    overhang = len(seq_three)
    #print overhang
    RNA.pf_fold(seq_ins)
    bp_probability = [[RNA.get_pr(i,j) for i in range(0,len(seq_ins)+1)] for j in range(0,len(seq_ins)+1)]
    Iexists = 0
    IIexists = 0
    pI = 1
    pII = 1
    length = shift + overhang
    seq_length = len(seq_five)
    full_length = len(seq_ins)
    for k in range(seq_length-length, seq_length+shift-1):
        if (bp_probability[k+1][full_length-(k-(seq_length-length))] > 0.1) and (bp_probability[k+1][full_length-(k-(seq_length-length))] < 0.2):
            Iexists += 1
            pI *= bp_probability[k][full_length-(k-(seq_length-length))]
    for k in range(seq_length - overhang, seq_length+shift-1):
        if (bp_probability[k+1][full_length-(k-(seq_length-overhang))] > 0.1) and (bp_probability[k+1][full_length-(k-(seq_length-overhang))] < 0.2):
            IIexists += 1
            pII *= bp_probability[k][full_length-(k-(seq_length-overhang))]
    if length > seq_length:
        corr_shift = shift
    else:
        corr_shift = 0
    if (Iexists >= length - corr_shift) and (IIexists >= length - shift):
        #if abs(pI-pII) < 0.2:
            return seq_ins
        #else:
        #    return ""
    else:
        return ""
    
import re

def c_iter(ter):
    c_array = []
    for a in ter:
        c_array += [[a.start(),len(a.group())]]
    return c_array

def o_iter(ter):
    o_array = []
    for a in ter:
        o_array += [[a.start(),len(a.group())]]
    return o_array

def build_clozed(array_nested):
    array = []
    for elem in array_nested:
        array = [elem[0]+elem[1]-i for i in range(1,elem[1]+1)] + array
    return array

def build_open(array_nested):
    array = []
    for elem in array_nested:
        array = array + [elem[0]+i for i in range(elem[1])]
    return array

def bond(structure):
    open_are = "(\(\(*)"
    clozed_are = "(\)\)*)"
    open = re.compile(open_are)
    clozed = re.compile(clozed_are)
    citer = c_iter(clozed.finditer(structure))
    oiter = o_iter(open.finditer(structure))
    t_prime = build_clozed(citer)
    f_prime = build_open(oiter)
    return f_prime, t_prime

def full_hairpin(seq_five, seq_three, aptamer, shift):
    RNAbet = ['a','c','g','u']
    DNAbet = ['a','c','g','t']
    final_seq = seq_five
    for k in range(0,shift):
            final_seq += random.choice(RNAbet)
    final_seq += aptamer
    for k in range(0,shift):
            final_seq += random.choice(RNAbet)
    final_seq += seq_three
    return final_seq

def prep_sec1(seq_five, seq_three, seq_apta, shift):
    (take, dump) = RNA.fold(seq_five+"N"+"G"*100+"N"*4+"C"*100+"N"+seq_three)
    seq = len(take.split("."+"("*100+"."*4+")"*100+".",1)[0])*"." +"("*shift+"."*(len(seq_apta))+")"*shift+ len(take.split("."+"("*100+"."*4+")"*100+".",1)[1])*"."
    return seq

def prep_sec2(seq_five, seq_three, seq_apta, shift, rand):
    (take, dump) = RNA.fold(seq_five+"N"+"G"*100+"N"*4+"C"*100+"N"+seq_three)
    (take1, dump) = RNA.fold(seq_apta)
    seq = len(take.split("."+"("*100+"."*4+")"*100+".",1)[0])*"." +"("*rand+take1+")"*(rand)+"."*shift+ len(take.split("."+"("*100+"."*4+")"*100+".",1)[1])*"."
    return seq

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]
    
def full_hairpin(seq_five, seq_three, aptamer, shift):
    RNAbet = ['a','c','g','u']
    DNAbet = ['a','c','g','t']
    take, free_E = RNA.fold(aptamer)
    final_seq = seq_five
    for k in range(0,shift):
            final_seq += random.choice(DNAbet)
    final_seq += aptamer
    for k in range(0,shift):
            final_seq += random.choice(DNAbet)
    final_seq += seq_three
    active_seq = final_seq.split(aptamer)[0]+take.replace(".","N")+final_seq.split(aptamer)[1]
    return final_seq, active_seq

def pseudo_revcomp(sequence):
    seq = sequence.lower()
    base_map = { "a":"u",
               "g":"u",
               "c":"g",
               "u":"a"}
    result = ""
    for elem in seq:
        result += base_map[elem]
    result = result[::-1]
    return result

def full_hairpin_comp(seq_five, seq_three, aptamer, comp, rand):
    RNAbet = ['a','c','g','u']
    DNAbet = ['a','c','g','t']
    take, free_E = RNA.fold(aptamer)
    final_seq = seq_five
    for k in range(0,rand):
            final_seq += random.choice(RNAbet)
    final_seq += aptamer
    for k in range(0,rand-comp):
            final_seq += random.choice(RNAbet)
    final_seq += pseudo_revcomp(seq_five[-comp:])
    final_seq += seq_three
    active_seq = final_seq.split(aptamer)[0]+take.replace(".","N")+final_seq.split(aptamer)[1]
    return final_seq, active_seq

def prep_sec1(seq_five, seq_three, seq_apta, shift):
    (take, dump) = RNA.fold(seq_five+"N"+"G"*100+"N"*4+"C"*100+"N"+seq_three)
    seq = len(take.split("."+"("*100+"."*4+")"*100+".",1)[0])*"." +"("*shift+"."*(len(seq_apta))+")"*shift+ len(take.split("."+"("*100+"."*4+")"*100+".",1)[1])*"."
    return seq

def prep_sec1_comp(seq_five, seq_three, seq_apta, shift):
    (take, dump) = RNA.fold(seq_five+"N"+"G"*100+"N"*4+"C"*100+"N"+seq_three)
    seq = len(seq_five)*"." +"("*shift+"."*(len(seq_apta))+")"*shift+ len(seq_three)*"."
    return seq

def prep_sec2(seq_five, seq_three, seq_apta, shift, rand):
    (take, dump) = RNA.fold(seq_five+"N"+"G"*100+"N"*4+"C"*100+"N"+seq_three)
    seq = len(take.split("."+"("*100+"."*4+")"*100+".",1)[0])*"." +"("*rand+"."*(len(seq_apta)-shift)+")"*(rand)+"."*shift+ len(take.split("."+"("*100+"."*4+")"*100+".",1)[1])*"."
    return seq

def prep_sec2_comp(seq_five, seq_three, seq_apta, comp, rand):
    (take, dump) = RNA.fold(seq_five+"N"+"G"*100+"N"*4+"C"*100+"N"+seq_three)
    seq = (len(seq_five)-comp)*"." +"("*(rand+comp)+"."*(len(seq_apta)-comp)+")"*(rand+comp) + len(seq_three)*"."
    return seq

def prep_sec1_left(seq_five, seq_three, seq_apta, shift):
    (take, dump) = RNA.fold(seq_five+"N"+"C"*100+"N"*4+"G"*100+"N"+seq_three)
    seq = len(seq_five)*"." +"("*shift+"."*(len(seq_apta))+")"*shift+ len(seq_three)*"."
    return seq

def prep_sec2_left(seq_five, seq_three, seq_apta, shift, rand):
    (take, dump) = RNA.fold(seq_five+"N"+"C"*100+"N"*4+"G"*100+"N"+seq_three)
    seq = (len(seq_five)-shift)*"." +"("*(rand+shift)+"."*(len(seq_apta)-shift)+")"*(rand)+")"*shift+ len(seq_three)*"."
    return seq

def diff(a, b):
    b = set(b)
    return [aa for aa in a if aa not in b]
    
def catch0(number):
    if number == 0:
        return 1
    else:
        return number

def pair_entropy(seq,base_idx):
    RNA.pf_fold(seq)
    prob = [[RNA.get_pr(i,j) for i in range(1,len(seq)+1)] for j in range(1,len(seq)+1)]
    s = sum([-elem*math.log(catch0(elem)) for elem in prob[base_idx]])
    return s

def full_hairpin_stem(seq_five, seq_three, aptamer, shift, sec1, sec2):
    global _ACTIVE_TOLERANCE
    global _INACTIVE_TOLERANCE
    global _P_THRESHOLD
    global _MAX_ENERGY_DIFFERENCE
    global _ENERGY_THRESHOLD
    global _ENTROPY
    seq, aseq = full_hairpin(seq_five, seq_three, aptamer, shift)
    seq_ins = seq
    overhang = len(seq_three)
    #print overhang
    dump, a_free_E = RNA.pf_fold(aseq)
    bp_probability = [[RNA.get_pr(i,j) for i in range(1,len(aseq)+1)] for j in range(1,len(aseq)+1)]
    sec1_f, sec1_t = bond(sec1)
    sec2_f, sec2_t = bond(sec2)
    #print bond(sec2)
    Iexists = 0
    IIexists = 0
    #print len(seq_ins)
    pI = (len(seq_ins)*(len(seq_ins)-1))/2-len(sec1_f)*(len(seq_ins)-1)
    #print pI
    pII = (len(seq_ins)*(len(seq_ins)-1))/2-len(sec2_f)*(len(seq_ins)-1)
    #print pII
    length = shift + overhang
    seq_length = len(seq_five)
    full_length = len(seq_ins)
    for k in range(len(sec1_t)):
        #print [sec1_f[k], sec1_t[k]]
        #pI *= bp_probability[sec1_f[k]][sec1_t[k]]
        if (bp_probability[sec1_f[k]][sec1_t[k]] > _P_THRESHOLD):
            #and (bp_probability[sec1_f[k]][sec1_t[k]-1] < 0.01):
            Iexists += 1
            pI *= bp_probability[sec1_f[k]][sec1_t[k]]
    #print "diff "+str(len(sec1_t)-Iexists)
            #print "1ex "+str(Iexists)+"\r"
    dump, b_free_E = RNA.pf_fold(seq_ins)
    bp_probability = [[RNA.get_pr(i,j) for i in range(1,len(seq_ins)+1)] for j in range(1,len(seq_ins)+1)]
    for k in range(len(sec2_t)):
        #print [sec2_f[k], sec2_t[k]]
        pII *= bp_probability[sec2_f[k]][sec2_t[k]]
        if (bp_probability[sec2_f[k]][sec2_t[k]] > _P_THRESHOLD):
            IIexists += 1
            #pII *= bp_probability[sec2_f[k]][sec2_t[k]]
    #print pII, IIexists
            #print pII
            #print IIexists
    #if length > seq_length:
    #    corr_shift = shift
    #else:
    #    corr_shift = 0
    #print pI, pII
    S = sum([pair_entropy(seq,i) for i in range(len(seq))]) 
    if (IIexists >= len(sec2_f)-_INACTIVE_TOLERANCE) and (Iexists >= len(sec1_f)-_ACTIVE_TOLERANCE) and (b_free_E < _ENERGY_THRESHOLD) and (S < _ENTROPY):
        #pI /= Iexists
        #pII /= IIexists
        if abs(a_free_E-b_free_E) < _MAX_ENERGY_DIFFERENCE:
            print "Free energy active: "+str(a_free_E), "Free energy inactive: "+str(b_free_E), "Free energy difference: "+str(a_free_E-b_free_E)
            print "Entropy: "+str(S)
	    print "Shift: "+str(shift)
            print "Active state base pairs: "+str(Iexists)
            print "Inactive state base pairs: "+str(IIexists)
            print "Sequence: "+str(seq_ins)
            print "good sequence"
            return seq_ins
        else:
            print "Free energy active: "+str(a_free_E), "Free energy inactive: "+str(b_free_E), "Free energy difference: "+str(a_free_E-b_free_E)
            print "Shift: "+str(shift)
            print "Active state base pairs: "+str(Iexists)
            print "Inactive state base pairs: "+str(IIexists)
            print "Sequence: "+str(seq_ins)
            print "less-than-good sequence"
            return seq_ins
        #else:
        #    return ""
    else:
        return ""

parser = argparse.ArgumentParser(description='JAWS - Joining Aptamers Without SELEX', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('fprime', metavar='5PRIME', help="5' sequence of the DNAzyme/ribozyme")
parser.add_argument('aptamer', help='aptamer sequence')
parser.add_argument('tprime', metavar='3PRIME', help="3' sequence of the DNAzyme/ribozyme")

parser.add_argument('-p', '--p-threshold', type=float, dest='P_THRESHOLD', default=1e-7, help='Minimal probability of base pairing at which a base-pair is considered to be present')
parser.add_argument('-a', '--active-tolerance', type=int, dest='ACTIVE_TOLERANCE', default=10, help='Number of nucleotides not conforming to the active structure to be tolerated')
parser.add_argument('-i', '--inactive-tolerance', type=int, dest='INACTIVE_TOLERANCE', default=2, help='Number of nucleotides not conforming to the inactive structure to be tolerated')
parser.add_argument('-e', '--energy-threshold', type=float, dest='ENERGY_THRESHOLD', default=-12, help='Maximal energy of structures allowed')
parser.add_argument('-m', '--max-energy-difference', type=float, dest='MAX_ENERGY_DIFFERENCE', default=9, help='Maximal energy difference between active and inactive state')
parser.add_argument('-s', '--stem-length', type=int, dest='STEM_LENGTH', default=10, help='Length of the stem connecting the aptamer to the DNAzyme or ribozyme')
parser.add_argument('-f', '--shift', type=int, dest='SHIFT', default=7, help='Number of nucleotides to be displaced upon binding of the ligand')
parser.add_argument('-r', '--params', dest='PARAMS', default='dna_mathews2004.par', help='ViennaRNA parameter set. Can be one of dna_matthews1999.par, dna_matthews2004.par, rna_turner1999.par, rna_turner2004.par, or rna_andronescu2007.par to use one of the parameter sets included with ViennaRNA or a path to a custom parameter set')
parser.add_argument('-v', '--vienna-path', dest='VIENNA_PATH', default='/usr/share/ViennaRNA', help='Directory containing ViennaRNA parameter files. Required only if using a parameter set included with ViennaRNA')
parser.add_argument('-y', '--entropy', type=float, dest='ENTROPY', default=10000)

args = parser.parse_args()


_P_THRESHOLD = args.P_THRESHOLD
_ACTIVE_TOLERANCE = args.ACTIVE_TOLERANCE
_INACTIVE_TOLERANCE = args.INACTIVE_TOLERANCE
_ENERGY_THRESHOLD = args.ENERGY_THRESHOLD
_MAX_ENERGY_DIFFERENCE = args.MAX_ENERGY_DIFFERENCE
_STEM_LENGTH = args.STEM_LENGTH
_SHIFT = args.SHIFT
_ENTROPY = args.ENTROPY

_5PRIME = args.fprime.upper()
_APTAMER = args.aptamer.upper()
_3PRIME = args.tprime.upper()

paramspath = os.path.join(args.VIENNA_PATH, args.PARAMS)
if not os.path.exists(paramspath):
    paramspath = args.PARAMS

RNA.read_parameter_file(paramspath)

print _5PRIME 
print _APTAMER 
print _3PRIME 
for i in range(100000):
    shift0 = _STEM_LENGTH
    shift1 = _SHIFT
    #random.choice(range(2,shift0))
    a = full_hairpin_stem(_5PRIME , _3PRIME , _APTAMER ,shift0,prep_sec1_left(_5PRIME , _3PRIME , _APTAMER ,shift0),prep_sec2_left(_5PRIME , _3PRIME , _APTAMER ,shift1,shift0))
    if a != "":
        print a

