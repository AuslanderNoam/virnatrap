"""
Model classes and functions to identify viral reads and assemble viral contigs from input fastq
"""
# Imports --------------------------------------------------------------------------------------------------------------
import os
import re
import sys
from multiprocessing import Pool,freeze_support
import multiprocessing as mp
import numpy as np
from pkg_resources import resource_filename
from collections import OrderedDict
import random
import tensorflow as tf
import ctypes
from ctypes import *
from tensorflow import get_logger, autograph
from tensorflow import keras
from tensorflow.keras import Sequential
from tensorflow.keras.layers import *
from tensorflow.keras.models import load_model
from tensorflow.keras.preprocessing.sequence import pad_sequences
import glob
import numpy as np
import argparse
from os.path import exists

DESCRIPTION = "Extract viral contigs from a directory with unmapped RNAseq reads fastq files and saves a file with contigs for each fastq in an output directory"

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


get_logger().setLevel('ERROR')
autograph.set_verbosity(0)


flatten = lambda l: [item for sublist in l for item in sublist]

PWD = os.getcwd()

# Constants ------------------------------------------------------------------------------------------------------------
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "T", "C", "G"])}
NUCLEOTIDES = sorted([x for x in DEFAULT_NUC_ORDER.keys()])
SEGMENT_LENGTH = 48
SEARCHSUBLEN = 24

# Functions ------------------------------------------------------------------------------------------------------------

def random_base():
    """
    Generate a random base.
    :return: Random base.
    """
    return random.choice(NUCLEOTIDES)


def handle_non_ATGC(sequence):
    """
    Handle non ATGCs.
    :param sequence: String input.
    :return: String output (only ATCGs), with randomly assigned bp to non-ATGCs.
    """
    ret = re.sub('[^ATCG]', random_base, sequence)
    assert len(ret) == len(sequence)
    return ret


def pad_sequence(sequence, source_sequence, length=SEGMENT_LENGTH):
    """
    Pad sequence by repeating it.
    :param sequence: Segmented sequence to pad.
    :param source_sequence: Original, complete sequence.
    :param length: Length of sequence to pad up to.
    :return: Padded sequence of lenth length.
    """
    assert len(sequence) < length, len(sequence)
    assert sequence == source_sequence or len(source_sequence) > length
    if len(source_sequence) > length:
        ret = source_sequence[len(source_sequence)-len(sequence)-(length-len(sequence)):
                              len(source_sequence)-len(sequence)]+sequence

    else:
        assert sequence == source_sequence
        ret = (source_sequence * (int(length / len(sequence)) + 1))[:length]
    assert len(ret) == length
    return ret


def encode_sequence(sequence, nuc_order=None):
    """
    Encode a sequence to integers for use in Python LSTM model.
    :param sequence: Sequence to encode.
    :param nuc_order: Order of nucleotides for encoding.
    :return: Encoded sequence as integers.
    """
    if nuc_order is None:
        nuc_order = DEFAULT_NUC_ORDER

    sequence = sequence.upper()[:SEGMENT_LENGTH]
    accepted_nucleotides = "".join(nuc_order.keys())

    assert re.match('^[{}]+$'.format(accepted_nucleotides), sequence) is not None, \
        "Only {} allowed".format(accepted_nucleotides)

    encoded_seq = [(nuc_order[x] + 1) for x in sequence]
    encoded_seq = np.array([encoded_seq])

    return encoded_seq


def encode_sequences(sequences, nuc_order=None, segment_length=SEGMENT_LENGTH):
    """
    Encode a sequence to integers for use in model.
    :param sequences: List of sequences to encode.
    :param nuc_order: Order of nucleotides for encoding.
    :param segment_length: Segments should be at most this length.
    :return: Encoded sequence as integers.
    """
    encoded_seqs = []
    for sequence in sequences:
        encoded_seqs.append(encode_sequence(sequence, nuc_order)[0])

    return np.array(pad_sequences(encoded_seqs, maxlen=segment_length, padding="post"))


def save_model(model,model_path):
    '''
    saves Keras model : lets save in a specific directory for saved models
    :param model:
    :param path_to_save:
    '''
    model.save(model_path+'.h5')

def load_model_keras(model_path):
    """
    Loads Keras model.
    :param model_path: Path to H5 model.
    :return: Keras model.
    """
    model_loaded = load_model(model_path)
    return model_loaded


def assemble_right(read, read_list,score_list,score_read = 1,sc_thr = 0.5,runs=5000,sublen=SEARCHSUBLEN):
    '''
    goes forward
    :param read:
    :param read_list:
    :param score_list:
    :param score_read:
    :param sc_thr:
    :param runs:
    :param sublen:
    :return:
    '''
    scl = [score_read]
    contig = read
    sb0 = read[len(read)-sublen:]
    flag = True
    cnt = 0
    while flag and cnt<runs:

        ids = [i for i in range(len(read_list)) if sb0 in read_list[i]]
        strings_with_substring = [read_list[i] for i in ids]
        score0 = [score_list[i] for i in ids]
        pos = [sb.find(sb0) for sb in strings_with_substring]
        if len(strings_with_substring)>0 and min(pos) < len(read)-sublen:
            newsub = strings_with_substring[np.argmin(pos)][min(pos)+len(sb0):]
            contig =  contig + newsub
            sb0 = strings_with_substring[np.argmin(pos)][len(read)-sublen:]
            scl+= [score0[np.argmin(pos)]]
            cnt += 1
            read_list = [read_list[i] for i in range(len(read_list)) if i not in ids]
            score_list = [score_list[i] for i in range(len(score_list)) if i not in ids]

            if np.mean(scl)<sc_thr:
                flag = False

        else:
            flag = False

    return contig,read_list,score_list,scl


def assemble_left(read, read_list, score_list, scores_read, sc_thr=0.5, runs=5000, sublen=SEARCHSUBLEN):
    '''
    goes backwards
    :param read:
    :param read_list:
    :param score_list:
    :param scores_read:
    :param sc_thr:
    :param runs:
    :param sublen:
    :return:
    '''
    scl = scores_read
    contig = read
    sb0 = read[:sublen]
    flag = True
    cnt = 0
    while flag and cnt < runs:
        ids = [i for i in range(len(read_list)) if sb0 in read_list[i]]
        strings_with_substring = [read_list[i] for i in ids]
        score0 = [score_list[i] for i in ids]
        pos = [sb.find(sb0) for sb in strings_with_substring]
        if len(strings_with_substring) > 0 and max(pos) > 0:
            newsub = strings_with_substring[np.argmax(pos)][:max(pos)]
            contig = newsub + contig
            sb0 = strings_with_substring[np.argmax(pos)][:sublen]
            scl += [score0[np.argmax(pos)]]
            cnt += 1
            read_list = [read_list[i] for i in range(len(read_list)) if i not in ids]
            score_list = [score_list[i] for i in range(len(score_list)) if i not in ids]
            if np.mean(scl) < sc_thr:
                flag = False

        else:
            flag = False

    return contig, read_list, score_list, scl


def assemble_read(read,read_list,score_list,score_read):
    '''

    :param read: seed read to star with
    :param read_list: all reads to search for assembling
    :param score_list: model scores per read
    :param score_read: model score of a read
    :return: assembled contig, new read and score list (removes reads already used) and the mean pred score of the contig
    '''
    contigr,read_list,score_list,scr = assemble_right(read, read_list,score_list,score_read)
    contigl,read_list,score_list,scl = assemble_left(read, read_list,score_list,scr)
    return contigl.replace(read, '') + contigr,read_list,score_list,np.mean(scl)


def assemble_read_loop(readsv, reads0, scores0, scoresv, filen, lenthr = 48):

    conts = []

    cnts=0

    f = open(filen, 'w')
    while len(readsv)>0:
        contig,reads0,scores0,msc = assemble_read(readsv[0], reads0,scores0,scoresv[0])

        if len(contig)>=lenthr and msc>0.5:
            f.write('>'+'contig'+str(cnts)+'['+str(msc)+']'+'\n'+str(contig.replace('\n',''))+'\n')
            conts.append(contig)
            cnts+=1

        incv = [i for i in range(len(readsv)) if readsv[i] not in contig]
        readsv = [readsv[i] for i in incv]
        scoresv = [scoresv[i] for i in incv]

    f.close()
    return 1


def load_virus_model(model_path):
    model = load_model_keras(model_path)
    return model

def filter_sequences(seqs):
    seq_remove = ['A'*SEARCHSUBLEN,'C'*SEARCHSUBLEN,'G'*SEARCHSUBLEN,'T'*SEARCHSUBLEN]
    seqs = [i for i in seqs if seq_remove[0] not in i and seq_remove[1] not in i and seq_remove[2] not in i and seq_remove[3] not in i]
    return seqs

def proc_fastq(infile):
    f = open(infile)
    l = f.readlines()
    x = [i for i in range(1, len(l), 4)]
    seqs = [handle_non_ATGC(l[i].replace('\n', '')) for i in x]
    seqs = list(np.unique([handle_non_ATGC(i) for i in seqs]))
    seqs = filter_sequences(seqs)
    medsize = np.median([len(i) for i in seqs])
    seqs = [i[:int(medsize)] for i in seqs]
    encoded_c = encode_sequences(seqs)
    f.close()
    return encoded_c,seqs

def make_clist(lst):
    return (c_char_p * len(lst))(*[x.encode() for x in lst])

def assemble_read_call_c(readsv,reads0,scores0,scoresv,filen):
    librd = CDLL(PWD+"/src/assemble_read_c.so")
    librd.connect()

    arr_f = (ctypes.c_float * len(scores0))(*list(scores0))
    arr_fv = (ctypes.c_float * len(scoresv))(*list(scoresv))
    filen = c_char_p(filen.encode())
    numelv = len(scoresv)
    arr_ch = make_clist(reads0)
    arr_chv = make_clist(readsv)
    m = c_int(len(reads0))
    n = c_int(len(readsv[0]))
    librd.assemble_read.restype = ctypes.c_char_p
    result = librd.assemble_read_loop(arr_f,arr_fv, arr_ch,arr_chv, n, m, numelv,filen)

    return result



def extract_contigs(invars,large_file_thr=1000000):
    '''

    :param invars:
    :param lenthr:
    :param large_file_thr:
    :return:
    '''

    inpath = invars[0]
    outpath = invars[1]
    fastmode = invars[2]
    model_path = invars[3]

    # output file name
    fn = outpath + inpath.split('/')[-1].replace('_unmapped', '_contigs').replace('fastq', 'txt')
    if '_contigs' not in fn:
        fn = fn.replace('.txt', '_contigs.txt')

    file_exists = exists(fn)
    if file_exists:
        return 0
    else:

        ##Thats the trained neural network model
        model = load_virus_model(model_path)

        #Encode unpapped RNAseq reads for model
        encoded_c,seqs = proc_fastq(inpath)

        # Predict the viral unpapped RNAseq reads using the model
        scc = model.predict(encoded_c)

        scores = list(scc)

        # Select seeds for assembling contigs - reads scored more than 0.7 by model
        qq = 0.7
        vv = [seqs[i] for i in range(len(scores)) if scores[i] > qq]
        scv = [scores[i] for i in range(len(scores)) if scores[i] > qq]

        if fastmode:
            assemble_read_call_c(vv, seqs, scores, scv, fn)
        else:
            assemble_read_loop(vv, seqs, scores, scv, fn)

        return 1

def run_virna_pred(inpath,outpath,fastmode,multi_proc,model_path,num_threads):


    ##get input fastq files
    infastq = list(set([i for i in glob.glob(inpath + '*.fastq') if '.fastq' in i]))

    ##To make sure the outputs were not yet generated, get the input fastq that do not have output contigs yet
    outs = glob.glob(outpath + '/' + '*.txt')
    nmi = [i.split('/')[-1].replace('.fastq', '').replace('_unmapped', '') for i in infastq]
    nmo = [i.split('/')[-1].replace('.txt', '').replace('_contigs', '') for i in outs]
    infiles = [infastq[i] for i in range(len(nmi)) if nmi[i] not in nmo]

    print('starting_prediction...')
    if multi_proc:
        freeze_support()
        pool = mp.Pool(processes=num_threads)
        pool.map(extract_contigs, [[f, outpath, fastmode,model_path] for f in infiles])

    else:
        ins = [[f, outpath, fastmode,model_path] for f in infiles]
        for i in range(len(ins)):
            extract_contigs(ins[i])

    print("Done processing")


