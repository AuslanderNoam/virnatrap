import sys
import os
from os.path import isdir,isfile
from .virnatrap import run_virna_pred
import argparse

# Constants ------------------------------------------------------------------------------------------------------------

DESCRIPTION = "Extract viral contigs from a directory with unmapped RNAseq reads fastq files and saves a file with contigs for each fastq in an output directory"
PWD = os.getcwd()

# Terminal functions ---------------------------------------------------------------------------------------------------
def virnatrap_predict():

    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--input", type=str, help="input directory", required=True
    )
    parser.add_argument(
        "--output", type=str, help="output directory", required=True
    )
    parser.add_argument(
        "--fastmode", type=int, help="1 to run in fastmode and call C function to assemble reads", required=False
    )
    parser.add_argument(
        "--multi_proc", type=int, help="run with pool multi processing, if many input files", required=False
    )
    parser.add_argument(
        "--num_threads", type=int,
        help="number of threads to run with pool multi processing",
        required=False
    )
    parser.add_argument(
        "--model_path", type=int, help="path to Tensorflow model to predict whether reads of a fixed length come from viruses or not", required=False
    )

    args = parser.parse_args()

    inpath = args.input
    outpath = args.output

    if not isdir(inpath):
        print('input directory %s not found',inpath)
        sys.exit(1)

    if not isdir(outpath):
        print('output directory %s not found',outpath)
        sys.exit(1)


    if args.fastmode is None:
        fastmode = False
    else:
        fastmode = args.fastmode

    if args.num_threads is None:
        num_threads = 48
    else:
        num_threads = args.num_threads

    if args.multi_proc is None:
        multi_proc = True
    else:
        multi_proc = args.multi_proc


    if args.model_path is None:
        model_path = PWD+'/model/model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5'
    else:
        if not isfile(args.model_path):
            print('model %s not found, using default viRNAtrap model', args.model_path)
            model_path = PWD + '/model/model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5'
        else:
            model_path = args.model_path


    print("Reading fastq at {}...".format(inpath))
    run_virna_pred(inpath, outpath, fastmode, multi_proc, model_path,num_threads)
    print("done")
