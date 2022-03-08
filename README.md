![logo](https://user-images.githubusercontent.com/85499012/154709951-8a742d2b-17ee-4a95-aeab-ea277a698435.png)

viRNAtrap is a package to generate predicted viral contigs from unmapped RNAseq reads



## Overview 
This file describes the package to implement viRNAtrap, an alignment-free method to identify viral reads based on a deep learning model, and assemble predicted viral contigs.
Additional scripts to perform up or down-stream analyses are provided in [TBD].
Operating systems tested: Linux and MacOS,
Python versions tested:  python3.6, 3.7, 3.8 and 3.9

## Citation
[TBD]

## Installation

pip 
conda

## Usage

To use virnatrap with an input fastq file of unmapped reads, provide a path to an input directory where the input fastq are located (preferably unmapped reads will be named  *_unmapped.fastq)*, and a path to an output directory, where a fasta containing predicted viral contigs will be generated for each fastq in the input directory  

### Running example:
a. To run with an example input fastq file (```input_fastq/example_unmapped.fastq```) run
Using Linux or Windows

```
./virnatrap --input $PWD/../input_fastq/ --output $PWD/../output_contigs/
```

Using Mac OS

```
python3 virnatrap.py --input $PWD/../input_fastq/ --output $PWD/../output_contigs/
```

And evaluate the output file generated in ```output_contigs/``` using the expected output in ```expected_output/output_py.txt```

b. Run this command line bellow to run in fast mode. The fast mode calls a C function to assemble the viral contigs from the model-predicted viral reads. The C code is located at ```virnatrap/src/assemble_read_c.c``` and is compiled using 
```
gcc -o src/assemble_read_c.so -shared -fPIC -O3 src/assemble_read_c.c
```

Using Linux or Windows

```
python3 virnatrap.py --fastmode 1 --input $PWD/../input_fastq/ --output $PWD/../output_contigs/
```

Using Mac OS

```
./virnatrap --fastmode 1 --input $PWD/../input_fastq/ --output $PWD/../output_contigs/
```

And evaluate the output file generated in ```output_contigs/``` using the expected output in ```expected_output/output_c.txt```

### Parameters description:

| Parameter | type | description | default |
| :---: | :---: | :---: | :---: |
| input | path (txt) | path to directory where input fastq is located  | - |
| output | path (txt) | path to directory where output fasta will be generated  | - |
| fastmode | binary | run with fast mode (calling C function to assemble viral reads)  | False |
| multi_proc | binary | run with multi-processing (if multiple files are in the input directory)  | True |


## Model

The neural network to predict viral sequnces based on 48bp reads is found in the `model` directory (`model/model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5`)

## Contact

If you run into issues, please an issue on Github or email us at either aelbasir@wistar.org or nauslander@wistar.org

