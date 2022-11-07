![logo](https://user-images.githubusercontent.com/85499012/154709951-8a742d2b-17ee-4a95-aeab-ea277a698435.png)

viRNAtrap is a package to generate predicted viral contigs from unmapped RNAseq reads



## Overview 
This file describes the software package viRNAtrap [1], an alignment-free method to identify viral reads in RNAseq datasets based on a deep learning model and to assemble predicted viral contigs. There are four key steps in viRNAtrap: 
1. Building a Tensorflow model to predict whether reads of a fixed length come from viruses or not; this has been precomputed based on 48bp reads. 
2. Given one or more input files of RNAseq reads (could be paired or unpaired), map reads to the human genome; this is performed independently by the user.
3. Use the Tensorflow model to predict which of the unmapped reads are viral; this is performed by the included software.
4. Assemble the predicted viral reads into longer contigs; this is performed by the included software using either slow (native Python implementation) or fast (C implementation) modes.

Additional scripts to perform up- or down-stream analyses will be provided upon publication.

Operating systems tested: Linux (CentOS version 7 and Ubuntu 20.04 LTS) and MacOS.

Python versions tested: Python 3.6, 3.7, 3.8 and 3.9

## Citation
[1] Abdurrahman Elbasir, Ying Ye, Daniel E. Schäffer, Jayamanna Wickramasinghe, Paul M. Lieberman, Alejandro A. Schäffer, Noam Auslander. [Characterizing the landscape of viral expression in cancer by deep learning.](https://github.com/AuslanderLab/virnatrap) *In preparation*.

## Installation

Typical install time on a "normal" desktop computer: less than 30 minutes (depending on the number of packages already installed)

### Install with pip
To install the current version of this Github repo, run the following commands
```
git clone https://github.com/AuslanderLab/virnatrap.git
cd virnatrap
pip install .
```

### Install with conda 
1- Get anaconda (64 bit)installer python3.x for linux : https://www.anaconda.com/download/#linux <br />
2- Run the installer : bash Anaconda3-2021.11-Linux-x86_64.sh, and follow the instructions to install anaconda at your preferred directory

#### Run the following commands: <br />
```
git clone https://github.com/AuslanderLab/virnatrap.git
cd virnatrap
conda create --name virnatrap python=3.7 pip
conda activate virnatrap
pip install .
```

#### To deactivate viRNAtrap environment: <br />
```
conda deactivate virnatrap
```

### Python packages and modules used: <br />

Among the publicly available packages or modules used in viRNAtrap are:

* ctypes (https://docs.python.org/3/library/ctypes.html)
* keras (Chollet et al. 2015)
* numpy (Harris et al. 2020)
* tensorflow (Abadi et al. 2015)

Martín Abadi, Ashish Agarwal, Paul Barham, Eugene Brevdo, Zhifeng Chen, Craig Citro, Greg S. Corrado, Andy Davis, Jeffrey Dean, Matthieu Devin, Sanjay Ghemawat, Ian Goodfellow,
Andrew Harp, Geoffrey Irving, Michael Isard, Rafal Jozefowicz, Yangqing Jia, Lukasz Kaiser, Manjunath Kudlur, Josh Levenberg, Dan Mané, Mike Schuster, Rajat Monga, Sherry Moore, Derek Murray, Chris Olah, Jonathon Shlens, Benoit Steiner, Ilya Sutskever, Kunal Talwar, Paul Tucker, Vincent Vanhoucke, Vijay Vasudevan, Fernanda Viégas, Oriol Vinyals, Pete Warden, Martin Wattenberg, Martin Wicke, Yuan Yu, and Xiaoqiang Zheng. TensorFlow: Large-scale machine learning on heterogeneous systems, 2015. Software available from tensorflow.org.

Francois Chollet et al., keras, GitHub. (https://github.com/fchollet/keras) (2015)

Charles R. Harris, K. Jarrod Millman, Stefan J. van der Walt, Ralf Gommers, Pauli Virtanen, David Cournapeau, Eric Wieser, Julian Taylor, Sebastia Berg, Nathaniel J. Smith, Robert Kern, Matti Picus, Stephan Hoyer, Marten H. van Kerkwijk, Matthew Brett, Allan Haldane, Jaime Fernandez del Rio, Mark Wiebe, Pearu Peterson, Pierre Gerard-Marchant, Kevin Sheppard, Tyler Reddy, Warren Weckesser, Hameer Abbasi, Christoph Gohke, Travis E. Oliphant. Array programming with NumPy. Nature 585:357-362 (2020)

## Usage
To use viRNAtrap, a user must provide an input directory continaing one or more input FASTQ files of unmapped reads with file names ending in  *\*_unmapped.fastq*, and a path to an output directory, where a FASTA containing predicted viral contigs will be generated for each FASTQ of unmapped reads in the input directory. For pairs of files with paired reads from the same sample, which may be stored separately for other sequence analysis, the user is advised to concatenate the two files into one combined file for input to viRNAtrap because viRNAtrap treats each distinct input file as if it comes from a distinct sample.


### Running example (Demo):

Expected run time for demo on a "normal" desktop computer: less than 10 minutes

a. To run with an example input fastq file (```input_fastq/example_unmapped.fastq```) run

```
virnatrap-predict --input input_fastq/ --output output_contigs/ 
```

And evaluate the output file generated in ```output_contigs/``` using the expected output in  ```expected_output/output_py.txt```

There is a one-to-one correspondence between input files in directory input_fastq/ and output files in directory output_contigs (or whatever subdirectories the user specifies). If an input file leads to zero predicted viral contigs, then the corresponding output file will be created but will be empty.
The output files are in FASTA format but have the suffix .txt because experience has shown that Mac user prefer the suffix .txt
If one wants to rerun the command with the same  input files and the same output_contigs/ output directory, one should first remove the previous output files. virrnatrap-predict will not overwrite output files that already exist. 

The package comes with a small example that is intended to be used to test if one has installed viRNAtrap correctly. The expected output is in subdirectory expected_output. To test if the above command worked as expected, run the additional command

 ```
 diff expected_output/output_py.txt output_contigs/example_contigs.txt
 ```
 to compare the output in the new installation to the expected output. The installation is correct if the above diff command retruns either no differences or small differences in the less significant digits for the scores in brackets, such as:
 ```
 5c5
 < >contig2[0.8009399]
 ---
 > >contig2[0.8009398]
 ```


b. To run viRNAtrap in fast mode, The fast mode calls a C library to assemble the viral contigs from the model-predicted viral reads. The C library is located at ```virnatrap/src/assemble_read_c.c``` and must first be compiled using the command 

```
gcc -o src/assemble_read_c.so -shared -fPIC -O3 src/assemble_read_c.c
```

The library can also be compiled using an equivalent command for other C compilers.

Then, to run viRNatrap in fast mode, run the command as shown below:

```
virnatrap-predict --input input_fastq/ --output output_contigs/ --fastmode 1 
```

And evaluate the output file generated in ```output_contigs/``` using the expected output in ```expected_output/output_c.txt```

c. It is possible to run viRNAtrap using multiple threads, to process multiple input files in parallel. To run in parrallel using 28 threads run the command as shown below:


```
virnatrap-predict --input input_fastq/ --output output_contigs/ --multi_proc 1 --num_threads 28
```

In multitreaded mode, viRNAtrap will use one thread per file, up to the minimum of the number of available threads and num_threads, where the default num_threads is 48.


### Parameters description:

| Parameter | type | description | default |
| :---: | :---: | :---: | :---: |
| input | path (txt) | path to directory where input fastq is located  | - |
| output | path (txt) | path to directory where output fasta will be generated  | - |
| fastmode | present/absent | run with fast mode (calling C function to assemble viral contigs)  | False (no argument) |
| multi_proc | present/absent | run with multi-processing (if multiple files are in the input directory)  | False (no argument) |
| num_threads | integer | number of threads to use  | Integer (no argument or 48 if multi_proc is True) |
| model_path | path (txt) | path to Tensorflow model to predict viral reads | `/model/model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5` |

## Model

The neural network to predict viral sequnces based on 48bp reads is found in the `model` directory (`model/model_lr_0.005_pool_5_emb_25_l2_0.02_64.hdf5`). A user provided trained Tensorflow model may replace this model, but this will require modification of other functions in virnatrap. 

## Contact

If you have any questions or encounter any difficulties, please create an issue on Github or email us at either aelbasir@wistar.org or nauslander@wistar.org.

