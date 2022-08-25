# Imports --------------------------------------------------------------------------------------------------------------
import glob
import os

input_dir = os.getcwd()


files0 = glob.glob(input_dir+'/*.fastq')
files = [i for i in files0 if '.1.fastq' in i or '.2.fastq' in i]

uf = list(set([i.split('.')[0] for i in files]))

for fs in uf:
    fls = [i for i in files if fs in i]
    if len(fls)>1:
        os.system("cat "+fls[0]+" "+fls[1]+" > "+fs+".fastq")
        os.system("rm -rf " + fls[0] + " " + fls[1])
