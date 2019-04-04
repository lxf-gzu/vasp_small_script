#-*- coding: UTF-8 -*-
import os
import numpy as np
os.system("grep 'including local field effects in DFT' OUTCAR  -A 5 |tail -n 5 > dielect.dat")
os.system(" sed -i '1d' dielect.dat ")
os.system(" sed -i '$d' dielect.dat ")
os.system("sed -i 's/^[ \t]*//g' dielect.dat") #del whitespace of head of each line
file=open('dielect.dat')
lines=file.readlines()
die=[]
for line in lines:
    tem=line.split()
    die.append(tem)
die=np.array(die)
die_harmonic_mean=((float(die[0][0])**(-1)+float(die[1][1])**(-1)+float(die[2][2])**(-1))/3)**(-1)
die_harmonic_mean=round(die_harmonic_mean,2)

os.environ['die_harmonic_mean']=str(die_harmonic_mean)
os.system("echo $die_harmonic_mean >die_harmonic_mean.dat")



