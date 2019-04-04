#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import sys

#Getting the bond length and bond engle Written by Gang Tang
### Usage Commands ###

'''
In terminal, type: 

python dos_angle_v2.py  CONTCAR  or

python dos_angle_v2.py  POSCAR  
'''

print """
----------------------------------------------------------
Please modify the atom information in the code before use,

     and only supports POSCAR with "direct" format!
----------------------------------------------------------
"""

#Read CONTCAR, please do not change thses!
in_file = sys.argv[1]
file_poscar=open(in_file,'r')
line=file_poscar.readlines()
a1 = float(line[2].split()[0])
a2 = float(line[3].split()[0])
a3 = float(line[4].split()[0])
b1 = float(line[2].split()[1])
b2 = float(line[3].split()[1])
b3 = float(line[4].split()[1])
z1 = float(line[2].split()[2])
z2 = float(line[3].split()[2])
z3 = float(line[4].split()[2])
vector_a=np.array([a1,b1,z1])
vector_b=np.array([a2,b2,z2])
vector_c=np.array([a3,b3,z3])
vector_all=np.array([vector_a,vector_b,vector_c])

#Read the atomic coordinates. Note that you need to change these (e.g., atom_YY and line[XX]) based on your own system. 
#For example, the 1st atom, line[8]; the 2st atom, line[9] and so on (YY +7 = XX in line[XX])
#atom 1 e.g.,C
atom_1=np.array([float(line[8].split()[0]),float(line[8].split()[1]),float(line[8].split()[2])])
atom_1_C=np.dot(atom_1,vector_all)    #C represents the Cartesian coordinate
#atmo 5 e.g.,Si
atom_5=np.array([float(line[12].split()[0]),float(line[12].split()[1]),float(line[12].split()[2])])
atom_5_C=np.dot(atom_5,vector_all)
#atom 8 e.g.,Si
atom_8=np.array([float(line[15].split()[0]),float(line[15].split()[1]),float(line[15].split()[2])])
atom_8_C=np.dot(atom_8,vector_all)


#bond length, note that you need to change these (e.g., bong_xx_yy and atom_zz_C) based on your own system.
bond_5_1=np.sqrt(np.sum(np.square(atom_5_C-atom_1_C)))  #bond length between atom 5 and atom 1.
bond_8_1=np.sqrt(np.sum(np.square(atom_8_C-atom_1_C)))

#bond angle, note that which is the top of angle! 
#and you need to change these (e.g., bong_xx_yy, atom_zz_C, and vector_xx_yy) based on your own system.
vector_5_1=np.array(atom_5_C-atom_1_C)
vector_8_1=np.array(atom_8_C-atom_1_C)
cos_angle=vector_5_1.dot(vector_8_1)/(bond_5_1*bond_8_1)  #bond angle among atom 5,atom 1, and atom 8. And atom 1 is the top of the angle!
angle=np.arccos(cos_angle)*360/2/np.pi  #unit conversion

#you need to change these (e.g., bong_xx_yy) based on your own system.
print bond_5_1,bond_8_1,angle 
#output
#out_file = sys.argv[2]
#output=open(out_file,'w')
#output.write(str(bond_7_10)+'\t'+str(bond_3_10)+'\t'+str(bond_1_6)+'\t'+str(angle)+'\n')
#output.close()

file_poscar.close()



