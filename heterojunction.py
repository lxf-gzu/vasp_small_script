#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
@author: XueFeiLiu(201307129@gznu.edu.cn,qq 316187631@qq.com)
@file: heterojunction.py
@time: 2019/1/29 21:44

"""
import re
import sys
import os
from numpy import *
import numpy as np
import pprint

def lattice_write(file):
    lattice=[]
    poscar_tem=[]
    scale=[]
    if os.path.exists(file):
       poscar=open(file,'r')
    
    else:
        raise IOError('POSCAR does not exist！')
    for line in poscar.readlines():
                poscar_tem.append(line)
    poscar.close()
    scale=poscar_tem[1].split()
    lat=poscar_tem[2:5]
    for i in lat:
        lattice.append(i.split())
    scale_lat=[]
    scale_lat.append(scale)
    scale_lat.append(lattice)
    return scale_lat
def other_line(file):
    lines=[]
    poscar_tem=[]
    if os.path.exists(file):
       poscar=open(file,'r')
    
    else:
        raise IOError('POSCAR does not exist！')
    for line in poscar.readlines():
                poscar_tem.append(line)
    line1=poscar_tem[0]
    line6=poscar_tem[5]
    line7=poscar_tem[6]
    if poscar_tem[7].strip().upper().startswith("S"):
       line8=poscar_tem[7]
       line9=poscar_tem[8]
    else:
    
       line8=' '
       line9=poscar_tem[7]    
    lines.append(line1)
    lines.append(line6)
    lines.append(line7)    
    lines.append(line8)
    lines.append(line9)
    return lines    
def atom_coord_write(file):
    cord=[]
    poscar_tem=[]
    atom_num=[]
    if os.path.exists(file):
       poscar=open(file,'r')
   
    else:
        raise IOError('POSCAR does not exist！')
    for line in poscar.readlines():
                poscar_tem.append(line)
    poscar.close()

    atom_n=poscar_tem[6]
    atom_n=atom_n.split()
    atom_number=0
    for i in range(len(atom_n)):

        atom_number+=int(atom_n[i])
    if poscar_tem[7].strip().upper().startswith("S"):    
       coor=poscar_tem[9:9+atom_number]
    else:
       coor=poscar_tem[8:8+atom_number]
    for i in coor:
        
        cord.append(i.split()[0:3])
    reslut=[]
    reslut.append(atom_number)
    reslut.append(cord)
    return reslut

    
if __name__=='__main__':
   print("************By XueFei Liu****************")
   #poscar_a="POSCAR_mullayerBN"
   #poscar_b="POSCAR_mullayerZnO"
   poscar_a= sys.argv[1]
   poscar_b=sys.argv[2]
   scal_lata=lattice_write(poscar_a)
   z_length=float(scal_lata[0][0])*float(scal_lata[1][2][2])
   
   
   scal_latb=lattice_write(poscar_b)
   atom_cordsa=atom_coord_write(poscar_a)
   a_cords=atom_cordsa[1]  
   z_a=[]   
   for i in range(len(a_cords)):
       z_a.append(float(a_cords[i][2]))   
   z_a=list(set(z_a))   
   min_z_a=min(z_a)
   max_z_a=max(z_a)
   
   lay_distance=float(input("Input the layer distance for heterojunction: "))
   lay_distance+=z_length*(max_z_a-min_z_a)
   z_lengtha0=float(scal_lata[1][2][2])
   z_lengthb0=float(scal_latb[1][2][2])
   scal_lata[1][2][2]=str(float(scal_lata[1][2][2])+lay_distance)
   
   
   atom_cordsb=atom_coord_write(poscar_b)
   b_cords=atom_cordsb[1]
   for i in range(len(a_cords)):
       a_cords[i][2]=str(float(a_cords[i][2])*z_lengtha0/float(scal_lata[1][2][2]))
   for i in range(len(b_cords)):
       b_cords[i][2]=str(float(b_cords[i][2])*z_lengthb0/float(scal_lata[1][2][2]))
   lay_dis=round(lay_distance/(float(scal_lata[1][2][2])*float(scal_lata[0][0])),4)
   z_a=[]   
   for i in range(len(a_cords)):
       z_a.append(float(a_cords[i][2]))   
   z_a=list(set(z_a))   
   min_z_a=min(z_a)
   z_b=[]   
   for i in range(len(b_cords)):
       z_b.append(float(b_cords[i][2]))   
   z_b=list(set(z_b))   
   min_z_b=min(z_b)   
   delt_ab_min=min_z_a-min_z_b
   for i in range(len(b_cords)):
       #print(b_cords[i][2])       
       b_cords[i][2]=str(float(b_cords[i][2])+lay_dis+delt_ab_min)
   other_linesa=other_line(poscar_a)
   if len(other_linesa[3]) < 2:
       del other_linesa[3]
   other_linesb=other_line(poscar_b)
   if len(other_linesb[3]) < 2:
       del other_linesb[3]
   #other_linesa[3]='Selective Dynamics\n'  
   atom_ab=other_linesa[1]+other_linesb[1]
   atom_ab=atom_ab.replace('\n','')
   atom_ab=atom_ab+'\n'
   other_linesa[1]=atom_ab
   atom_num_ab=other_linesa[2]+other_linesb[2]
   atom_num_ab=atom_num_ab.replace('\n','')
   atom_num_ab=atom_num_ab+'\n'
   other_linesa[2]=atom_num_ab
   
   new_POSCAR=[]
   new_POSCAR.append(other_linesa[0])
   new_POSCAR.append(" ".join(scal_lata[0])+'\n')
   for i in range(3):
       new_POSCAR.append(" ".join(scal_lata[1][i])+'\n')
       
   for i in range(1,len(other_linesa)):
       new_POSCAR.append(other_linesa[i])
   for i in range(len(a_cords)):
       new_POSCAR.append("   ".join(a_cords[i])+'\n')
   for i in range(len(b_cords)):
       new_POSCAR.append("   ".join(b_cords[i])+'\n')
   poscar_new=open("POSCAR_NEW",'w+')
   for i in new_POSCAR:
       poscar_new.write(i)
   poscar_new.close()
   print("******************************************")
   print("See POSCAR_NEW!")
   
       
   