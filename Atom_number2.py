#!/bin/python
# -*- coding:utf-8 -*-
"""
@author: XueFeiLiu(201307129@gznu.edu.cn,qq 316187631@qq.com)
@file: Atom_number.py
@time: 2019/1/29 21:44
A script to obtain the atoms in POSCAR to extract data of PDOS with vaspkit -task 114
"""
import os
from numpy import *
import numpy as np
import pprint
def atom_write():
    cord=[]
    poscar_tem=[]
    atom_num=[]
    if os.path.exists("POSCAR"):
       poscar=open("POSCAR",'r')
    elif os.path.exists("CONTCAR"):
         poscar=open("CONTCAR",'r')
    else:
        raise IOError('CONTCAR OR POSCAR does not exist！')
    for line in poscar.readlines():
                poscar_tem.append(line)
    poscar.close()

    atom_n=poscar_tem[6]
    atom_n=atom_n.split()
    atom_number=0
    for i in range(len(atom_n)):

        atom_number+=int(atom_n[i])
    coor=poscar_tem[8:8+atom_number]
    for i in coor:
        cord.append(i.split()[0:3])
    reslut=[]
    reslut.append(atom_number)
    reslut.append(cord)
    return reslut
    
	
def write_atom_acording_layer():
    
     data=atom_write()
    
     atom_number=data[0]
     cord=data[1]
     key=list(range(1,atom_number+1))
     value=[]
     for i in range(atom_number):
         value.append(cord[i][2])
     layer=dict(zip(key,value))
     layer=sorted(layer.items(),key=lambda x:x[1])
     return layer
	
if __name__  == '__main__':
   
   print("You may contact 316187631@QQ.com if any questions")
   print(" ***************************************************")
   print("There are two methods to select atoms,one based on a Z coordination and another based on ")
   print(" layer order numbers. please enter '1' for method 1 and '2' for method 2 ")
   method=int(input("Input the method you choose: "))
   if method == 1:
      print("****************************************************************************************")
      print("You need firstly get a Z coordinate in a layer in POSCAR, e.g. 0.308030 ")
      print("Then,you need input a 'delt' value to include all the coordinates in a same layer ")
      print("e.g. 0.00001,it means all the atoms in a layer will be selected if Z[i]-0.308030 < delt ")
      print("****************************************************************************************")
      print("I will give a list to help you select a atom : " )
      print(" \n ")
      atom_sort=write_atom_acording_layer()
      pprint.pprint(atom_sort)
      data=atom_write()
      atom_number=data[0]
      cord=data[1]
      c_atom=float(input("Input an arbitrary coordinate in correspond layers: "))
      delt=float(input("Input the threshold value to seleclt a atom in a layer, e.g. 0.002 :"))
      atom=open("atom.dat",'w+')
      for i in range(atom_number):
           if abs(float(cord[i][2])-c_atom ) < delt:
              atom.write(str(i+1)+' ')
      atom.write("\n")

      atom.close()
      print("**********************************")
      print(" The atom data has been written in 'atom.dat' ")
      print(" The PDOS dat has been produced with vaspkit if it works ")
      print(" The log file by vaspkit has been write in 'log_vaspkit'" )
      os.system(" echo '11' >> select.dat ")
      os.system(" echo '114' >> select.dat")
      os.system(" cat atom.dat >> select.dat")
      os.system(" vaspkit < select.dat >log_vaspkit")
   elif method == 2:
        layer=write_atom_acording_layer()
        pprint.pprint(layer)
        data=atom_write()
        atom_number=data[0]
        delt=0.01
        k=1
        lay=[]
        lay_n=[]
        tmp=[]
        for i in range(1,atom_number):
            if abs(float(layer[i][1])-float(layer[i-1][1]))> delt:
               k+=1
               lay.append(i)
        lay.append(atom_number)
        key0=list(range(k))
        lay_num=dict(zip(key0,lay))
        print("\n")
        print("****************************************")
        print("There are total %.d layers are found"  %k)
        la_n=int(input("Please input the number of layer to be selected :"))
        la_n1=la_n-1
        la_n0=la_n-2
        sel_low=lay_num[la_n0]
        sel_up=lay_num[la_n1]
        atoms=[]
        for i in range(sel_low,sel_up):
            atoms.append(layer[i][0])
        atom=open("atom1.dat",'w+')
        for i in atoms:
            atom.write(str(i)+' ')
        atom.write("\n")
        atom.close()
        print("**********************************")
        print(" The atom data has been written in 'atom1.dat' ")
        print(" The PDOS dat has been produced with vaspkit if it works ")
        print(" The log file by vaspkit has been write in 'log_vaspkit'" )
        os.system(" echo '11' >> select.dat ")
        os.system(" echo '114' >> select.dat")
        os.system(" cat atom.dat >> select.dat")
   
        os.system(" vaspkit < select.dat >log_vaspkit")
   else:
        print("You have input wrong,please rerun this script !!")


