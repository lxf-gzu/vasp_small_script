#-*- coding: UTF-8 -*-
import re
import sys
import os
 
str1=[]
str2=[]
str_dump=[]
d_center=[]
fa=open("../../bulk/POSCAR",'r')
fb=open("POSCAR",'r')
fc=open("defect_center.dat",'w+')
 
for line in fa.readlines():
    str1.append(line.replace("\n",''))
for line in fb.readlines():
    str2.append(line.replace("\n",''))
 
for i in str1:
    if i in str2:
        str_dump.append(i)
 
str_all=set(str1+str2)
 
for i in str_dump:
    if i in str_all:
        str_all.remove(i)

for i in list(str_all):
    fc.write(i+'\n')
fa.close()
fb.close()
fc.close()
f_defect_r=open("defect_center.dat",'r')
f_defect_w=open("center_coord.dat",'w+')
#readline=f_defect_r.readlines()
for line in f_defect_r:
    d_center.append(line)
    #line=line.replace('\n','')
d_center.sort()
d_center1=d_center[0]
d_center2=d_center1[:-2]
print(d_center2)
d_center4=','.join(d_center2.split())


print("**********The coordinations of defect center*********** !")
f_defect_r=open("defect_center.dat",'r')
f_defect_w=open("center_coord.dat",'w+')
print (d_center4)
for i in list(d_center4):
    f_defect_w.write(i)
f_defect_w.write('\n')
f_defect_r.close()
f_defect_w.close()
