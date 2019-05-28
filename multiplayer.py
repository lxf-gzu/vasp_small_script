#!/bin/python
# -*- coding:utf-8 -*-
"""
@author: XueFeiLiu(201307129@gznu.edu.cn),Nxu
@file: multiplayer.py
@time: 2019/1/29 21:44
A script to produce multi-layers of POSCAR
"""
import os
from numpy import *
import numpy as np
import shutil
import pprint
class VASP(object):
    def __init__(self,layers,layerdistance,mode,lata=2.178,latb=2.178,disp_n=200000000000,disp_flag=True,agl=30):
        self.lattice=np.zeros((3,3))
        self.Selective_infomation=False
        self.Direct_mode=True
        #self.mode=mode0
        self.mode0=mode
        self.layers=layers # layers
        self.layerdistance=layerdistance
        self.displace=disp_n
        self.flag_displace=disp_flag
        self.rotate0=mat([[cos((agl/360)*2*np.pi),sin((agl/360)*2*np.pi),0],[-sin((agl/360)*2*np.pi),cos((agl/360)*2*np.pi),0],[0,0,1]])
        self.M=self.rotate0
        self.rotate1=mat([[cos((-agl/360)*2*np.pi),sin((-agl/360)*2*np.pi),0],[-sin((-agl/360)*2*np.pi),cos((-agl/360)*2*np.pi),0],[0,0,1]])
        self.M1=self.rotate1
    def xyz_read(self):
        print('Now reading vasp structures.')
        if os.path.exists("POSCAR"):
            poscar=open("POSCAR",'r')
        elif os.path.exists("CONTCAR"):
            poscar=open("CONTCAR",'r')
        else:
            raise IOError('CONTCAR OR POSCAR does not exist！')

        self.title=poscar.readline().rstrip('\r\n').rstrip('\n')
        self.scaling_factor=float(poscar.readline())
        for i in range(3):
            self.lattice[i]=np.array([float(j) for j in poscar.readline().split()])
        #self.lattice*=self.scaling_factor
        self.element_list=[j for j in poscar.readline().split()]
        try:
            self.element_amount=[int(j) for j in poscar.readline().split()]
        except ValueError:
            raise ValueError('VASP 5.x POSCAR is needed!')
        line_tmp=poscar.readline()
        if line_tmp.strip().upper().startswith("S"):
            self.Selective_infomation=True 
            line_tmp=poscar.readline()  
        else:# no atoms fixed
            self.Selective_infomation=False 
        if line_tmp.strip().upper().startswith("D"):
            self.Direct_mode=True
        elif line_tmp.strip().upper().startswith("C"):
            self.Direct_mode=False
        else:
            raise ValueError("POSCAR format is not correct!")
        self.total_atom=sum(self.element_amount)
        self.atomic_position=np.zeros((self.total_atom,3))
        self.Selective_TF=[]
        if self.Selective_infomation == True:   
            for i in range(self.total_atom):
                line_tmp=poscar.readline()
                self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])
                self.Selective_TF.append([j for j in line_tmp.split()[3:]]) 
        else:
            for i in range(self.total_atom):
                line_tmp=poscar.readline()
                self.atomic_position[i]=np.array([float(j) for j in line_tmp.split()[0:3]])  
        self.atomic_position*=self.scaling_factor       
        if self.Direct_mode == True:
            self.cartesian_position=np.dot(self.atomic_position,self.lattice)
            self.direct_position=self.atomic_position
        else:
            self.cartesian_position=self.atomic_position
            self.direct_position=np.dot(self.atomic_position,np.linalg.inv(self.lattice))
        poscar.close()
        return (self.cartesian_position,self.direct_position)
    
    def xyz_write(self):
        print('Now writing new vasp structures.')
        if os.path.exists("POSCAR1"):
            shutil.copyfile("POSCAR1","POSCAR2")
            os.remove("POSCAR1")
        writen_lines=[]
        writen_lines.append(self.title)
        writen_lines.append(str(self.scaling_factor))
        for i in range(3):
            writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                .format(self.lattice[i][0],self.lattice[i][1],self.lattice[i][2]))
        writen_lines.append('  '+'  '.join(self.element_list))
        writen_lines.append('  '+'  '.join([str(j) for j in self.element_amount]))
        if self.Selective_infomation == True:
            writen_lines.append("Selective") 
        if self.mode0.upper().startswith('C'):            
            writen_lines.append("Cartesian")
            position_to_be_writen=self.cartesian_position 
        else:
            writen_lines.append("Direct") 
            position_to_be_writen=self.direct_position        
        position_to_be_writen/=self.scaling_factor
        if self.Selective_infomation == True:
            for i in range(np.size(position_to_be_writen,0)):
                writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                    .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                        position_to_be_writen[i][2],self.Selective_TF[i][0],\
                            self.Selective_TF[i][1],self.Selective_TF[i][2]))                    
        else:
            for i in range(np.size(position_to_be_writen,0)):
                writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                    .format(position_to_be_writen[i][0],position_to_be_writen[i][1], 
                        position_to_be_writen[i][2]))
        writen_lines=[j+'\n' for j in writen_lines]
        poscar=open("POSCAR1",'w')
        poscar.writelines(writen_lines)
        poscar.close()        
    
    def energy_extract(self):
        print('Now reading vasp energies.')
        if os.path.exists("OUTCAR"):
            outcar=open("OUTCAR",'r')
            for index,line in enumerate(outcar):
                if "energy(sigma->0)" in line:
                    E0=float(line.split()[-1])
            outcar.close()
            self.energy=E0
            return E0         
        elif os.path.exists("OSZICAR"):
            oszicar=open("OSZICAR",'r')
            for index,line in enumerate(oszicar):
                if "E0=" in line:
                    E0=float(line.split()[4])
            oszicar.close()
            self.energy=E0
            return E0                         
        else:
            raise IOError('OSZICAR OR OUTCAR does not exist！')

    
    def multilayer(self):
        print('Now produce new multi layers vasp structures.')
        if os.path.exists("POSCAR_mullayer"+str(self.layers)):
            shutil.copyfile("POSCAR_mullayer"+str(self.layers),"POSCAR_mullayer"+str(self.layers)+".bak")
            os.remove("POSCAR_mullayer"+str(self.layers))
        writen_lines=[]
        writen_all_layers=[]
        writen_each_layers=[]
        writen_lines.append(self.title)
        writen_lines.append(str(self.scaling_factor))
        self.lattice[2][2]+=self.layerdistance*(self.layers-1)/self.scaling_factor
        self.lattice_a=sqrt(self.lattice[0][0]**2+self.lattice[0][1]**2+self.lattice[0][2]**2)
        self.lattice_b=sqrt(self.lattice[1][0]**2+self.lattice[1][1]**2+self.lattice[1][2]**2)

        for i in range(3):
            writen_lines.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
            .format(self.lattice[i][0],self.lattice[i][1],self.lattice[i][2]))
        writen_lines.append('  '+'  '.join(self.element_list))
        self.element_amount0=self.element_amount[:]
        self.element_amount=[self.element_amount0[i]*self.layers for i in range(len(self.element_amount0))]
        writen_lines.append('  '+'  '.join([str(j) for j in self.element_amount]))
        if self.Selective_infomation == True:
            writen_lines.append("Selective") 
        if self.mode0.upper().startswith('C'):            
            writen_lines.append("Cartesian")
            position_to_be_writen=self.cartesian_position 
        else:
            writen_lines.append("Direct") 
            position_to_be_writen=self.direct_position        
        position_to_be_writen/=self.scaling_factor
        if self.Selective_infomation == True:
           if self.mode0.upper().startswith('D'):
              for i in range(np.size(position_to_be_writen,0)):
                  position_to_be_writen[i][2]-=(self.layers-1)*self.layerdistance/self.scaling_factor/2/self.lattice[2][2]
           else:
              pass
           for num_layer in range(self.layers):
               if num_layer < 1:
                  for i in range(np.size(position_to_be_writen,0)):
                      if self.mode0.upper().startswith('D'):
                         
                         writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                         .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                         position_to_be_writen[i][2],self.Selective_TF[i][0],\
                         self.Selective_TF[i][1],self.Selective_TF[i][2])) 
                        
                      else:
                         
                         writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                         .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                         position_to_be_writen[i][2],self.Selective_TF[i][0],\
                         self.Selective_TF[i][1],self.Selective_TF[i][2])) 
               else:         
                  if self.flag_displace == True:
                     for i in range(np.size(position_to_be_writen,0)):
                         if self.mode0.upper().startswith('D'):
                            position_xyz=mat(position_to_be_writen[i])*self.M
                            position_to_be_writen[i]=position_xyz

                            position_to_be_writen[i][0]+=self.lattice_a/self.displace/self.scaling_factor

                            position_to_be_writen[i][1]+=self.lattice_b/self.displace/self.scaling_factor

                            position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]
                            writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                            .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                            position_to_be_writen[i][2],self.Selective_TF[i][0],\
                            self.Selective_TF[i][1],self.Selective_TF[i][2])) 
                            position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]
                         else:
                            position_xyz=mat(position_to_be_writen[i])*self.M
                            position_to_be_writen[i]=position_xyz

                            position_to_be_writen[i][0]+=self.lattice_a/self.displace/self.scaling_factor
                            position_to_be_writen[i][1]+=self.lattice_b/self.displace/self.scaling_factor
                            position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor
                            writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                            .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                            position_to_be_writen[i][2],self.Selective_TF[i][0],\
                            self.Selective_TF[i][1],self.Selective_TF[i][2])) 
                            position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor
                     self.flag_displace = False
                  else:
                     for i in range(np.size(position_to_be_writen,0)):
                         if self.mode0.upper().startswith('D'):
                            position_xyz=mat(position_to_be_writen[i])*self.M1
                            position_to_be_writen[i]=position_xyz

                            position_to_be_writen[i][0]-=self.lattice_a/self.displace/self.scaling_factor

                            position_to_be_writen[i][1]-=self.lattice_b/self.displace/self.scaling_factor

                            position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]
                            writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                            .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                            position_to_be_writen[i][2],self.Selective_TF[i][0],\
                            self.Selective_TF[i][1],self.Selective_TF[i][2])) 
                            position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]
                         else:
                            position_xyz=mat(position_to_be_writen[i])*self.M1
                            position_to_be_writen[i]=position_xyz

                            position_to_be_writen[i][0]-=self.lattice_a/self.displace/self.scaling_factor
                            position_to_be_writen[i][1]-=self.lattice_b/self.displace/self.scaling_factor
                            position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor
                            writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}   {3:4}{4:4}{5:4}" \
                            .format(position_to_be_writen[i][0],position_to_be_writen[i][1], \
                            position_to_be_writen[i][2],self.Selective_TF[i][0],\
                            self.Selective_TF[i][1],self.Selective_TF[i][2])) 
                            position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor
                     self.flag_displace = True  
           writen_all_layers=[j+'\n' for j in writen_all_layers]
           writen_each_layers=[writen_all_layers[k:k+self.total_atom] for k in range(0,len(writen_all_layers),self.total_atom)]            
        else:
            if self.mode0.upper().startswith('D'):
               for i in range(np.size(position_to_be_writen,0)):
                   position_to_be_writen[i][2]-=(self.layers-1)*self.layerdistance/self.scaling_factor/2/self.lattice[2][2]
            else:
                pass
            for num_layer in range(self.layers):
                if num_layer < 1:
                   for i in range(np.size(position_to_be_writen,0)):
                           if self.mode0.upper().startswith('D'):

                              position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]

                              writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                              .format(position_to_be_writen[i][0],position_to_be_writen[i][1],position_to_be_writen[i][2]))
                              position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]

                           else:

                              position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor
                              writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                              .format(position_to_be_writen[i][0],position_to_be_writen[i][1],position_to_be_writen[i][2]))
                              position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor

                else:
                    
                    if self.flag_displace == True:
                       for i in range(np.size(position_to_be_writen,0)):
                           if self.mode0.upper().startswith('D'):
                              position_xyz=mat(position_to_be_writen[i])*self.M
                              position_to_be_writen[i]=position_xyz

                              position_to_be_writen[i][0]+=self.lattice_a/self.displace/self.scaling_factor

                     
                              position_to_be_writen[i][1]+=self.lattice_b/self.displace/self.scaling_factor


                              position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]

                              writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                              .format(position_to_be_writen[i][0],position_to_be_writen[i][1],position_to_be_writen[i][2]))
                              position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]

                           else:
                              position_xyz=mat(position_to_be_writen[i])*self.M
                              position_to_be_writen[i]=position_xyz
                          
                              position_to_be_writen[i][0]+=self.lattice_a/self.displace/self.scaling_factor

                              position_to_be_writen[i][1]+=self.lattice_b/self.displace/self.scaling_factor

                         
                              position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor
                              writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                              .format(position_to_be_writen[i][0],position_to_be_writen[i][1],position_to_be_writen[i][2]))
                              position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor
                       self.flag_displace = False
                    else:
 
                       for i in range(np.size(position_to_be_writen,0)):
                           if self.mode0.upper().startswith('D'):
                              position_xyz=mat(position_to_be_writen[i])*self.M1
                              position_to_be_writen[i]=position_xyz
                              position_to_be_writen[i][0]-=self.lattice_a/self.displace/self.scaling_factor


                              position_to_be_writen[i][1]-=self.lattice_b/self.displace/self.scaling_factor


                              position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]

                              writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                              .format(position_to_be_writen[i][0],position_to_be_writen[i][1],position_to_be_writen[i][2]))
                              position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor/self.lattice[2][2]

                           else:
                              position_xyz=mat(position_to_be_writen[i])*self.M1
                              position_to_be_writen[i]=position_xyz

                              position_to_be_writen[i][0]-=self.lattice_a/self.displace/self.scaling_factor
                              position_to_be_writen[i][1]-=self.lattice_b/self.displace/self.scaling_factor
  
                              position_to_be_writen[i][2]+=num_layer*self.layerdistance/self.scaling_factor
                              writen_all_layers.append("{0:>15.8f}{1:>15.8f}{2:>15.8f}" \
                              .format(position_to_be_writen[i][0],position_to_be_writen[i][1],position_to_be_writen[i][2]))



                              position_to_be_writen[i][2]-=num_layer*self.layerdistance/self.scaling_factor
                       self.flag_displace = True
                   
            writen_all_layers=[j+'\n' for j in writen_all_layers]
            #pprint.pprint(writen_all_layers)
            writen_each_layers=[writen_all_layers[k:k+self.total_atom] for k in range(0,len(writen_all_layers),self.total_atom)]             
            #pprint.pprint(writen_each_layers)
        writen_lines=[j+'\n' for j in writen_lines]
        poscar=open("POSCAR_mullayer"+str(self.layers),'w')
        poscar.writelines(writen_lines)
        j=0
        
        #print(self.element_amount0)
        for i in self.element_amount0:
                        
            for num_layer in range(self.layers):
                poscar.writelines(writen_each_layers[num_layer][j:j+i])
            j+=i    
       
        poscar.close()        

if __name__ == "__main__":
    a=str("Given the number of layers ,distance beween two layers and type of atoms coordination you want to build:\n")
    words=len(a)
    print("*"*words)
    print(a)
    print("*"*words)

    #n_layer,ldistance=input("Input the number of layers and distance:\n")
    n_layer=int(input("Input the number of layers:\n"))
    ldistance=float(input("Input the distance of two layers:\n"))
    coordtype=str(input("Input the type of atoms coordination(C or D):\n"))
    angle0=float(input("Input the angle you want to rotate:\n"))
    #distance_n=float(input("Input the displace number according to supercell:\n")) 
    poscar=VASP(layers=n_layer,layerdistance=ldistance,mode=coordtype,agl=angle0,disp_n=1000000000000000)
   
    cartesian_position,direct_position=poscar.xyz_read()
    #print(cartesian_position)
    poscar.multilayer()
    




