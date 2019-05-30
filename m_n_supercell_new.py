#!/opt/miniconda3/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:16:57 2019

@author: XueFei Liu 316187631@qq.com
"""
import  numpy as np
def m_n(c,alpha_v,k1,k2,gamma0,theta):
    k=k1**2+k2**2-k1*k2
    m_n1=[]
    angle=[]
   
    
    for alpha in alpha_v:
       a1=[round(k1*np.cos(alpha),4),round(k1*np.sin(alpha),4)]
       b1=[round(-np.sin(theta)*k2*np.cos(alpha)-np.cos(theta)*k2*\
        np.sin(alpha),4),round(-np.sin(theta)*k2*np.sin(alpha)\
        +np.cos(theta)*k2*np.cos(alpha),4)]
       ab1=a1+b1
       a1b1=round(np.dot(a1,b1),6) #a1 b1 点乘
       a1_v=round(np.linalg.norm(a1),6)
       b1_v=round(np.linalg.norm(b1),6) #取 a1 b1 
       a_v=round(np.linalg.norm(a),6)
       aa1=round(np.dot(a,a1),6)              
       delta=0.0003       
       for m in range(-c,c):
           for n in range(-c,c):
               mn=m**2+n**2-m*n
               ab2=[m*a[0]+n*b[0],m*a[1]+n*b[1]] #由ma+nb 构成
               ab2_v=round(np.linalg.norm(ab2),4) #ab2 模
               ab2b1=round(np.dot(ab2,b1),4) #ab2 与b1 点乘
               c_1=round(ab2b1/(b1_v*ab2_v),4)
               c_aa1=round(aa1/(a_v*b1_v),4)
               aa1_angle=round(np.arccos(c_aa1)*360/(2*np.pi),2)
               
               np.seterr(invalid='ignore')               
               if abs(ab2[0]-(ab1[0])) < 0.003\
                 and abs(ab2[1]-(ab1[1])) <0.003 \
                 and abs(mn-k) < delta and abs(c_1 + np.cos(gamma0))\
                <abs(np.cos(gamma0))*2+0.00003:
               
                   m_n1.append([m,n])
                   #m_n1.append([n,m])
                   angle.append(aa1_angle)
                   
    #print(angle)            
    m_n2=[]
    for i in m_n1:
        if i not in m_n2:
           m_n2.append(i)    
    return m_n2,angle
def find_k(theta,m,n,a):    
    ab=[[a,0],[-a*np.sin(theta),a*np.cos(theta)]]
    mn=[m,n]
    mn_ab=np.mat(mn)*np.mat(ab)
    lat_lenth=round(np.linalg.norm(mn_ab),4)
    return lat_lenth    
def k_pick(a2=3.2740,b2=2.556,mismatch=0.05,n=5):
    k0=[]
    ktmp0=[]
    k1=[]
    k0_pick=[]
    k1_pick=[]
    ktmp1=[]
    theta0=(30/360)*np.pi*2    
    for i in range(1,n):  
        for j in range(1,n):            
            ktmp0.append(find_k(theta0,i,j,a2))
            ktmp1.append(find_k(theta0,i,j,b2))
    for i in ktmp0:
        if i not in k0:
            k0.append(i)
            ktmp0=[]
    for i in ktmp1:
        if i not in k1:
            k1.append(i)
            ktmp0=[]    
    for i in k0:
        for j in k1:
            if abs((i-j)/j) < mismatch:
               k0_pick.append(round((i/a2)**2,4))
               k1_pick.append(round((j/b2)**2,4))
    print("***The possible new constants***")
    print("----------------------------------")
    print("{:5s} {}".format("Sqrt(lat_a):",k0_pick))
    print("{:5s} {}".format("Sqrt(lat_b):",k1_pick)) 
    print("----------------------------------")          
    return k0_pick,k1_pick

if __name__=='__main__':
    mismatch=0.05                 
    gamma0=120 # ab 夹角
    a=[1,0]
    b=[-1/2,np.sqrt(3)/2]  # a b crystal vector\
    # of hexagonal system of primitive cell
    alpha_i=-360 #angle rot range
    evalue=5 # a value almost < sqrt(k1) 
    #if you are always can't find m n ,
    #you may increas '5' to a bigger value
    alpha2=np.linspace(0,(alpha_i/360)*2*np.pi,1000)   
    print("************By XueFei Liu****************")
    a2=float(input("The lattice constant for structure a:"))
    b2=float(input("The lattice constant for structure b:"))
    ka,kb=k_pick(a2=a2,b2=b2,mismatch=mismatch,n=evalue)
    print("*****************************************")
    print("Produce m n groups for structure a...")
    for i in range(len(ka)):
        k=round(ka[i],4)
        m_na,anglea=m_n(evalue,alpha2,round(np.sqrt(k),4),\
            round(np.sqrt(k),4),(gamma0/360)*2*np.pi,\
            theta=(30/360)*np.pi*2)        
        if len(m_na) == 0:
           print("{:4s}{} {} {}".format("m_n:"\
                 ,i+1,m_na,"No m n group found!"))
           print("{:4s}{} {} {}".format("twist_angle:"\
                 ,i+1,anglea,"No angle group found!"))
        else:
           print("{:4s}{} {}".format("m_n:"\
                 ,i+1,m_na))   
           print("{:4s}{} {}".format("twist_angle:"\
                 ,i+1,anglea))
    print("*****************************************")      
    print("Produce m n groups for structure b...")
    for i in range(len(kb)):
        k=kb[i]
        m_nb,angleb=m_n(evalue,alpha2,round(np.sqrt(k),4),\
            round(np.sqrt(k),4),(gamma0/360)*2*np.pi,\
            theta=(30/360)*np.pi*2)        
        if len(m_nb) == 0:
            print("{:4s}{} {} {}".format("m_n:"\
                 ,i+1,m_nb,"No m n group found!"))
            print("{:4s}{} {} {}".format("twist_angle:"\
                 ,i+1,angleb,"No angle group found!"))
        else:
        
           print("{:4s}{} {}".format("m_n:"\
                 ,i+1,m_nb))
           print("{:4s}{} {}".format("twist_angle:"\
                 ,i+1,angleb))
