# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:16:57 2019

@author: XueFei Liu 316187631@qq.com
"""

import  numpy as np

def m_n(c,alpha_v,k1,k2,gamma0):
    k=k1**2+k2**2-k1*k2
    m_n1=[]
    for alpha in alpha_v:
       a1=[round(k1*np.cos(alpha),4),round(k1*np.sin(alpha),4)]
       b1=[round(-0.5*k2*np.cos(alpha)-np.sqrt(3)*k2/2*\
        np.sin(alpha),4),round(-0.5*k2*np.sin(alpha)\
        +np.sqrt(3)*k2/2*np.cos(alpha),4)]
       ab1=a1+b1
       a1b1=round(np.dot(a1,b1),6) #a1 b1 点乘
       a1_v=round(np.linalg.norm(a1),6)
       b1_v=round(np.linalg.norm(b1),6) #取 a1 b1 
       
       
       delta=0.0003
       
       for m in range(-c,c):
           for n in range(-c,c):
               mn=m**2+n**2-m*n
               ab2=[m*a[0]+n*b[0],m*a[1]+n*b[1]] #由ma+nb 构成
               ab2_v=round(np.linalg.norm(ab2),4) #ab2 模
               ab2b1=round(np.dot(ab2,b1),4) #ab2 与b1 点乘
               c_1=round(ab2b1/(b1_v*ab2_v),4)
               np.seterr(invalid='ignore')
               
               if abs(ab2[0]-(ab1[0])) < 0.003\
                 and abs(ab2[1]-(ab1[1]))\
                 <0.003 and abs(c_1 + np.cos(gamma0))\
                 <abs(np.cos(gamma0))*2+0.00003 and abs(mn-k) <delta:
               
                 m_n1.append([m,n])           
    m_n2=[]
    for i in m_n1:
        if i not in m_n2:
           m_n2.append(i)
    
    return m_n2
if __name__=='__main__':
    a=[1,0]
    b=[-1/2,np.sqrt(3)/2]  # a b crystal vector\
    # of hexagonal system of primitive cell
    alpha_i=-360 #angle rot range
    evalue=4 # a value almost < sqrt(k1)
    k1=7
    k2=7 #sqrt(k1)^2
    gamma0=120 # ab 夹角
    alpha2=np.linspace(0,(alpha_i/360)*2*np.pi,1000)
    
    print("************By XueFei Liu****************")
    k1=int(input("Input sqrt(k1):  "))
    k2=int(input("Input sqrt(k1):  "))
    m_n=m_n(evalue,alpha2,round(np.sqrt(k1),4),\
            round(np.sqrt(k1),4),(gamma0/360)*2*np.pi)
    print("m n groups as below:")
    if len(m_n) == 0:
        print("No m n group found,please check it ")
    else:
        
        print(m_n)