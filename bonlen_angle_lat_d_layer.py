import sys
import numpy as np
global name,lat,num_atom_types,atom_type_list,num_atoms,positions,layer0,layer1,bucklheight,bondlength1,bondlength2,atom0,atom1,atom2,name0,vector1,vector2,cos_angle,angle,sel,lattice,atom3,atom4,atom5,bondlength3,atom6
atom0=0;atom1=0;atom2=0;atom3=0;atom4=0;atom5=0;atom6=0;lattice=0;layer0=0;layer1=0 ;cos_agnle=0;angle=0;      
def input_vasp_file(file):
    f = open(file,'r')
    #read in structure name
    line = f.readline().split()
    name = ''
    for i in range(len(line)):
        name = name + line[i] + ' '
    #read in scale factor
    scale = float(f.readline().split()[0])
    #read in unit cell lattice vectors
    lat = np.zeros((3,3))
    for i in range(3):
        line = f.readline().split()
        for j in range(3):
            lat[i,j]=float(line[j])
    #read in number of atom types, number of atoms of each type,
    #and total number of atoms
    line=f.readline()
    line = f.readline().split()
    num_atom_types = len(line)
    atom_type_list = []
    for i in range(len(line)):
        atom_type_list.append(int(line[i]))
    num_atoms = 0
    for i in range(num_atom_types):
        num_atoms += atom_type_list[i]
    #read in atom coordinate convention
    convention = f.readline().split()[0]
    #read in atomic positions
    positions = np.zeros((num_atoms,3))
    for i in range(num_atoms):
        line = f.readline().split()
        for j in range(3):
            positions[i,j] = float(line[j])
    #convert atomic positions to cart coords in not already in them
    cees = ['c','C']
    if not cees.count(convention[0]):
        for i in range(num_atoms):
            positions[i] = np.dot(lat.transpose(),positions[i])
    #scale atomic positions and lattice vectors by scale factor
    lat = lat*scale
    f0=open(file,'r')
    f1=f0.readlines()
    f2=f1[0].split()[1:]
    firstline=[int(i) for i in f2]
    lata= np.sqrt((lat[0][0])**2+(lat[0][1])**2+(lat[0][2])**2)
    latb= np.sqrt((lat[1][0])**2+(lat[1][1])**2+(lat[1][2])**2)
    latc= np.sqrt((lat[2][0])**2+(lat[2][1])**2+(lat[2][2])**2)
    lata=lata/firstline[0]
    latb=latb/firstline[1]
    latc=latc/firstline[2]
    positions = positions*scale
    bucklheigth=0 #abobve need not change
    bucklheight=positions[layer1][2]-positions[layer0][2]  #need change accordingly
    if atom0 | atom1 |atom2 |atom3|atom4|atom5|atom6 > num_atoms:
       print(" you input a number large than a max_atom in POSCAR!")
       print("             Please rerun this code !")
       sys.exit(0)
    else:
       a1=(positions[atom0][0]-positions[atom1][0])**2
       b1=(positions[atom0][1]-positions[atom1][1])**2
       c1=(positions[atom0][2]-positions[atom1][2])**2
       a2=(positions[atom3][0]-positions[atom2][0])**2
       b2=(positions[atom3][1]-positions[atom2][1])**2
       c2=(positions[atom3][2]-positions[atom2][2])**2
       a3=(positions[atom3][0]-positions[atom4][0])**2
       b3=(positions[atom3][1]-positions[atom4][1])**2
       c3=(positions[atom3][2]-positions[atom4][2])**2
       bondlength1=np.sqrt(a1+b1+c1)
       bondlength2=np.sqrt(a2+b2+c2)
       bondlength3=np.sqrt(a3+b3+c3)
       vector1=positions[atom2]-positions[atom3]
       vector2=positions[atom4]-positions[atom3]
       if sel != 1:
          if sel !=4:
             if sel !=2:
                cos_angle=vector1.dot(vector2)/(bondlength2*bondlength3)
                angle=np.arccos(cos_angle)*360/2/np.pi
             else:
                angle=0
          else:
             angle=0
       else:
           angle=0
    if sel ==1:
       print("%4s, %9.3f" %(name0,bucklheight))
    elif sel==2:
         print("%4s,%9.3f" %(name0,bondlength1))
    elif sel==3:
         print("%4s,%9.3f" %(name0,angle))
    elif sel==4:
         print("%4s,%9.3f,%9.3f,%9.3f" %(name0,lata,latb,latc))
    elif sel==0:
         print("%4s,%12.3f,%12.3f,%12.3f,%8.3f,%8.3f,%8.3f" %(name0,bucklheight,bondlength1,angle,lata,latb,latc))
    #print("bondlength1: %6.3f" %(bondlength1))
    f.close()
print("\n")
print("#----------------------------split line------------------------------#")
print("| This code is used to calculate a bondlength,angle,lattice constant|")
print("|   and layer distance wit POSCAR                                    | ")
print("|   Written by XueFei Liu (Email:201307129@gznu.edu.cn)              |")
print("#-----------------------------split line-----------------------------#")

print("#-----------------------------split line-----------------------------#")
print("|1)For layer distance calculation,you need input the number of atom  |")
print("|  in different layers,e.g. atom N in layer 0, atom B in layer 1,    |")
print("|  you need input 0 1 respectively                                   |")
print("|2)For bondlenth,you need choose two atoms in POSCAR,the order of    |")
print("| atom number is sorted according with POSCAR,e.g.If you want to     |")
print("| calculate the bondlength between first and second atom,you need    |")
print("| input 0 and 1.etc.                                                 |") 
print("|3)For bondangle,you need choose three atoms in POSCAR,the order of  |")
print("|  atom number is sorted according with POSCAR,e.g.If you want to    |")
print("| calculate the bondanle between first  and second bondlength which  |")
print("|is constructed by atom 0 ,1 and 2,please input with the order 0,1,2 |")
print("|                           1                                        |" )
print("|                           /\                                       |")
print("|                          /  \                                      |")
print("|                         0    2                                     |")
print("|   Noting that for bondangles only superscell is surppoted !        |") 
print("|4) For lattice constant,Supercell of POSCAR must be produced by     |")
print("| vaspkit or set the first line as:e.g. system  4 4 1                |")
print("#-----------------------------split line-----------------------------#")

print("#-----------------------------split line-----------------------------#")
print("|Input 1 for layer distance , 2 for bondlength , 3 for bondangle and |")
print("|4 for lattice constant,if you want to calculate both please input 0 |")
print("#-----------------------------split line-----------------------------#")
sel= int(input("                           Input your chooice: "))
if sel == 1:
   layer0 = int(input("                    Input layer0 for layer distance: "))
   layer1 = int(input("                    Input layer1 for layer distance: "))
elif sel == 2:
     atom0=int(input("                     Input atom0 for bondlenth: "))
     atom1=int(input("                     Input atom1 for bondlenth: "))
elif sel == 3:
     atom2=int(input("                     Input atom0 for bondangle: "))
     atom3=int(input("                     Input atom1 for bondangle: "))
     atom4=int(input("                     Input atom2 for bondangle: "))
elif sel == 4:
     print( " Supercell of POSCAR must be produced by vaspkit or set the")
     print( " first line as:e.g. system  4 4 1   ")
elif sel==0:
     layer0 = int(input("           Input layer0 for layer distance: "))
     layer1 = int(input("           Input layer1 for layer distance: "))
     atom0=int(input("              Input atom0 for bondlenth: "))
     atom1=int(input("              Input atom1 for bondlenth: "))
     atom2=int(input("              Input atom0 for bondangle: "))
     atom3=int(input("              Input atom1 for bondangle: "))
     atom4=int(input("              Input atom2 for bondangle: "))


else:
    print(" You have input wrong ! Please rerun this code !" )
    sys.exit(0)	 
if sel==1:   
   print("        layer distance ")
elif sel==2:
     print("          bondlenth")
elif sel==3:
     print("          angle")
elif sel==4:
     print("             lata    ,latb    ,latc   ")
elif sel==0:
     print("          layer distance, bondlenth,   angle      ,lata   ,latb   ,latc   ")
el1=[''] # you need change here according to yourself's case,e.g.If your file name is "POSCAR",you need change\
el2=[''] # el1=[''],el2=['']

#el1=['N','P','As','Sb','Bi']
#el2=['B','Al','Ga','In','Tl']
name0=''
for i in range(len(el1)):
    for j in range(len(el2)):
        name0=el2[i]+el1[j]+'POSCAR'#to avoid making wrong,you'd better use a POSCAR
        #file1=sys.argv[name0]

        input_vasp_file(name0)








