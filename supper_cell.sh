#author: XueFeiLiu(201307129@gznu.edu.cn,qq 316187631@qq.com)
echo "************By XueFei Liu****************"
read -p "Input the two m n groups produced by m_n_supercell.py: 
      " m1 n1 m2 n2 
(echo 4;echo 400;echo $m1 $n1 0;echo $m2 $n2 0;echo 0 0 1)|vaspkit	  
cp SUPERCELL.vasp POSCAR

echo " The POSCAR_REV has been produced by vaspkit !"
echo " Run multiplayer.py to obtain multiple-layer POSCAR if necessary !"
(echo 92;echo 921)|vaspkit 

cp POSCAR_REV POSCAR
(echo 92;echo 921)|vaspkit 
cp POSCAR_REV POSCAR