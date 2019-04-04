#!bin/bash
# prepare defet.sh defect_coord_center.py dielect.py first in your current directory such as "./BP"
bulk_dir="/home/lxf/test/defect/defect_cal/GaN.bak/bulk"
dielectric_dir="/home/lxf/test/defect/defect_cal/GaN.bak/dielectric/" #The above two lines need change appendly
path=`pwd`

#calculation_dir=" as* vac* "
calculation_dir="vac_1_Ga  vac_2_N"
#read -p "
#         Input the encut in INCAR: " encut
#read -p " Input the defect center coordination in direct type: " center
cp dielect.py $dielectric_dir
cd $dielectric_dir
python dielect.py #obtain dielect constant
eps=`cat die_harmonic_mean.dat`
cd $path
for i in $calculation_dir
    do
	 
	#read -p " Input the defect center coordination in direct type: " center
	cd $i
	
	for j in *
	    do
	    cd $j
		 ulimit -s unlimited

		ecut=`grep ENCUT INCAR`
		encut=${ecut#*=}
		encut=$(echo "scale=2;$encut /13.606" |bc);#Ry and eV conversion
	    cp $path/defect_coord_center.py ./
		echo "$i\ $j-data producing...."
	    python defect_coord_center.py #obtain the defect_center_coordination
		center=`cat center_coord.dat`
		charge=`echo ${j#*_}|awk '{print int($0)}'`
		charge=$[-$charge]
		echo $encut $charge $eps $center $bulk_dir
		sxdefectalign --ecut $encut --charge $charge --eps $eps --center $center --relative --vdef LOCPOT  --vref $bulk_dir/LOCPOT --vasp  >sxdefectalign_dat
		
		cd ../
	done
	cd ../
done