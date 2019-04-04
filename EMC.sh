#!/bin/bash
path=`pwd`
echo "
      **************************************************************************
        The VPKIT.in for EMC calculation will be modified by you with the follow
      steps. The KPOINTS and VPKIT.in produced by vaspkit will be displayed for 
      you to setup the VPKIT.in in src_file folder.
      **************************************************************************
      
      "
 
cat ./KPOINTS
cat ./VPKIT.in
read -p  "
        *************************************************************************
        Please set the number for second line in VPKIT.in (for EMC calculation) 
        *************************************************************************
             Now ,please input the number: "  second_line
expr 1 + $second_line > /dev/null 2>&1
b=`echo $?`
if [ $b == 0 ]
then
  sed -i "2c  $second_line    !Typical value 4-8, number of points to fit second-order polynomial E=(h*k)^2/2m  " ./VPKIT.in
else
echo "
     ****************************************************
      you have input wrong, please rerun this program ! 
     ****************************************************
      "
 read -t 1     
exit

fi          

 read -p "
         *************************************************************************
          Please set the number for third line in VPKIT.in (for EMC calculation) 
         *************************************************************************
             Now ,please input the number: "  third_line 
b=$(echo " "$third_line" > 0.01" |bc )

if [ $b == 0 ]
then
   sed -i "3c $third_line    !Typical value 0.005-0.01 (in 1/A), kpoint spacing to determine the fitting precision of effective-mass " ./VPKIT.in
else

echo "
     ****************************************************
      you have input wrong, please rerun this program ! 
     ****************************************************
      "
 read -t 1     
exit

fi          
 
 
 read -p "
        *************************************************************************
         Please set the number for fourth line in VPKIT.in (for EMC calculation) 
        *************************************************************************
             Now ,please input the number: " fourth_line 
expr 1 + $fourth_line > /dev/null 2>&1
b=`echo $?`
if [ $b == 0 ]
then

   sed -i "4c $fourth_line   ! Number of task " ./VPKIT.in
 else
echo "
     ****************************************************
      you have input wrong, please rerun this program ! 
     ****************************************************
      "
 read -t 3     
exit

fi            
   
cat ./KPOINTS
cat ./VPKIT.in
   for (( k=1;k<=$fourth_line;k++))
  	do
  		q=`expr 4 + $k`
read -p  "
     *************************************************************************
     Please set the k-path of  "$q"th line in VPKIT.in (for EMC calculation) 
     according to your KPOINTS.
     *************************************************************************
             Now ,please input the values: "    k1 k2 k3 k4 k5 k6 k7		
    
   sed -i ""$q"c  $k1 $k2 $k3 $k4 $k5 $k6 $k7  !Calculate the effective-mass of K point along the "$k7" direction " ./VPKIT.in
done

cat ./VPKIT.in

 
sed -i "1c   1                                                     ! 1 for preprocess, 2 for postprocess " VPKIT.in

 vaspkit -task 913 -kpr 0.02 -kps G
sed -i "1c   2                                                     ! 1 for preprocess, 2 for postprocess " VPKIT.in
vaspkit -task 101 -inp ST
sed - i "s/EDIFF  =  1E-08/EDIFF  =  1E-05/g" INCAR

read -p  " **************************************************
           Input the number of cores you want to calculate:
	  **************************************************
		  "  ncore
mpirun -np $ncore vasp_std.base 
vaspkit -task 913 >sum.emc
sed -n '/^.*Summary/{:a;N;${s/\(.*KPATH No[^\n]*\).*/\1/p};Ta}' sum.emc>sum.emc1 #obtain the summary of EMC
mv sum.emc1 sum.emc
echo "
     ***********************************************************************************
      The  EMC have been calculated,please check the sum.emc.
     ***********************************************************************************
     "

