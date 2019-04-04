#/bin/bash
read -p "
         Iput the file you want to copy :

          "     a   b   c  d  e  f g h
dir="/home/lxf/bin"
        for i in $a $b   $c  $d  $e  $f $g $h
        do
    chmod +x $i
        cp -f $i ${dir}
echo "  You have copied $i into ${dir} ! "
        done

    	