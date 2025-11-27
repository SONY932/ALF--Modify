Variable=(Run Lx_R Ly_R  UH_R J_Red_R J_Rose_R J_Yellow_R J_Green_R  Beta_R Dtau_R NW_R NS_R NB_R CPU_R )
Name=(Run Lx_R Ly_R  UH_R J_Red_R J_Rose_R J_Yellow_R J_Green_R  Beta_R Dtau_R NW_R NS_R NB_R CPU_R )
#     0   1    2     3    4       5        6          7           8     9      10   11   12   13
{
read -a Variable
echo ${Variable[@]}
echo ${Variable[0]}
while [ ! ${Variable[0]} = "stop" ];   do
    if [ ${Variable[0]} = "Y" ]; then
	
	echo "Hi"
	export B_R_dir=`echo ${Variable[8]} | sed s/"\.0"//`
	export U_R_dir=`echo ${Variable[3]} | sed s/"\.0"//`
	export J1_R_dir=`echo ${Variable[4]} | sed s/"\.0"//`
	export J2_R_dir=`echo ${Variable[5]} | sed s/"\.0"//`
	export J3_R_dir=`echo ${Variable[6]} | sed s/"\.0"//`
	export J4_R_dir=`echo ${Variable[7]} | sed s/"\.0"//`
	
	export Dir="Lx"${Variable[1]}"_Ly"${Variable[2]}"_JRe"$J1_R_dir"_JRo"$J2_R_dir"_JYe"$J3_R_dir"_JGr"$J4_R_dir"_U"$U_R_dir"_B"$B_R_dir"_dtau"${Variable[9]}
	echo $Dir
	if [ ! -e $Dir ]; then
	    mkdir $Dir
	    cd  $Dir
	    cp  ../Start/* . 
	    cd ..
	else
	    cd $Dir
	    cp ../Start/parameters .
	    bash out_to_in.sh
	    cd ..
	fi
	cd $Dir
	
	#echo $Dir
	#echo  "nohup mpirun  -np 4 /home/debian/Kondo_KN_Code/Prog/KN_Kondo.out & " > run.sh  
	
	let i=1
	while [  $i -lt 14 ]; do
	    sed s/${Name[$i]}/${Variable[$i]}/    parameters  > tmp
	    mv tmp parameters
	    let i=i+1 
	done
	
	sed s/"DIR_R"/${Dir}/    job.sh  > tmp
        mv tmp job.sh
	sed s/"Name_R"/${Dir}/    job.sh  > tmp
        mv tmp job.sh

	sbatch job.sh
	
	cd ..
    fi
    read -a Variable
done
}<Sims

