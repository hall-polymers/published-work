#!/bin/bash

ckpt_file=checkpoint_file.dat
# Check if the "checkpoint file" exists
if [ -e ${ckpt_file} ]; then
    startIdx=`cat ${ckpt_file}`
else
    startIdx=0
fi

TEMPLATES=../../2.templates

trajnumber=-11
for i in `seq ${startIdx} 9`; do
    let trajnumber+=1
    echo $trajnumber
    cp ${TEMPLATES}/dump2conf.py .
    sed -i "s/<traj>/$trajnumber/" dump2conf.py
    wait
    python3 dump2conf.py ${TEMPLATES}/in_100_6.lammps traj.lammpstrj in_100_6-"$i".lammps
    wait
    rm dump2conf.py
    wait
done

cp ${TEMPLATES}/orientation.py .
for j in `seq ${startIdx} 9`; do
    python3 orientation.py in_100_6-"$j".lammps
    wait
done

cp ${TEMPLATES}/std_dev.py .
wait
python3 std_dev.py status
wait
rm orientation.py in_* std_dev.py
