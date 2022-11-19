#!/bin/bash

target="pdb/$1"

for f in "$target"*
do
    ./PRAS -i $f -f xxxx -r yes -m yes -l no -h no -p no -s no
    #python3 PlotStructure.py jet 
    echo "fixed" $(basename $f)
    echo " "
done

