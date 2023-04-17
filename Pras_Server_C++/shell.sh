#!/bin/bash

target="pdb/$1"

for f in "$target"*
do
    ./PRAS -i $f -f xxxx -r yes -m yes -l no -h yes -p no -s yes -o $(basename $f) -k yes
    #python3 PlotStructure.py jet 
    echo "fixed" $(basename $f)
    echo " "
done

