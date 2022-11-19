#!/bin/bash

target="cif/$1"

for f in "$target"*
do
    ./PRAS -i $f -f xxxx -r yes -m yes -l no -h yes -p no -s yes
    python3 PlotStructure.py jet 
    echo "fixed" $(basename $f)
    printf $(basename $f)"\n">> rec.txt
    echo " "
done

