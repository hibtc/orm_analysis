#! /usr/bin/env zsh

tsp -S 4


data=../data/2018-10-20-orm_measurements/M8-E108-F1-I9-G1_
args=( -b $data/base/model_0 -m $data/ )

tsp python plot_folders.py $args $data/params.pos/*(/)
tsp python plot_folders.py $args $data/params.neg/*(/)
tsp python plot_folders.py $args $data/params.flip/*(/)

tsp python plot_folders.py $args $data/ealign.pos/*(/)
tsp python plot_folders.py $args $data/ealign.neg/*(/)
