#!/bin/sh

prog=dflp
limit=4000
conf="../../data/d30l5/set1/input.txt"
seed=1023

exec_time=`time (./$prog $conf $seed $limit > $prog.out) 2>&1` 
exec_time=`echo "${exec_time}" | grep real | awk '{print $2}'\
                                   | sed 's/m/ /' | sed 's/s//' | awk '{printf "%3.2f", ($1 * 60 + $2)}'` 
echo ${exec_time} >> perf.txt