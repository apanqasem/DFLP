#!/bin/sh

limits="5000 10000 15000 20000 25000 30000 35000 40000 45000 50000 85000"

conf="../../data/d30l5/set1/input.txt"
seed=1023

# rm -rf omp_perf.txt omp_dyn_rc.txt

# make clean
# make dflp_dyn_rc


# for limit in $limits; do
#     echo -n $limit " " >> omp_dyn_rc.txt
#     get_primary.sh "./dflp $conf $seed $limit" >> omp_dyn_rc.txt
# done

# make clean
# make dflp
prog=dflp
echo "./$prog $conf $seed $limit > $prog.out" 
exit 0
for limit in $limits; do
    exec_time=`time (./$prog $conf $seed $limit > $prog.out) 2>&1` 
    exec_time=`echo "${exec_time}" | grep real | awk '{print $2}'\
                                   | sed 's/m/ /' | sed 's/s//' | awk '{printf "%3.2f", ($1 * 60 + $2)}'` 
    echo ${exec_time} >> omp_perf.txt
done

