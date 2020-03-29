
ls steps_data/ > list.tetra
head -n 3 steps_data/AAAA/bp_step_all.par | tail -n 1 | awk '{for(i=1;i<=NF;i++)print $i}' > list.params
