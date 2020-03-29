awk '{print $61}' results/patmatrix.txt |sort -rgk 1 | awk '{if($1<0.2)print NR}' |head -n 1
