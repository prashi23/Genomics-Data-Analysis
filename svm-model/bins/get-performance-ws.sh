

awk '{if(NR>1)for(i=1;i<=5;i++)print $1"-bin-"i}' ../params_stats/results/breaks5bins.txt > names.tmp

awk '{if(NR>1)print $(3),$NF}' results/performance-loo-ws1.txt > tmp1.txt
awk '{if(NR>1)print $NF}' results/performance-loo-ws3.txt > tmp2.txt
awk '{if(NR>1)print $NF}' results/performance-loo-ws5.txt > tmp3.txt
awk '{if(NR>1)print $NF}' results/performance-loo-ws7.txt > tmp4.txt

echo "SD MAE(WS1) MAE(WS3) MAE(WS5) MAE(WS7)" > results/performance-ws.txt
paste names.tmp tmp1.txt tmp2.txt tmp3.txt tmp4.txt  >> results/performance-ws.txt

rm -f tmp[1-3].txt names.tmp

