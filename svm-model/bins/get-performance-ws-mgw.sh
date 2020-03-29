
i=1
rm -f names.tmp
for i in 1 2 3 4 5 
do
echo "mgw-bin.$i" >> names.tmp
done

awk '{if(NR>1)print $(3),$NF}' results/performance-loo-ws1-mgw.txt > tmp1.txt
awk '{if(NR>1)print $NF}' results/performance-loo-ws3-mgw.txt > tmp2.txt
awk '{if(NR>1)print $NF}' results/performance-loo-ws5-mgw.txt > tmp3.txt
awk '{if(NR>1)print $NF}' results/performance-loo-ws7-mgw.txt > tmp4.txt

echo "SD MAE(WS1) MAE(WS3) MAE(WS5) MAE(WS7)" > results/performance-ws-mgw.txt
paste names.tmp tmp1.txt tmp2.txt tmp3.txt tmp4.txt  >> results/performance-ws-mgw.txt

rm -f tmp[1-3].txt names.tmp

