#Please cite: Tran-Nguyen, V. K., Junaid, M., Simeon, S. & Ballester, P. J. A practical guide to machine-learning scoring for structure-based virtual screening. Nat. Protoc. (2023)

#Provide the file name of the hit list:
echo "Please provide the file name of your hit hist, make sure that the hit list has been sorted from the best (highest) to the worst (lowest):"
read hitlistname

#Count the number of test molecules:
i=$(wc -l $hitlistname | awk '{print $1}')-1
N=$((i))

#Calculate the number of test molecules in the top 1%:
n=$(($N/100))

#Count the number of true active molecules (true hits) in the top 1%:
m=$(($n+1))
a=$(sed -n 2,${m}p $hitlistname | grep -c 'Active')

#Count the number of true active molecules (true hits) in the whole test set:
A=$(grep -c 'Active' $hitlistname)

#Compute the EF1% value:
EF1=$(echo "scale=2; 100*$a/$A" | bc)
echo "EF1% = $EF1"

#Compute the maximal EF1% value:
if [ "$A" -le "$n" ]; then
	EF1max=100
	echo "EF1%max = $EF1max"
elif [ "$A" -gt "$n" ]; then
	EF1max=$(echo "scale=2; 100*$n/$A" | bc)
	echo "EF1%max = $EF1max"
fi

#Compute the NEF1% value:
NEF1=$(echo "scale=3; $EF1/$EF1max" | bc)
echo "NEF1% = $NEF1"

