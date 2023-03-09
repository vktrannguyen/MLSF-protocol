#Please cite: Tran-Nguyen, V. K., Junaid, M., Simeon, S. & Ballester, P. J. A practical guide to machine-learning scoring for structure-based virtual screening. Nat. Protoc. (2023)

#Provide the file name of the hit list:
echo "Please provide the file name of your hit hist, make sure that the hit list has been sorted from the best (highest) to the worst (lowest):"
read hitlistname

#Provide the name of the csv data file:
echo "Please provide the name of your csv data file:"
read datafilename

#Type a name with which the file containing the top 1% population is called:
echo "Please type a name with which you want to call the file containing the top 1% population (csv format):"
read toplistname

#Type a name with which the file containing the true hits in the top 1% is called:
echo "Please type a name with which you want to call the file containing the true hits in the top 1% population (csv format):"
read truehitsname

#Count the number of test molecules:
i=$(wc -l $hitlistname | awk '{print $1}')-1
N=$((i))

#Calculate the number of test molecules in the top 1%:
n=$(($N/100))

#Extract the top 1%-ranked population from the hit list:
m=$(($n+1))
sed -n 2,${m}p $hitlistname >& $toplistname

#Process the list of top 1%:
sed -i 's/,/   /g' $toplistname

#Extract the list of active molecules (true hits) in the top 1%:
grep 'Active' $toplistname | awk '{print $1}' >& $truehitsname

#Process the original data file (csv):
cp $datafilename data.csv
sed -i 's/,/   /g' data.csv
awk '{print $1}' data.csv >& col1
awk '{print $3}' data.csv >& col3
paste col1 col3 >& data.csv
rm col1 col3
