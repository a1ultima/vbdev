echo "Error Messages: " > ./log.txt


for CONTINENT in Asia Americas Africa Europe

#echo python /home/ab108/software/oldiff.py --file1 /home/ab108/popbio/data_andy/isatabs/andy_0411_MAP/old/$CONTINENT/a_collection.txt --file2 /home/ab108/0VB/isatab/malariaAtlas/data/isatab/$CONTINENT/a_collection.txt > ./$CONTINENT.txt

do 

echo $CONTINENT

python /home/ab108/software/oldiff.py --file1 /home/ab108/popbio/data_andy/isatabs/andy_0411_MAP/old/$CONTINENT/a_collection.txt --file2 /home/ab108/0VB/isatab/malariaAtlas/data/isatab/$CONTINENT/a_collection.txt > ./$CONTINENT.txt


COUNT1=`grep \< $CONTINENT.txt | wc -l`
COUNT2=`grep \> $CONTINENT.txt | wc -l`
CHECK=`expr $COUNT1 - $COUNT2` 

if [ $CHECK -ne 0 ]; then
	echo "Suspicious differences: $CONTINENT";
	echo "Suspicious differences: $CONTINENT" >> ./log.txt;
	FAIL=1
fi 

done


