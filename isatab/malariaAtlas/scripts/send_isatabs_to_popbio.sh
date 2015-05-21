# copy isatabs: a_collection.txt s_samples.txt a_species.txt from ../data/isatabs/$CONTINENT/* over to popbio/

for CONTINENT in Asia Americas Africa Europe

do
	echo $CONTINENT
	cp ~/0VB/isatab/malariaAtlas/data/isatab/$CONTINENT/* ~/popbio/data_andy/isatabs/andy_0411_MAP/datesFixed_linkColumnAdded/$CONTINENT/
done	

