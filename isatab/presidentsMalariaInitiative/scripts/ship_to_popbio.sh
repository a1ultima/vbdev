#
# Ships all the local isatab sheets to the popbio locale, ready for CHADO-loading....
#

# Check if target directory exists
if [ -d /home/ab108/popbio/data_andy/isatabs/andy_2015-07-17_PMI ]; then
    echo "..."
else
    echo "creating popbio target location..."
    mkdir /home/ab108/popbio/data_andy/isatabs/andy_2015-07-17_PMI
fi

# File conversions from ISO-8859-1 -to- UTF-8 (to allow them to be compatible w/ the postgres thing)

echo -e "\tConverting sheets that are encoded in 'ISO-8859-1' into 'UTF-8'..."

cd /home/ab108/0VB/isatab/presidentsMalariaInitiative/data/isatab/

echo -e "\ta_collection.txt"

iconv -f ISO-8859-1 -t UTF-8 a_collection.txt > a_collection_utf8.txt

mv a_collection_utf8.txt a_collection.txt

echo -e "\ta_IR_BA.txt"

iconv -f ISO-8859-1 -t UTF-8 a_IR_BA.txt > a_IR_BA_utf8.txt

mv a_IR_BA_utf8.txt a_IR_BA.txt

echo -e "\tp_IR_BA.txt"

iconv -f ISO-8859-1 -t UTF-8 p_IR_BA.txt > p_IR_BA_utf8.txt

mv p_IR_BA_utf8.txt p_IR_BA.txt

echo -e "\ta_IR_WHO.txt"

iconv -f ISO-8859-1 -t UTF-8 a_IR_WHO.txt > a_IR_WHO_utf8.txt

mv a_IR_WHO_utf8.txt a_IR_WHO.txt

echo -e "\tp_IR_WHO.txt"

iconv -f ISO-8859-1 -t UTF-8 p_IR_WHO.txt > p_IR_WHO_utf8.txt

mv p_IR_WHO_utf8.txt p_IR_WHO.txt

# Ship files over to popbio ready for loading

echo -e "\tShipping files over to popbio..."

cp /home/ab108/0VB/isatab/presidentsMalariaInitiative/data/isatab/* /home/ab108/popbio/data_andy/isatabs/andy_2015-07-17_PMI 
