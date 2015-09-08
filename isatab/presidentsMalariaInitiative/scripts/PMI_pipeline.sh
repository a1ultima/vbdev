# README v 2015/07/24
#
# Runs the ISA-Tab generator then ships to PopBio
#
#
cd /home/ab108/0VB/isatab/presidentsMalariaInitiative/scripts

echo "generating foreign keys (mapping RAW PMI values -to- VB freindly ontology terms)..."

python ./ontologyKeys_to_dict.py

echo "generating ISA-Tabs at: ../data/isatabs/*, but does not include i_investigations.txt..."

python ./dict_to_isatabGen.py           # data gets generated and put here: ~/0VB/isatab/presidentsMalariaInitiative/data/isatab/

echo "shipping ISA-Tabs to: /popbio/data_andy/isatabs/andy_2015-07-17_PMI/..."

bash ./ship_to_popbio.sh                # data then gets shipped here: /popbio/data_andy/isatabs/andy_2015-07-17_PMI/* where it will be loaded into Chado

# then follow the admin notes to load into chado: https://docs.google.com/document/d/1w_3nhphdkEJCzu97wqLor57k9G8m2SGGu3schWSHy4A/edit#

echo -e "\nSUCCESS!!!\n"

echo "data should be ready to load into Chado.."

echo "...you can find the data that needs loading into Chado from:     ~/popbio/data_andy/isatabs/andy_2015-07-17_PMI/"

echo "...and you can find the scripts to do this loading from:   ~/popbio/VBPopBio/api/Bio-Chado-VBPopBio/bin/"

echo -e "\n"

echo "e.g. cd ~/popbio/; cd VBPopBio/api/Bio-Chado-VBPopBio; bin/load_project.pl --dry-run ~/popbio/data_andy/isatabs/andy_2015-07-17_PMI/"

echo "here are clear instructions: https://docs.google.com/document/d/1w_3nhphdkEJCzu97wqLor57k9G8m2SGGu3schWSHy4A/edit"

# then run ... 