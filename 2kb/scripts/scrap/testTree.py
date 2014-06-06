
import speciesManage

#moo = ('PHUMU','RPROC','AAEGL','CPIPJ','AALBS','ADARC','AARAD','AQUAS','AMELC','AMERM','AGAMP','ACHRA','AEPIE','ACULA','AMINM','AFUNF','AMACM','ASTEI','ASTES','ADIRW','AFARF','AATRE','ASINS','LLONJ','PPAPI','DANA1','DERE1','DYAK1','DMEL5','DSEC1','DSIM1','DPER1','DPSE3','DWIL1','DGRI1','DMOJ1','DVIR1','GMORY','BMORI','DPLEX','TCAST','AMELL','LHUMI')

species_list = speciesManage.generate_list('./species_list.txt') # generates a list of species names corresponding to EnsEMBl MySQL db names
species_to_motifs_to_pwms = {}


print [i.split('_')[0][0].upper()+i.split('_')[1][0:3].upper() for i in species_list]

#'anopheles_gambiae'.split('')