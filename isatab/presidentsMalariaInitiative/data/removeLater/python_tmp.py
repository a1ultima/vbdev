
import pdb

fi = open("./locationsOnly_pmi_dataset.txt")


locations = []

count = 0

while True:

    line = fi.readline()

    count += 1

    if line == "":
        break
    elif line[0] == "#":
        headers = line.split("\t")

        #                                                   Priority:
        # 0: '#Province             1st admin level',       6 [0] 1
        # 1: 'District 2nd admin level'                     5 [1] 3
        # 2: 'District Update'                              4 [2] 2
        # 3: 'Commune 3rd admin level'                      1 [3] 6  <-- 
        # 4: 'Village or Locality (site) ORIG'              3 [4] 5
        # 5: 'Village or Locality UPDATED (Co-ordinates)\n' 2 [5] 4

        continue

    split_line = line.split("\t")


    province            = split_line[0]          # 6
    district            = split_line[1]          # 5
    district_update     = split_line[2]          # 4
    commune             = split_line[3]          # 1  <-- 
    village             = split_line[4]          # 3
    village_update      = split_line[5].rstrip() # 2
    

    # VB IDEAL: leads to len(locations) == 338... too many
    # if not commune == "":
    #     locations.append(commune)             
    # elif not village_update == "":
    #     locations.append(village_update)
    # elif not village == "":
    #     locations.append(village)
    # elif not district_update == "":
    #     locations.append(district_update)
    # elif not district == "":
    #     locations.append(district)
    # elif not province == "":
    #     locations.append(province)

    if not province == "":
        locations.append(province)

    elif not district_update == "":
        locations.append(district_update)
    elif not district == "":
        locations.append(district)

    elif not commune == "":
        locations.append(commune)

    elif not village_update == "":
        locations.append(village_update)
    elif not village == "":
        locations.append(village)
    else:
        raise Exception("There are no localities for row: "+str(count))

fi.close()





