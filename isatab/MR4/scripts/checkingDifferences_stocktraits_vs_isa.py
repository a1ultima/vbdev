raw=sorted(['MRA-121','MRA-913','MRA-114','MRA-860','MRA-763','MRA-765','MRA-186','MRA-112','MRA-115','MRA-762','MRA-105','MRA-861','MRA-111','MRA-334','MRA-698','MRA-594','MRA-339','MRA-764','MRA-856','MRA-1155','MRA-803','MRA-1156','MRA-493','MRA-130','MRA-891','MRA-139','MRA-1027','MRA-729','MRA-700','MRA-489','MRA-1154','MRA-128','MRA-126','MRA-726','MRA-735','MRA-730CDC','MRA-734','MRA-804','NR-43025','NR-43026','MRA-1163','MRA-1164','MRA-1165','MRA-1167','MRA-1166','MRA-978','MRA-728','MRA-893','MRA-314','MRA-862','MRA-863','MRA-864'])

isa=sorted(['MRA-105','MRA-111','MRA-112','MRA-114','MRA-115','MRA-121','MRA-126','MRA-128','MRA-130','MRA-139','MRA-186','MRA-314','MRA-334','MRA-339','MRA-489','MRA-493','MRA-594','MRA-698','MRA-700','MRA-726','MRA-728','MRA-729','MRA-730-CDC','MRA-734','MRA-735','MRA-762','MRA-763','MRA-764','MRA-765','MRA-803','MRA-804','MRA-854','MRA-856','MRA-860','MRA-861','MRA-862','MRA-863','MRA-864','MRA-891','MRA-893','MRA-913','MRA-978','MRA-1027','MRA-1154','MRA-1155','MRA-1156','MRA-1163','MRA-1164','MRA-1165','MRA-1166','MRA-1167'])


print 'isa is missing the following, i.e. raw has the following that isa does not have...'
print list( set(raw) - set(isa) )
# ['MRA-730CDC', 'NR-43025', 'NR-43026']

print 'raw is missing the following, isa has the following that raw does not have...'
print list( set(isa) - set(raw) )
# ['MRA-854', 'MRA-730-CDC']
