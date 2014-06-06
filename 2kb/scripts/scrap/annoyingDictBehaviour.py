
dic  = {
    'please_make_me_Foo':{}, 
    'dont_make_him_Bar':{}
    }

appendic= {  # intended to be appended indepdendantly to each key of the above
    'Python_made_me':'',
}

dic['please_make_me_Foo']= appendic
dic['dont_make_him_Bar'] = appendic

dic['please_make_me_Foo']['Python_made_me'] = 'Foo' 

dic['dont_make_him_Bar']['Python_made_me']  = 'Bar' 

print(dic['please_make_me_Foo']['Python_made_me']) # Phew...
print(dic['please_make_me_Foo']['Python_made_me']) # What...?

