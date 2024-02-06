import os


acc_list = []
for file in os.listdir():
    if file.startswith('GCA'):
        if file.endswith('.fna'):
            working_string = file.split('_', 2)
            acc = working_string[0]+'_'+working_string[1]
            acc_list.append(acc)
            new_string = working_string[0]+'_'+working_string[1]+'.fna'

            os.rename(file, new_string)

with open('files.txt','w') as f:
    for item in acc_list:
        f.write(item)
