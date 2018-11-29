#!/usr/bin/python3.6

from glob import glob
import os

files = glob('T*_*')

with open('parents.csv', 'wt') as outfile:
    outfile.write('Target,Domain,Model,Parent,Chain\n')
    for f in files:
        name = os.path.splitext(f)[0]
        target=name[0:5]
        model=name[5:12]
        domain=name.split('-')[1]
        lines = None
        with open(f, 'rt') as infile:
            lines = infile.read().split('\n')
        line = lines[0]
        i=0
        while 'PARENT' not in line:
            i+=1
            line = lines[i]
        if 'N/A' in line:
            continue
        line = line.split()
        for entry in line[3:]:
            parent = chain = ''
            if entry[0].isdigit() and len(entry) >=4:
                parent = entry[:4].lower()
                remainder = entry[4:]
                if "_" in remainder:
                    remainder = remainder.split('_')[1]
                if len(remainder) > 1:
                    remainder = remainder[0]
                    # raise RuntimeError("Remainder '{}' for entry {} with parent {} in directory {} has wrong length".format(
                    #     remainder, entry, parent, os.curdir
                    # ))
                if len(remainder):
                    chain = remainder[0]

                outfile.write('{},{},{},{},{}\n'.format(target, domain, model, parent, chain))
