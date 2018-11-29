#!/bin/python3.6

if __name__ == '__main__':
    import os
    from glob import glob
    files = glob('*_?')
    for f in files:
        with open(f, 'rt') as infile:
            lines = infile.read().split('\n')
        i=0
        line = lines[0]
        try:
            with open(f+'.pdb', 'wt') as outfile:
                while 'MODEL' not in line:
                    outfile.write('REMARK   6 '+line+'\n')
                    i+=1
                    line=lines[i]
                # CASP PDB files have a line starting "PARENT " after the "MODEL "
                # line. Need to push this back to before the "MODEL " entry and
                # turn it into a REMARK
                model_line = line
                i+=1
                parent_line = lines[i]
                i+=1
                outfile.write('REMARK   6 '+parent_line+'\n')
                outfile.write(model_line+'\n')
                for line in lines[i:]:
                    outfile.write(line+'\n')
        except IndexError:
            print('Failed on {}'.format(f))
            os.remove(f+'.pdb')
