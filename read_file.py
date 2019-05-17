#read evidence.txt from maxquant
import re

def read_file(fi_name):
    # transfer file to intensity matrix
    print('start reading file...')
    print('...')
    fi_trans = open(fi_name, 'r',encoding='UTF-8')
    matrix = {}
    line1 = fi_trans.readline()
    lines1 = line1.strip('\n').split('\t')
    con = rev = e = 0
    for i, value in enumerate(lines1):
        if re.search(r'Experiment', value):
            e = i
        elif re.search(r'Reverse', value):
            rev = i
        elif re.search(r'contaminant', value):
            con = i
    for line in fi_trans:
        li = line.strip('\n').split('\t')
        if line.startswith('Sequence') or li[con] == '+' or li[rev] == '+':
            continue
        seq = li[0]
        exp = li[e]
        intensity = float(li[-11]) if li[-11] != '' else 0
        if exp not in matrix:
            matrix[exp] = {}
        matrix[exp][seq] = matrix[exp][seq] + intensity if seq in matrix[exp] else intensity
    fi_trans.close()
    print('file reading finished')
    return matrix







