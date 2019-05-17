#probe information

import map
import mysql.connector
import re
import time
from read_file import read_file

mydb = mysql.connector.connect(host="localhost",user="root",passwd="mym19983280",database="protein_array_test1")
mycursor = mydb.cursor()

def probe(matrix, fo_name):
    fo = open(fo_name, 'w')
    print('start writing probe information')
    print('...')
    pro_dir = {}
    exp_list = []
    line1 = 'probe_id\tprotein\t'
    for exp in matrix:
        if exp not in exp_list:
            exp_list.append(exp)
        for pro in matrix[exp]:
            if pro not in pro_dir:
                pro_dir[pro] = {}
                for seq in matrix[exp][pro]:
                    seqid = matrix[exp][pro][seq]['id']
                    seqI = matrix[exp][pro][seq]['intensity']
                    pro_dir[pro][seqid] = {}
                    pro_dir[pro][seqid][exp] = seqI
            else:
                for seq in matrix[exp][pro]:
                    seqid = matrix[exp][pro][seq]['id']
                    seqI = matrix[exp][pro][seq]['intensity']
                    if seqid not in pro_dir[pro]:
                        pro_dir[pro][seqid] = {}
                        pro_dir[pro][seqid][exp] = seqI
                    else:
                        pro_dir[pro][seqid][exp] = seqI
    for e in exp_list:
        line1 = line1 + e + '\t'
    fo.write(line1[:-1] + '\n')
    for p in pro_dir:
        for seq in pro_dir[p]:
            line = seq + '\t' + p + '\t'
            for exp in exp_list:
                if exp in pro_dir[p][seq]:
                    line = line + str(pro_dir[p][seq][exp]) + '\t'
                else:
                    line = line + 'NA' + '\t'
            line = line[:-1] + '\n'
            fo.write(line)
        fo.write('\n')
    fo.close()
    print('probe information writed')
    return pro_dir






