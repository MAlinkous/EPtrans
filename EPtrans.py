#main program

import map
import mysql.connector
import re
import time
from read_file import read_file
from probe import probe
from annotate import annotate
from ESSTopInten import TopInten
import sys
import getopt
import os

print(time.asctime(time.localtime(time.time())))
mydb = mysql.connector.connect(host="localhost",user="root",passwd="mym19983280",database="protein_array_test1")
mycursor = mydb.cursor()


def usage():
    print("==========")
    print("Usage:")
    print("-i Input file")
    print("-p output FilePath")
    print("-t output file label")
    print("-h Help Message")
    print("-m quantificaton method")
    print("-n top n probes")
    print("-r Sorting criteria")
    print("==========")

def getargv():
    argx = {}
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hp:t:m:n:r:i:")
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python EPtrans.py -i <inputfile> -p <outputfilepath> -t <outputfilelabel> -m <quantificaton method> -n <top n probes> -r <Sorting criteria>')
            usage()
            sys.exit()
        elif opt in '-i':
            argx['fi'] = arg
        elif opt in '-p':
            argx['fopath'] = arg + '/'
        elif opt in '-t':
            argx['label'] = arg
        elif opt in '-m':
            argx['method'] = arg
        elif opt in '-n':
            argx['n'] = int(arg)
        elif opt in '-r':
            argx['rank'] = arg
    if 'fopath' not in argx:
        path = './output/'
        if not os.path.exists(path):
            os.mkdir(path)
        argx['fopath'] = path
    if 'fopath' in argx:
        if not os.path.exists(argx['fopath']):
            print('path not exit!')
            sys.exit(2)

    return argx

def writematrix(matrix,argsx, name):
    pro_dir = {}
    for exp in matrix:
        for pro in matrix[exp]:
            if pro not in pro_dir:
                sql = "SELECT uniprot_ids FROM protein_new WHERE id ='{}'"
                sql = sql.format(pro)
                mycursor.execute(sql)
                myresult = mycursor.fetchall()
                pro_names = re.sub('_', ',', myresult[0][0])
                pro_dir[pro] = pro_names

    fo = open(argsx['fopath'] + argsx['label'] + '_' + name + '_protein_matrix.txt', 'w') if 'label' in argsx else open(
        argsx['fopath']  + name +'_protein_matrix.txt', 'w')
    line1 = 'id' + '\t' + 'proteins' + '\t'
    for exp in matrix:
        line1 = line1 + exp + '\t'
    fo.write(line1.strip('\t') + '\n')
    for pro in pro_dir:
        fo.write(pro + '\t' + pro_dir[pro] + '\t')
        line = ''
        for exp in matrix:
            if pro in matrix[exp]:
                line = line + str(matrix[exp][pro][0]) + '\t'
            else:
                line = line + 'NA' + '\t'
        fo.write(line.strip('\t') + '\n')
    fo.close()
    print('matrix file writed.')

def main():
    args = getargv()

    #quantify output
    matrix_map = map.map(args['fi'])

    matrix_inten_Top = TopInten(matrix_map,args['method'],args['n'],args['rank'])

    tag = args['method'] + '_' + str(args['n']) + '_' + args['rank']

    writematrix(matrix_inten_Top, args, tag)

    #probe information

    probe_fo = args['fopath'] + args['label'] + '_probe_information.txt' if  'label' in args else args['fopath'] + 'probe_information.txt'
    annot_dir = probe(matrix_map, probe_fo)

    annot_file = args['fopath'] + args['label'] + '_annotation.txt' if  'label' in args else args['fopath'] + 'annotation.txt'
    label = args['label'] if 'label' in args else ''
    annotate(annot_dir, annot_file, label)

    print('Done!')




if __name__ == '__main__':
    main()
print(time.asctime(time.localtime(time.time())))








