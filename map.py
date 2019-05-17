#map peptides in file to database

import mysql.connector
import re
from read_file import read_file

mydb = mysql.connector.connect(host="localhost",user="root",passwd="mym19983280",database="protein_array_test1")
mycursor = mydb.cursor()

sql = "SELECT id, type, ESS, PSS, is_predict, proteins, proteins_id FROM peptide_new WHERE sequence ='{}'"

def exp_map(exp):
    eDir = {}
    for seq in exp:
        sql_sel = sql.format(seq)
        mycursor.execute(sql_sel)
        myresult = mycursor.fetchall()
        if len(myresult) == 0:
            continue
        pro_name = myresult[0][5].split('|')
        pro_id = myresult[0][6].split('|')
        id = myresult[0][0]
        type = myresult[0][1]
        ESS = float(myresult[0][2]) if myresult[0][2] != None else None
        PSS = float(myresult[0][3]) if myresult[0][3] != None else None
        is_predict = myresult[0][4]

        uni = 1 if len(pro_name) == 1 else 0
        for pro in pro_id:
            if pro not in eDir:
                eDir[pro] = {}
                eDir[pro][seq] = {}
                eDir[pro][seq]['intensity'] = exp[seq]
                eDir[pro][seq]['id'] = id
                eDir[pro][seq]['type'] = type
                eDir[pro][seq]['ESS'] = ESS
                eDir[pro][seq]['PSS'] = PSS
                eDir[pro][seq]['is_predict'] = is_predict
            else:
                eDir[pro][seq] = {}
                eDir[pro][seq]['intensity'] = exp[seq]
                eDir[pro][seq]['id'] = id
                eDir[pro][seq]['type'] = type
                eDir[pro][seq]['ESS'] = ESS
                eDir[pro][seq]['PSS'] = PSS
                eDir[pro][seq]['is_predict'] = is_predict
    return eDir

def map(fi_name):
    matrix = read_file(fi_name)
    print('start mapping...')
    print('...')
    map_matrix = {}
    for e in matrix:
        map_matrix[e] = exp_map(matrix[e])
    print('mapping finished.')
    return map_matrix
