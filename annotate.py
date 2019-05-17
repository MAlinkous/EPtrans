#annotation file

#probe information

import map
import mysql.connector
import re
import time
from read_file import read_file
from probe import probe

mydb = mysql.connector.connect(host="localhost",user="root",passwd="mym19983280",database="protein_array_test1")
mycursor = mydb.cursor()

sql_pro = "select * from protein_new where id = '{}'"
sql_pep = "select * from peptide_new where id = '{}'"
def annotate(pro_dir, fo_name, title):
    print('start writing annotation file')
    print('...')

    T = time.strftime("%Y-%m-%d", time.localtime())

    fo = open(fo_name, 'w')
    fo.write('protein microarray v_1.0\n\n' + 'date = ' + T + '\n\n' + 'title = ' + title + ' annotation\n\n' + 'organism = Homo sapiens\n\n' + '...\n\n')
    fo.write('#ID\tSpecies\tUniprot Ids\tDate\tType\tDescription\tProtein(s) Title\tGene(s) Title\tGene(s) Ontology Biological Process\tGene(s) Ontology Cellular  Component\tGene(s) Ontology Molecular Function\tpeptide ids\tpeptide counts\tunique peptide counts\n\n')
    fo.write('\t\t!ID\tPeptideAtlas Id\tDate\tProteins\tProteins Id\tSequence\tType\tESS\tOBS\tEOS\tN expriments\tPSS\tIs Predict\n\n')
    for pro in pro_dir:
        sqlp1 = sql_pro.format(pro)
        mycursor.execute(sqlp1)
        myresult = mycursor.fetchall()[0]
        fo.write('#' + pro + '\t')
        for i in myresult[1:-1]:
            fo.write(str(i) + '\t')
        fo.write(str(myresult[-1]) + '\n\n')
        for pep in pro_dir[pro]:
            sqlp2 = sql_pep.format(pep)
            mycursor.execute(sqlp2)
            myresult2 = mycursor.fetchall()[0]
            fo.write('\t\t!' + pep + '\t')
            for j in myresult2[1:-1]:
                fo.write(str(j) + '\t')
            fo.write(str(myresult2[-1]) + '\n')
        fo.write('\n')
    fo.close()
    print('annotation file writed')










