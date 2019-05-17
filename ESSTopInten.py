#quantify using top 3 ESS or intensity peptides median or mean

import map
import mysql.connector
import re
import time
import numpy as np
from read_file import read_file


def TopInten(matrix, method = 'mean', pnum = 10, rank = 'intensity'):
    print('start writing matrix file')
    print('...')
    matrix_inten = {}
    if method == 'median':
        me = np.median
    elif method == 'mean':
        me = np.mean
    else:
        me = sum

    for exp in matrix:
        matrix_inten[exp] = {}
        for pro in matrix[exp]:
            matrix_inten[exp][pro] = []
            if pro.endswith('au'):
                if len(matrix[exp][pro]) <= pnum:
                    intenli = [matrix[exp][pro][seq]['intensity'] for seq in matrix[exp][pro]]
                    inten = me(intenli)
                    matrix_inten[exp][pro].append(inten)
                else:
                    seq_all = [seq for seq in matrix[exp][pro]]
                    seq_all_rank = {}
                    for seq in seq_all:
                        seq_all_rank[seq] = matrix[exp][pro][seq][rank]
                    seq_all_rank_sorted = sorted(seq_all_rank.items(), key=lambda item: item[1], reverse=True)
                    seq_list = [s[0] for s in seq_all_rank_sorted[:pnum]]
                    inten = me([matrix[exp][pro][seq]['intensity'] for seq in seq_list])
                    matrix_inten[exp][pro].append(inten)

            if pro.endswith('an'):
                uFlag = False
                for seq in matrix[exp][pro]:
                    if matrix[exp][pro][seq]['type'] == 'uo':
                        uFlag = True
                        break
                if uFlag == True:
                    uniSeq = {}
                    for seq in matrix[exp][pro]:
                        if matrix[exp][pro][seq]['type'] == 'uo':
                            uniSeq[seq] = matrix[exp][pro][seq]
                    if len(uniSeq) <= pnum:
                        inten = me([uniSeq[seq]['intensity'] for seq in uniSeq])
                        matrix_inten[exp][pro].append(inten)
                    else:
                        seq_all = [seq for seq in uniSeq]
                        seq_all_rank = {}
                        for seq in seq_all:
                            seq_all_rank[seq] = uniSeq[seq][rank]
                        seq_all_rank_sorted = sorted(seq_all_rank.items(), key=lambda item: item[1], reverse=True)
                        seq_list = [s[0] for s in seq_all_rank_sorted[:pnum]]
                        inten = me([uniSeq[seq]['intensity'] for seq in seq_list])
                        matrix_inten[exp][pro].append(inten)
                else:
                    if len(matrix[exp][pro]) < 2:
                        delete = matrix_inten[exp].pop(pro)
                        continue
                    else:
                        if len(matrix[exp][pro]) <= pnum:
                            inten = me([matrix[exp][pro][seq]['intensity'] for seq in matrix[exp][pro]])

                            matrix_inten[exp][pro].append(inten)

                        else:
                            seq_all = [seq for seq in matrix[exp][pro]]
                            seq_all_rank = {}
                            for seq in seq_all:
                                seq_all_rank[seq] = matrix[exp][pro][seq][rank]
                            seq_all_rank_sorted = sorted(seq_all_rank.items(), key=lambda item: item[1], reverse=True)
                            seq_list = [s[0] for s in seq_all_rank_sorted[:pnum]]
                            inten = me([matrix[exp][pro][seq]['intensity'] for seq in seq_list])

                            matrix_inten[exp][pro].append(inten)
            if pro.endswith('m'):
                if len(matrix[exp][pro]) < 2:
                    delete = matrix_inten[exp].pop(pro)
                    continue

                numPep = pnum * 2
                uniSeq = {}
                nSeq = {}
                for seq in matrix[exp][pro]:
                    if matrix[exp][pro][seq]['type'] == 'um':
                        uniSeq[seq] = matrix[exp][pro][seq]
                nNum = numPep - len(uniSeq)
                for seq in matrix[exp][pro]:
                    if matrix[exp][pro][seq]['type'] == 'n':
                        nSeq[seq] = matrix[exp][pro][seq]
                if len(nSeq) <= nNum:
                    inten = me( [uniSeq[seq]['intensity'] for seq in uniSeq]  +  [nSeq[seq]['intensity'] for seq in nSeq])

                    matrix_inten[exp][pro].append(inten)

                else:

                    seq_all_rank = {}
                    for seq in nSeq:
                        seq_all_rank[seq] = nSeq[seq][rank]
                    seq_all_rank_sorted = sorted(seq_all_rank.items(), key=lambda item: item[1], reverse=True)
                    seq_list = [s[0] for s in seq_all_rank_sorted[:nNum]]
                    inten = me([uniSeq[seq]['intensity'] for seq in uniSeq] + [nSeq[seq]['intensity'] for seq in nSeq])
                    matrix_inten[exp][pro].append(inten)
            if pro.endswith('pr'):
                uoNum = uPnum = nNum = nPnum = 0
                for seq in matrix[exp][pro]:
                    type = matrix[exp][pro][seq]['type']
                    if type == 'uo':
                        uoNum += 1
                    elif type == 'uP':
                        uPnum += 1
                    elif type == 'n':
                        nNum += 1
                    else:
                        nPnum += 1
                if (uoNum+uPnum) >= 1:
                    seq_all = {}
                    for seq in matrix[exp][pro]:
                        if matrix[exp][pro][seq]['type'] == 'uo' or matrix[exp][pro][seq]['type'] == 'uP':
                            seq_all[seq] = matrix[exp][pro][seq]
                    inten = me([seq_all[seq]['intensity'] for seq in seq_all])

                    matrix_inten[exp][pro].append(inten)

                else:
                    if nNum >= pnum:
                        seq_all_rank = {}
                        for seq in matrix[exp][pro]:
                            if matrix[exp][pro][seq]['type'] == 'n':
                                seq_all_rank[seq] = matrix[exp][pro][seq][rank]
                        seq_all_rank_sorted = sorted(seq_all_rank.items(), key=lambda item: item[1], reverse=True)
                        seq_list = [s[0] for s in seq_all_rank_sorted[:pnum]]
                        inten = me([matrix[exp][pro][seq]['intensity'] for seq in seq_list])
                        matrix_inten[exp][pro].append(inten)

                    else:
                        if (nNum+nPnum) >= pnum:
                            seq_n = {}
                            for seq in matrix[exp][pro]:
                                if matrix[exp][pro][seq]['type'] == 'n':
                                    seq_n[seq] = matrix[exp][pro][seq]
                            seq_nP = {}
                            for seq in matrix[exp][pro]:
                                if matrix[exp][pro][seq]['type'] == 'nP':
                                    seq_nP[seq] = matrix[exp][pro][seq]
                            seq_nP_rank = {}
                            rankP = 'PSS' if rank == 'ESS' else 'intensity'
                            for seq in seq_nP:
                                seq_nP_rank[seq] = seq_nP[seq][rankP]
                            seq_nP_rank_sorted = sorted(seq_nP_rank.items(), key=lambda item: item[1], reverse=True)
                            seq_list = [s[0] for s in seq_nP_rank_sorted[:(pnum-nNum)]]
                            inten = me([seq_nP[seq]['intensity'] for seq in seq_list] + [seq_n[seq]['intensity'] for seq in seq_n])

                            matrix_inten[exp][pro].append(inten)

                        else:
                            seq_all = {}
                            for seq in matrix[exp][pro]:
                                if matrix[exp][pro][seq]['type'] == 'n' or matrix[exp][pro][seq]['type'] == 'nP':
                                    seq_all[seq] = matrix[exp][pro][seq]
                            inten = me([seq_all[seq]['intensity'] for seq in seq_all])

                            matrix_inten[exp][pro].append(inten)
    return matrix_inten




