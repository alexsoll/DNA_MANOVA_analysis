import numpy as np
import pandas as pd
import sys
import os
from statsmodels.multivariate.manova import MANOVA
from numpy.testing import assert_raises
import random

def get_dict():
    dict = {}
    file = open("annotations.txt", 'r')
    tmp = file.readline()
    while True:
        tmp = file.readline()
        if (len(tmp) == 0):
            break
        tmp = tmp.split("\t")
        dict[tmp[0]] = tmp[1:]
    print(len(dict))
    print(sys.getsizeof(dict) / (1024*1024))
    file.close()

def correct_data():  

    #        This function returns a dictionary with the original data, 
    #         among which data relating to the sex chromosomes 
    #          and the geo type Island and Shore were excluded.

    file = open("annotations.txt", 'r')
    dict = {}
    tmp = file.readline()
    while True:
        tmp = file.readline()
        if (len(tmp) == 0):
            break
        #if (tmp.find("Island") == -1 and tmp.find("Shore") == -1):
        if(tmp.find("Shelf") == -1 ):
            tmp = tmp.split("\t")
            if (tmp[1] != "X" and tmp[1] != "Y" and tmp[1] != '' and tmp[9] != '' and tmp[5] != '' and tmp[10] != 'classC' and tmp[10] != 'classD'):
                dict[tmp[0]] = tmp[1:]
            else:
                continue
        else:
            continue
    return dict

def get_dict_cpg_gene(dict):    #Creating a dictionary, where the key - cpg, data - genes
    dict_cpg_genes = {}
    for key in dict:
        dict_cpg_genes[key] = dict[key][4]
    return dict_cpg_genes

def dict_cpg_g_without_repetition(dict):
    dict_cpg_genes = {}
    for key in dict:
        genes = dict[key][4].split(";")
        tmp = ''
        for i in range(len(genes)):
            if genes[i] in tmp:
                continue
            else:
                tmp += genes[i] + ';'
        dict_cpg_genes[key] = tmp
    return dict_cpg_genes

def get_dict_gene_cpg(dict):    #Creating a dictionary, where the key - gene, data - cpg
    #sorted(dict.items(), key = lambda item: item[1][4])
    dict_genes_cpg = {}
    for key in dict:
        tmp = dict[key][4].split(";")
        for j in range(len(tmp)):
            if tmp[j] not in dict_genes_cpg:
                dict_genes_cpg[tmp[j]] = key
            else:
                dict_genes_cpg[tmp[j]] += ";" + key
    return dict_genes_cpg

################## Dictionaty bop-cpg & cpg-bop ###################

def get_dict_cpg_bop(dict):
    file = open("CPG_BOP.txt", "w")
    d = {}
    for key in dict:
        file.write(str(key) + " " + str(dict[key][10]) + "\n")
        d[key] = dict[key][10]
    file.close()
    return d

def get_dict_bop_cpg(dict_cpg_bop):
    d = {}
    file = open("BOP_CPG.txt", "w")
    for key in dict_cpg_bop:
        bop = dict_cpg_bop[key]
        if dict_cpg_bop[key] in d:
            d[bop] += ";" + key
        else:
            d[bop] = key
    for key in d:
        file.write(str(key) + " " + str(d[key]) + "\n")

    file.close()
    return d
################################################################

def result_table():
    print("Wait, the operation is in progress ...")
    d = correct_data()
    #cpg_gene = get_dict_cpg_gene(d)  # Dictionary with repetitions genes
    cpg_gene = dict_cpg_g_without_repetition(d) # Dictionary without repetitions genes
    result_tab = {}
    file = open("average_beta.txt", 'r')
    file.readline()
    for line in file:
        line = line.split()
        if line[0] in cpg_gene:
            genes = cpg_gene.get(line[0])
            genes = genes[:len(genes) - 1] # Used for dictionary without gene repetition
            genes = genes.split(";")
            for i in range(len(genes)):
                if genes[i] not in result_tab:
                    result_tab[genes[i]] = line[1:]
                else:
                    for k in range(1, len(line)):
                        result_tab[genes[i]][k - 1] += ' ' + line[k]
        else:
            continue
    for key in result_tab:
        num_of_cpg = len(result_tab[key][0].split(" "))
        for i in range(len(result_tab[key])):
            sum = 0
            for j in range(num_of_cpg):
                tmp = result_tab[key][i].split(" ")
                sum += float(tmp[j])
            result_tab[key][i] = float(sum / num_of_cpg)
        

    file.close()
    return result_tab

def print_result_dict(dict):
    file = open("result_new.txt", "w")
    for key in dict:
        file.write(str(key) + '\t')
        for i in range(len(dict[key])):
            file.write(str(dict[key][i]) + " ")
        file.write("\n")
    clear = lambda: os.system('cls')
    clear()
    print("Done!  The data is saved in the file result.txt.")
    file.close()


def get_ages():
    file_age = open("attributesGSE87571.txt", "r")
    ages = []
    file_age.readline()
    i = 0
    while(True):
        line = file_age.readline()
        if (len(line) == 0):
            break
        line = line.split()
        age = int(line[3])
        ages.append(age)
        i += 1
    return ages

def MANOVA_analysis(dict_cpg_bop, dict_bop_cpg):
    dict_BoP_PValue = {}
    age = get_ages()
    file = open("average_beta.txt", "r")
    file.readline()
    for line in file:
        line = line.split()
        name_cpg = line.pop(0)
        if name_cpg in dict_cpg_bop:
            bop = dict_cpg_bop[name_cpg]
            l = dict_bop_cpg[bop].split(";")
            if len(l) < 3:
                continue
            else:
                if bop in dict_BoP_PValue:
                    dict_BoP_PValue[bop].append(line)
                else:
                    dict_BoP_PValue[bop] = []
                    dict_BoP_PValue[bop].append(line)
    file = open("DataFrame.txt", "w")
    print(len(dict_BoP_PValue))
    num = 0
    for key in dict_BoP_PValue:
        print(num)
        num += 1
        dict = {}
        pVal = []
        l = len(dict_BoP_PValue[key])
        for i in range(0, l - 2):
            cpg1 = []
            cpg2 = []
            cpg3 = []

            cpg1 = list(np.float_(dict_BoP_PValue[key][i]))
            cpg2 = list(np.float_(dict_BoP_PValue[key][i+1]))
            cpg3 = list(np.float_(dict_BoP_PValue[key][i+2]))

            #for j in range(len(dict_BoP_PValue[key][i])):
            #    cpg1.append(float(dict_BoP_PValue[key][i][j]))
            #    cpg2.append(float(dict_BoP_PValue[key][i+1][j]))
            #    cpg3.append(float(dict_BoP_PValue[key][i+2][j]))

            DatFrame = pd.DataFrame({'age': age,
                              'cpg1' : cpg1,
                              'cpg2' : cpg2,
                              'cpg3' : cpg3
                              })
            #DatFrame.to_csv(file, header=None, index = None, sep=' ', mode='a') 
            #DatFrame.to_csv(file, sep=' ', mode='a') 
            model = MANOVA.from_formula('cpg1 + cpg2 + cpg3 ~ age', data=DatFrame)
            test = model.mv_test()
            pVal.append(test.results['age']['stat'].values[3,4])
        pVal.sort()
        min_pVal = pVal[0]
        dict_BoP_PValue[key] = min_pVal
    '''
    age = get_ages()
    for i in range(1):
        dict = {}
        cpg1 = []
        cpg2 = []
        cpg3 = []
        for i in range(728):
            tmp = random.random()
            cpg2.append(tmp)
            cpg3.append(tmp)
            cpg1.append(tmp)
        tmp = 0.954697456795
        cpg2.append(tmp)
        cpg3.append(tmp+0.000001)
        cpg1.append(tmp+0.00001)


        #DatFrame = pd.DataFrame({'age': age,
        #                         'cpg1': cpg1,
        #                         'cpg2': cpg2,
        #                         'cpg3': cpg3
        #                        })
        dict['age'] = age;
        dict['cpg1'] = cpg1
        dict['cpg2'] = cpg2
        dict['cpg3'] = cpg3
        #print(DatFrame)
        model = MANOVA.from_formula('cpg1 + cpg2 + cpg3 ~ age', data=dict)
        test = model.mv_test()
        res = test.results['age']['stat'].values[1,4]
        print(res)
    '''
    return dict_BoP_PValue



def main():
    dict = correct_data()
    cpg_bop = get_dict_cpg_bop(dict)
    bop_cpg = get_dict_bop_cpg(cpg_bop)
    print(len(bop_cpg))
    res = MANOVA_analysis(cpg_bop, bop_cpg)
    
    print(res)
    print(len(res))

    file = open("result.txt", 'w')
    for key in res:
        file.write(str(key) + " " + str(res[key]) + "\n")
    file.close()

main()
