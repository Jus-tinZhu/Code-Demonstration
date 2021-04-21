# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
from Bio.Blast import NCBIWWW
from Bio import SearchIO
import datetime
import pandas
from multiprocessing import Pool

sequences = ["YP_004452975", "YP_004452943", "YP_004452869", "YP_004452118", "YP_004453094", 
             "YP_004454342", "YP_004455163", "YP_004452446", "YP_004453545", "YP_004453539", 
             "YP_004453347", "YP_004452304", "YP_004452569", "YP_004454374", "YP_004454836", 
             "YP_004453828", "YP_004454314", "YP_004452080", "YP_004455183", "YP_004454121"]

seqdone = ["YP_004451722", "YP_004454730", "YP_004454773"]

organisms = ["Cellulomonas flavigena" , "Mycobacterium tuberculosis", 
        "Mycobacterium smegmatis", "Corynebacterium glutamicum", 
        "Corynebacterium diphtheriae", "Rhodococcus jostii", 
        "Rhodococcus erythropolis", "Streptomyces", 
        "Synechococcus sp", "Anabaena sp", "Bradyrhizobium japonicum"]
            
def Generate(index, sequence):
    seq = sequence
    organism = '"' + organisms[index] + '" [organism]'
    filename = "blastp " + organism + " "+ seq
    filename_xml = filename + ".xml"
    
    try:
        file = open(filename_xml, 'r')
        file.close()
        return None
    except IOError:
        return (seq, organism, filename)
    

def Query(seq, organism, filename):
    filename_xml = filename + ".xml"
    success = False
    
    while success == False:
        start = datetime.datetime.now()
        try:
            results = NCBIWWW.qblast("blastp", "nr", seq, expect = 40.0, entrez_query = organism)
            success = True
            print(filename + " successful")
        except Exception as e:
            print(e)
            print(filename + " unsuccessful")
        end = datetime.datetime.now()
        print(end-start)
        
    with open(filename_xml, "w") as out_handle:
        out_handle.write(results.read())
        results.close()


def Extract(filename):
    filename_xml = filename + ".xml"
    filename_pkl = filename + ".pkl"
    results = open(filename_xml, "r") 
    results_open = SearchIO.read(results, "blast-xml")
 
    value = []
    score = []
    names = []

    for i in range(len(results_open)):
        for j in range(len(results_open[i])):
            temp = results_open[i][0]
            
            value.append(temp.evalue)
            score.append(temp.bitscore)
            names.append(results_open[i].id)

    data = {"E-values": value, "Bitscores": score}
    dataframe = pandas.DataFrame (data, names)
    dataframe.to_pickle(filename_pkl)
    results.close()

    
def Worker(sequence):
    for i in range(len(organisms)):
        search = Generate(i, sequence)
        if search != None:
            print("Querying", search[0], search[1])
            Query(search[0], search[1], search[2])
            Extract(search[2])
            print(str(i+1) + " of " + str(len(organisms)) + " complete in " + str(sequence))
        else:
            print("skipping " + sequence + " " + organisms[i])
    
        
if __name__ == "__main__":
    with Pool(20) as pool:
        pool.map(Worker, sequences)
    print("complete")
