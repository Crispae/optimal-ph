from Pro_descriptor import BioPythDescriptor
import numpy as np 
from Bio import SeqIO
from sklearn.preprocessing import OneHotEncoder 
from sklearn.ensemble import RandomForestRegressor

def main(): 

    file_dir = "../data/28protien.fasta"
    seq_np = fasta_to_np(file_dir)

    ph_dir = "../data/28pdbid_Ph.txt"
    ph_values = np.array(read_ph(ph_dir)).ravel()
    print(ph_values)
    

    one_hot = OneHotEncoder()
    seq_list = one_hot.fit_transform(seq_np)
    print(seq_list)

    rf_reg = RandomForestRegressor()
    rf_reg.fit(seq_list, ph_values)
    
    test = rf_reg.predict(seq_list)

    print(test)

def fasta_to_np(file_name): 
    seqs = SeqIO.parse(file_name, "fasta")
    seq_list = [str(x.seq) for x in seqs]
    return np.array(seq_list).reshape(-1, 1)

def read_ph(file_name): 
    with open(file_name, "r") as file: 
        data = file.read().split()

    print(data)
    return data[1::2]

if __name__ == "__main__":
    main()