

__all__ = ["fasta_file_validator",
            "sequence_extractor",
            "Data_frame",
            "predict_sequence",
            "Predict_frame",]

import os
import pandas as pd

def fasta_file_validator(file):
    if os.path.isfile(file):
        return True
    return TypeError("File doesn't exist.")

def sequenceRepair(seq_list):
    repaired = []
    for index,seq in enumerate(seq_list):
        if ("X" in seq):
            seq = seq.replace("X","")
        if ("B" in seq):
            seq = seq.replace("B","")
        if ("Z" in seq):
            seq = seq.replace("Z","")
        
        repaired.insert(index,seq)
    return repaired
    

def predict_sequence(file):
    try:
        data = pd.read_csv(file)
    except:
        raise FileNotFoundError("File handling error occured")
    sequence = sequenceRepair(data["sequence"].to_list())
    
    return sequence


def sequence_extractor(file):
    try:
        data = pd.read_csv(file)
    except:
        raise FileNotFoundError("File handling error occured")
    sequence = sequenceRepair(data["sequence"].to_list())
    Target = data["mean_growth_PH"]

    return [sequence,Target.to_list()]

def Predict_frame(data,):
    descriptor_data = pd.DataFrame(data)
    return descriptor_data
    

def Data_frame(data,target):
    descriptor_data = pd.DataFrame(data)
    if descriptor_data.shape[0] == len(target):
        descriptor_data.insert(descriptor_data.shape[1],column="Target",value=target)
    return descriptor_data
    