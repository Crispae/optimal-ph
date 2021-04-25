

__all__ = ["DescriptorsFrame",
            "PredictFrame",
            ]

import sys
from .utils import *
from .Pro_descriptor import *
from tqdm import tqdm
import pandas as pd

class DescriptorsFrame:

    """

        Create dataframe of descriptors
    
    """

    def __init__(self,sequence,Type="fasta",descriptors_list=[]):
        if Type == "fasta":
            if fasta_file_validator(sequence):
                self.file_path = sequence
                self.sequences,self.Target= sequence_extractor(self.file_path)

        else:
            print("error")
            self.sequences = []
            self.sequences.append(sequence.strip("\n"))
        self.descriptors_list = descriptors_list
        self.data = None
        self.rejected_data = []



    def _loader(self):
        """
        
        Load different descriptors from all different classes.
        
        
        """
        descriptors_name = ["BioPythDescriptor",
                    "AminoAcidComposition",
                    "CTDD",
                    "AutoCorrelations",
                    "AAC",
                    "ADP",
                    "CC",
                    "TT",
                    "DD",
                    "geary",
                    "moran",
                    "nbmoran"]

        descriptors =[BioPythDescriptor,
                    AminoAcidComposition,
                    CTDD,
                    AutoCorrelations,
                    AminoAcidSingle,
                    AminoAcidDipeptide,
                    composition,
                    Transition,
                    Distribution,
                    Geary,
                    Moran,
                    NBmoran
                    ]
        des_dict = dict(zip(descriptors_name,descriptors))

        if len(self.descriptors_list) == 0:
            descriptors_to_be_calc = descriptors
        else:
            descriptors_to_be_calc = [ des_dict[i] for i in descriptors_name for j in self.descriptors_list if i.strip() == j.strip() ]
            if len(descriptors_to_be_calc) == 0:
                raise ValueError("Enter valid Descriptor name")

        all_descriptors = []
        for index,seq in tqdm(enumerate(self.sequences)):
            try:
                for descrip in descriptors_to_be_calc:
                    all_descriptors.insert(index,descrip(seq).total_descriptor())
            except:
                self.rejected_data.append((index,seq))
                del self.Target[index]
                continue


        #all_descriptors = [descrip(seq).total_descriptor() for seq in tqdm(self.sequences) for descrip in descriptors_to_be_calc]
        self.data = all_descriptors
        Target = self.Target
        if len(Target) == len(self.data):
            return [self.data,Target]
    
    

    def show(self):
        if self.data == None:
            data,target = self._loader()
        dataframe = Data_frame(data,target)
        return dataframe


    

class PredictFrame:

    def __init__(self,sequence,Type="fasta",descriptors_list=[]):
        if Type == "fasta":
            if fasta_file_validator(sequence):
                self.file_path = sequence
                self.sequences =    predict_sequence(self.file_path)
        else:
            print("error")
            self.sequences = []
            self.sequences.append(sequence.strip("\n"))
        self.descriptors_list = descriptors_list
        self.data = None


    def _loader(self):
        """
        Load different descriptors from all different classes.
        
        """
        descriptors_name = ["BioPythDescriptor",
                    "AminoAcidComposition",
                    "CTDD",
                    "AutoCorrelations",
                    "AAC",
                    "ADP",
                    "CC",
                    "TT",
                    "DD",
                    "geary",
                    "moran",
                    "nbmoran"]

        descriptors = [BioPythDescriptor,
                    AminoAcidComposition,
                    CTDD,
                    AutoCorrelations,
                    AminoAcidSingle,
                    AminoAcidDipeptide,
                    composition,
                    Transition,
                    Distribution,
                    Geary,
                    Moran,
                    NBmoran
                    ]
        des_dict = dict(zip(descriptors_name,descriptors))

        if len(self.descriptors_list) == 0:
            descriptors_to_be_calc = descriptors
        else:
            descriptors_to_be_calc = [ des_dict[i] for i in descriptors_name for j in self.descriptors_list if i.strip() == j.strip() ]
            if len(descriptors_to_be_calc) == 0:
                raise ValueError("Enter valid Descriptor name")

        all_descriptors = []
        for index,seq in tqdm(enumerate(self.sequences)):
            try:
                for descrip in descriptors_to_be_calc:
                    all_descriptors.insert(index,descrip(seq).total_descriptor())
            except:
                self.rejected_data.append((index,seq))
                continue


        #all_descriptors = [descrip(seq).total_descriptor() for seq in tqdm(self.sequences) for descrip in descriptors_to_be_calc]
        self.data = all_descriptors
        return [self.data,]

    def raw(self):
        if self.data == None:
            self._loader()
        return self.data
        

    def show(self):
        if self.data == None:
            self._loader()
        dataframe = Predict_frame(self.data,)
        return dataframe
