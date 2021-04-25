__all__ = ["BioPythDescriptor",
        "AminoAcidComposition",
        "CTDD",
        "AutoCorrelations",
        "AminoAcidSingle",
        "AminoAcidDipeptide",
        "composition",
        "Transition",
        "Distribution",
        "Geary",
        "Moran",
        "NBmoran"]

from PyBioMed.PyProtein import AAComposition
from PyBioMed.PyProtein import CTD
from PyBioMed.PyProtein import Autocorrelation as AutoCorrelation
from PyBioMed.PyProtein import  AAIndex
from Bio.SeqUtils.ProtParam import ProteinAnalysis


class Descriptors:


    """
    
    Descriptor base class all class inherit from this class
    
    """

    def __init__(self,sequence):
        self.sequence = sequence

    def total_descriptor(self):
        pass
    




class AminoAcidComposition(Descriptors):
    """
    
    Amino acid composition descriptors base class.

    calculate total amino acid composition:
        +AmnioAcid Composition
        +Dipeptide compostion
        +Tripeptide compostion

    """

    def __init__(self,sequence):
        super().__init__(sequence)

        self.AA = None
        self.DP = None
        self.TP = None
        self.compostion = None

    
    def _AAC(self):
        if self.AA == None:
            comp = AAComposition.CalculateAAComposition(self.sequence)
            self.AA = comp
    
        
    def _DPC(self,):
        if self.DP == None:
            comp =AAComposition.CalculateDipeptideComposition(self.sequence)
            self.DP = comp


    def _ADPC(self,):
        if self.TP == None:
            comp= AAComposition.CalculateAADipeptideComposition(self.sequence)
            self.TP = comp

        
    def total_descriptor(self):
        self._AAC()
        self._DPC()
        self._ADPC()
        self.compostion = [self.AA,self.DP,self.TP]
        final_dict = {key:value for x in self.compostion for key,value in x.items()}
        return final_dict



class AminoAcidSingle(AminoAcidComposition):

    """

    class for extracting only single amino acid composition.

    """
    
    def __init__(self,sequence):
        super().__init__(sequence)


    def total_descriptor(self):
        self._AAC()
        return self.AA


class AminoAcidDipeptide(AminoAcidComposition):

    """
    class for extracting only single amino acid composition.
    
    """

    def __init__(self,sequence):
        super().__init__(sequence)

    
    def total_descriptor(self):
        self._DPC()
        return self.DP


class AminoAcidTripeptide(AminoAcidComposition):
    
    """

    class for extracting only single amino acid composition.
    
    """

    def __init__(self,sequence):
        super().__init__(sequence)

    
    def total_descriptor(self):
        self._ADPC()
        return self.TP





class CTDD(Descriptors):

    """
    Base class for CTDD descriptors
    composition, Decompostion, transition descriptors.
    
    
    """

    def __init__(self,sequence):
        super().__init__(sequence)
        self.cc = None
        self.tt = None
        self.dd = None
        self.total = None

        
    def _Compostion(self):
        if self.cc == None:
            comp = CTD. CalculateC(self.sequence)
            self.cc= comp


    def _Transition(self):
        if self.tt == None:
            comp = CTD. CalculateT(self.sequence)
            self.tt= comp


    def _Distribution(self):
        if self.dd == None:
            comp = CTD. CalculateT(self.sequence)
            self.dd= comp



    def total_descriptor(self):
        if self.total == None:
            comp = CTD.CalculateCTD(self.sequence)
            self.total = comp
        return self.total

class composition(CTDD):

    def __init__(self,sequence):
        super().__init__(sequence)

    def total_descriptor(self):
        self._Compostion()
        return self.cc

class Transition(CTDD):

    def __init__(self,sequence):
        super().__init__(sequence)

    def total_descriptor(self):
        self._Transition()
        return self.tt

class Distribution(CTDD):

    def __init__(self,sequence):
        super().__init__(sequence)

    def total_descriptor(self):
        self._Distribution()
        return self.dd



class AutoCorrelations(Descriptors):

    """

    NormalizedMB --> added
    MoranAC  --> added
    GearyAC --> added
    TotalAC --> added
    
    """

    def __init__(self,sequence):
        super().__init__(sequence)
        self.geary = None
        self.Moran = None
        self.NBmoran = None
        self.total = None
    
    def _NormalizedMB(self):
        if self.NBmoran == None:
            comb = AutoCorrelation.CalculateNormalizedMoreauBrotoAutoTotal(self.sequence)
            self.NBmoran = comb


    def _MoranAC(self):
        if self.Moran == None:
            comb = AutoCorrelation.CalculateMoranAutoTotal(self.sequence)
            self.Moran = comb
        

    def _GearyAC(self):
        if self.geary == None:
            comb = AutoCorrelation.CalculateGearyAutoTotal(self.sequence)
            self.geary = comb
    

    def _total_descriptor(self):
        if self.total == None:
            comb = AutoCorrelation.CalculateAutoTotal(self.sequence)
            self.total = comb

    
    def total_descriptor(self):
        self._total_descriptor()
        return self.total



class  Geary(AutoCorrelations):

    def __init__(self,sequence):
        super().__init__(sequence)

    def total_descriptor(self):
        self._GearyAC()
        return self.geary



class Moran(AutoCorrelations):

    def __init__(self,sequence):
        super().__init__(sequence)

    def total_descriptor(self):
        self._MoranAC()
        return self.Moran



class NBmoran(AutoCorrelations):

    def __init__(self,sequence):
        super().__init__(sequence)

    def total_descriptor(Self):
        self._NormalizedMB()
        return self.NBmoran



class Aaindex(Descriptors):

    def __init__(self,sequence):
        super().__init__(sequence)
        self.aa1 = None
        self.aa2 = None

    
    def _AAindex1(self):
        if self.aa1 == None:
            comb = AAIndex.GetAAIndex1(self.sequence)
            self.aa1 = comb


    def _AAindex2(self):
        if self.aa2 == None:
            comb = AAIndex.GetAAIndex23(self.sequence)
            self.aa2 = comb


    def total_descriptor(self):
        pass

    

class AAindex(Aaindex):
    def __init__(self,sequence):
        super().__init__(sequence)




class BioPythDescriptor(Descriptors):

    """

    Descriptors derived from the protein module of Biopython library.
    
    
    """

    def __init__(self,sequence):
        super().__init__(sequence)
        self._descrip = None


    def _Proteinanalysis(self):
        protein_molecule = ProteinAnalysis(self.sequence)
        AAC = protein_molecule.count_amino_acids()
        ARO = {"Aromaticity":protein_molecule.aromaticity()}
        gravy = {"gravy":protein_molecule.gravy()}
        instability = {"instability_index":protein_molecule.instability_index()}
        isoelectric_point = {"isoP":protein_molecule.isoelectric_point()}
        molecular_weight = {"MW":protein_molecule.molecular_weight()}
        secondary_structure = dict(zip(("Helix","Turn","Sheet"),protein_molecule.secondary_structure_fraction()))
        molar_extinction_coeff = dict(zip(("cysteine","cystines"),protein_molecule.molar_extinction_coefficient()))

        final_disc = {**AAC,**ARO,**gravy,**instability,**isoelectric_point,**molecular_weight,**secondary_structure,**molar_extinction_coeff}
        self._descrip = final_disc
        

    
    
    def total_descriptor(self):
        if self._descrip == None:
            self._Proteinanalysis()
        return self._descrip

