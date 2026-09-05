from ..lib import re
from ..utils import Seq

class Scenario:

    ''' NAGNAG splice scenarios '''

    def __init__(self, phase:int=None, up:str=None, motif:str=None, down:str=None, allele:str="ref"):
        
        ''' NAGNAG splice scenario object

        Parameters
        ----------
        phase : int {-1, 0, 1, 2}
            position of the motif in the reading frame
        up : str
            last 3 bases of the upstream exon
        motif : str
            6 splice site bases
        down : str
            first 3 bases following the motif
        '''

        self.reset()
        self.set_scenario(phase, up, motif, down, allele)

    def reset(self):

        # Default Attributes
        self.allele = "ref"
        self.empty = True
        
        # scenario information
        self.phase = None
        self.up_seq = None
        self.motif = None
        self.down_seq = None
        self.canon = False
        
        # flanking amino acids
        self.aa_up = None
        self.aa_down = None

        # variable codons
        self.vc_ps = None
        self.vc_ds = None

        # amino acids
        self.aa_ps = None
        self.aa_ds = None

        # transition type
        self.aatt = None

        # insertions, deletions
        self.ins_aa = None
        self.del_aa = None

    def set_scenario(self, phase:int=None, up:str=None, motif:str=None, down:str=None, allele:str="ref"):

        ''' Set splice scenario information

        Parameters
        ----------
        phase : int {-1, 0, 1, 2}
            position of the motif in the reading frame
        up : str
            last 3 bases of the upstream exon
        motif : str
            6 splice site bases
        down : str
            first 3 bases following the motif
        '''

        self.allele = allele if allele in ["ref", "alt"] else "ref"

        # Verify arguments
        phase = phase if (phase != None) else self.phase
        if phase not in [-1, 0, 1, 2]:
            return
        
        up = up if (up != None) else self.up_seq
        if not Seq.valid_seq(up) or (len(up) < 2):
            return
        
        motif = motif if (motif != None) else self.motif
        if not Seq.valid_seq(motif) or (len(motif) != 6):
            return
        
        down = down if (down != None) else self.down_seq
        if not Seq.valid_seq(down) or (len(down) < 2):
            return

        self.empty = False

        # Scenario Information
        self.phase = phase
        self.up_seq = up
        self.motif = motif
        self.down_seq = down
        self.canon = bool(re.fullmatch(".AG.AG", self.motif))
        
        # Flanking Amino Acids
        self.aa_up, self.aa_down = Scenario.flank_aa(self.phase, self.up_seq, self.down_seq)

        # Variable Codons
        self.vc_ps, self.vc_ds = Scenario.variable_codons(phase, up, motif, down)

        # Amino Acids
        self.aa_ps, self.aa_ds = Scenario.variable_amino_acids(self.vc_ps, self.vc_ds)

        # Transition Type
        self.aatt = Scenario.transition_type(self.aa_ps, self.aa_ds)

        # Inserted Amino Acid(s)
        self.ins_aa = Scenario.inserted_aa(self.aa_ps, self.aa_ds)
        
        # Deleted Amino Acid
        self.del_aa = Scenario.deleted_aa(self.aa_ps, self.aa_ds)

    @staticmethod
    def flank_aa(phase:int, up_seq:str, down_seq:str) -> tuple[str, str]:

        ''' Get the flanking amino acids

        Parameters
        ----------
        phase : int {-1, 0, 1, 2}
            position of the motif in the reading frame
        up_seq : str
            last 9+ bases of the upstream exon
        down_seq : str
            first 9+ bases following the motif

        Returns
        -------
        up_aa : str
            amino acid upstream (N-terminus) of variable amino acids
        ds : str
            amino acid downstream (C-terminus) of variable amino acids
        '''

        # Non-Coding
        if phase == -1:
            return "_", "_"

        # Coding
        if phase == 0:
            up = up_seq[-3:]
            down = down_seq[:3]
        elif phase == 1:
            up = up_seq[-4:-1]
            down = down_seq[2:5]
        elif phase == 2:
            up = up_seq[-5:-2]
            down = down_seq[1:4]
        
        # translate
        up_aa = Seq.translate(up) if up else None
        down_aa = Seq.translate(down) if down else None
        
        return up_aa, down_aa

    @staticmethod
    def variable_codons(phase:int, up:str, motif:str, down:str) -> tuple[str, str]:

        ''' Get the variable codons in the mRNA

        Distal: the codon that spans the splice junction during distal splicing
        * phase -1: "_" (non-coding)
        * phase 0: "." (in-frame)
        * phase 1: last base of the upstream exon + first 2 bases of the downstream exon
        * phase 2: last 2 bases of the upstream exon + first base of the downstream exon

        Proximal: insertion of the second NAG between the upstream and downstream exon bases
        
        Parameters
        ----------
        phase : int {-1, 0, 1, 2}
            position of the motif in the reading frame
        up : str
            last 3+ bases of the upstream exon
        motif : str
            6 splice site bases
        down : str
            first 3+ bases following the motif

        Returns
        -------
        ps : str
            variable codon bases included when motif is spliced proximally
        ds : str
            variable codon bases included when motif is spliced distally
        '''
        
        up = up.upper()
        motif = motif.upper()
        down = down.upper()

        # Non-Coding
        if phase == -1:
            return "_", "_"

        # Phase 0
        if phase == 0:
            ps = motif[-3:]
            ds = "."
        
        # Phase 1
        elif phase == 1:
            ps = f"{up[-1]}{motif[-3:]}{down[:2]}"
            ds = f"{up[-1]}{down[:2]}"
        
        # Phase 2
        elif phase == 2:
            ps = f"{up[-2:]}{motif[-3:]}{down[0]}"
            ds = f"{up[-2:]}{down[0]}"

        return ps, ds

    @staticmethod
    def variable_amino_acids(vc_ps:str, vc_ds:str) -> tuple[str, str]:

        ''' Get the variable codons in the mRNA

        Distal: the amino acid encoded by the codon that spans the splice junction during distal splicing
        * phase -1: "_" (non-coding)
        * phase 0: "." (in-frame)
        * phase 1 & 2: 2 amino acids (including STOP)

        Proximal: the amino acid encoded by the proximal variable codon(s)
        
        Parameters
        ----------
        vc_ps : str
            variable codon bases included when motif is spliced proximally
        vc_ds : str
            variable codon bases included when motif is spliced distally

        Returns
        -------
        aa_ps : str
            variable amino acid(s) included when motif is spliced proximally
        aa_ds : str
            variable amino acid included when motif is spliced distally
        '''
        
        vc_ps = vc_ps.upper()
        vc_ds = vc_ds.upper()
        
        # Non-Coding
        if (vc_ps == "_") and (vc_ds == "_"):
            return "_", "_"

        # Coding
        aa_ps = Seq.translate(vc_ps)

        # Phase 0
        if vc_ds == ".":
            aa_ds = "."
        
        # Phases 1 & 2
        else:
            aa_ds = Seq.translate(vc_ds)

        return aa_ps, aa_ds

    @staticmethod
    def transition_type(ps:str, ds:str) -> str:

        ps = ps.upper()
        ds = ds.upper()
        
        ''' Get the amino acid transition type
        
        Phase 0 Transition Types
        * E the second NAG is 'GAG', causing the insertion of a glutamic acid during proximal splicing 
        * Q the second NAG is 'CAG', causing the insertion of a glutamine during proximal splicing 
        * K the second NAG is 'AAG', causing the insertion of a lysine during proximal splicing 
        * *: the second NAG is 'TAG', causing the insertion of a STOP during proximal splicing 

        Phase 1 & 2 Transition Types
        * Duplicate Indel (DID): both proximal amino acids are the same as the distal
        * Pre-Indel (NID): only the N-terminal (upstream) proximal amino acid is the same as the distal
        * Post-Indel (CID): only the C-terminal (downstream) proximal amino acid is the same as the distal
        * Indel + Replacement (IDR): neither proximal amino acid is the same as the distal
        * No Change (NC): the N-terminal proximal and distal amino acids are both STOPs
        * Elongation/Truncation (ET): there exists a STOP for one splice form but not the other        

        Parameters
        ----------
        ps : str
            variable amino acid(s) included when motif is spliced proximally
        ds : str
            variable amino acid included when motif is spliced distally

        Returns
        -------
        _ : str
            amino acid transition type
        '''

        # Non-Coding
        if (ps == "_") and (ds == "_"):
            return "_"
        
        # Phase 0
        if ds == ".":
            return ps
        
        # Phases 1 & 2
        pn, pc = ps

        # No Change (NC)
        if (pn == "*") and (ds == "*"):
            return "NC"
        
        # Elongation/Truncation (ET)
        if (((pn != "*") and (pc != "*") and (ds == "*")) or 
            (((pn == "*") or (pc == "*")) and (ds != "*"))):
            return "ET"
        
        # Duplicate Indel (DID)
        if (pn == ds) and (pc == ds):
            return "DID"

        # Pre-Indel (NID)
        if (pn != ds) and (pc == ds):
            return "NID"

        # Post-Indel (CID)
        if (pn == ds) and (pc != ds):
            return "CID"

        # Indel + Replacement (IDR)
        return "IDR"

    @staticmethod
    def inserted_aa(ps:str, ds:str) -> list[str]:

        ''' Get the amino acids that are present in the proximal splice form but absent in the distal

        Parameters
        ----------
        ps : str
            variable amino acid(s) included when motif is spliced proximally
        ds : str
            variable amino acid included when motif is spliced distally

        Returns
        -------
        _ : list of str
            'inserted' amino acid(s), from distal to proximal splicing
        '''

        ps = ps.upper()
        ds = ds.upper()

        # Non-Coding
        if (ps == "_") and (ds == "_"):
            return ["_"]
        
        # Phase 0
        if ds == ".":
            return list(ps)
        
        pn, pc = ps
        
        # Duplicate Indel
        if (pn == ds) and (pc == ds):
            return list(ds)
        
        # Doubling
        if pn == pc:
            return list(ps)
        
        # STOP
        if pn == "*":
            
            # No Change
            if ds == "*":
                return ["."]
            
            # Truncation
            return ["*"]

        # Pre-Indel
        if (pn != ds) and (pc == ds):
            return [pn]
        
        # Post-Indel
        if (pn == ds) and (pc != ds):
            return [pc]
        
        # Indel + Replacement
        return list(ps)

    @staticmethod
    def deleted_aa(ps:str, ds:str) -> list[str]:

        ''' Get the amino acids that are present in the distal splice form but absent in the proximal

        Parameters
        ----------
        ps : str
            variable amino acid(s) included when motif is spliced proximally
        ds : str
            variable amino acid included when motif is spliced distally

        Returns
        -------
        _ : list of str
            'deleted' amino acid(s), from distal to proximal splicing
        '''

        ps = ps.upper()
        ds = ds.upper()
        
        # Non-Coding
        if (ps == "_") and (ds == "_"):
            return ["_"]
        
        # Phase 0
        if ds == ".":
            return ["."]
        
        # remove anything after STOP
        ps = ps.split("*")[0]
        
        # Duplicate, Pre-, Post-Indel; No Change
        if ds in ps:
            return ["."]
        
        # Indel + Replacement, Elongation/Truncation
        return list(ds)
