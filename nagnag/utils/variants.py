from ..lib import re
from ..utils import Seq, Scenario

class Variant:

    ''' NAGNAG-affecting variants '''

    def __init__(self, chrom:str=None, var_coord:int=None, alt:str=None, ss5_coord:int=None, strand:str=None, index:int=0):

        ''' NAGNAG-affecting variant

        Parameters
        ----------
        chrom : str
            chromosome, formatted as chr#
        ss5_coord : int
            0-index position of the end of the upstream exon;
            if on the plus strand, this is the first base of the intron;
            if on the minus strand, this is the last base of the exon
        var_coord : int
            0-index position of variant
        strand : str {"+", "-"}
            genomic strand
        alt : str {"A", "C", "G", "T"}
            alternate allele base (on the proper strand)
        index : int {0, 1} (default = 0)
            first position on chromosome
        '''

        self.reset()
        self.set_variant(chrom, var_coord, alt, ss5_coord, strand, index)

    def reset(self):

        # Default Attributes
        self.empty = True
        self.not_nagnag = True
        self.not_scenario = True

        # coordinates
        self.chrom = None
        self.ss5_coord = None
        self.motif_start = None
        self.motif_end = None
        self.var_coord = None
        self.down_start = None
        self.strand = None

        # sequences
        self.var_pos = None
        self.ref_motif = None
        self.alt_motif = None
        self.up_seq = None
        self.down_seq = None
        self.eflank_5ss = None
        self.iflank_5ss = None
        self.iflank_3ss = None
        self.eflank_3ss = None

        self.ref = None
        self.alt = None
        self.effect = None

        # proteome effects
        self.ref_scen0 = Scenario()
        self.ref_scen1 = Scenario()
        self.ref_scen2 = Scenario()

        self.alt_scen0 = Scenario()
        self.alt_scen1 = Scenario()
        self.alt_scen2 = Scenario()
    
    def set_variant(self, chrom:str=None, var_coord:int=None, alt:str=None, ss5_coord:int=None, strand:str=None, index:int=0):

        ''' Set variant information
        
        Parameters
        ----------
        chrom : str
            chromosome, formatted as chr#
        ss5_coord : int
            position of the end of the upstream exon;
            if on the plus strand, this is the first base of the intron;
            if on the minus strand, this is the last base of the exon
        var_coord : int
            position of variant
        strand : str {"+", "-"}
            genomic strand
        alt : str {"A", "C", "G", "T"}
            alternate allele base (on the proper strand)
        index : int {0, 1} (default = 0)
            first position on chromosome
        '''

        # Verify Variant arguments
        chrom = Seq.format_chrom(chrom) if (chrom != None) else self.chrom
        if not chrom:
            return
        
        var_coord = var_coord if (var_coord != None) else self.var_coord
        if not Seq.valid_coord(var_coord, chrom):
            return
        var_coord -= index

        strand = strand if (strand != None) else self.strand
        if strand not in ["+", "-"]:
            return
        
        alt = alt if (alt != None) else self.alt
        if alt.upper().replace("U", "T") not in ["A", "C", "G", "T"]:
            return

        # Valid Variant
        self.empty = False

        # Coordinates
        self.chrom = chrom
        self.var_coord = var_coord
        self.strand = strand

        self.ref = Seq.gseq(chrom=chrom, start=var_coord, strand=strand)
        self.alt = alt

        # Find NAGNAG
        nn = Variant.affected_nagnag(chrom, var_coord, strand, alt)
        
        # No NAGNAG
        if not nn:
            return
        
        # NAGNAG-Affecting
        self.not_nagnag = False

        # Location
        self.motif_start = nn[0]
        self.motif_end = nn[1]
        self.down_start = nn[1]

        # Motif Sequences
        self.var_pos = nn[2]
        self.ref_motif = nn[3]
        self.alt_motif = nn[4]

        self.effect = Variant.nagnag_effect(self.ref_motif, self.alt_motif)

        # Verify 5' splice site
        ss5_coord = ss5_coord if (ss5_coord != None) else self.ss5_coord
        if not Seq.valid_coord(ss5_coord, chrom):
            return
        ss5_coord -= index

        # Valid Scenario
        self.not_scenario = False
        self.ss5_coord = ss5_coord

        # Region Sequences
        flank = 50
        if strand == "+":
            self.up_seq = Seq.gseq(chrom=chrom, start=ss5_coord-3, end=ss5_coord, strand=strand)
            self.down_seq = Seq.gseq(chrom=chrom, start=self.down_start, end=self.down_start+3, strand=strand)

            self.eflank_5ss = Seq.gseq(chrom=chrom, start=ss5_coord-flank, end=ss5_coord, strand=strand)
            self.iflank_5ss = Seq.gseq(chrom=chrom, start=ss5_coord, end=ss5_coord+flank, strand=strand)
            self.iflank_3ss = Seq.gseq(chrom=chrom, start=self.down_start-flank, end=self.down_start, strand=strand)
            self.eflank_3ss = Seq.gseq(chrom=chrom, start=self.down_start, end=self.down_start+flank, strand=strand)
        elif strand == "-":
            self.up_seq = Seq.gseq(chrom=chrom, start=ss5_coord, end=ss5_coord+3, strand=strand)
            self.down_seq = Seq.gseq(chrom=chrom, start=self.down_start-3, end=self.down_start, strand=strand)

            self.eflank_5ss = Seq.gseq(chrom=chrom, start=ss5_coord, end=ss5_coord+flank, strand=strand)
            self.iflank_5ss = Seq.gseq(chrom=chrom, start=ss5_coord-flank, end=ss5_coord, strand=strand)
            self.iflank_3ss = Seq.gseq(chrom=chrom, start=self.down_start, end=self.down_start+flank, strand=strand)
            self.eflank_3ss = Seq.gseq(chrom=chrom, start=self.down_start-flank, end=self.down_start, strand=strand)
        
        # Reference Allele Proteome Effects
        self.ref_scen0 = Scenario(0, self.eflank_5ss, self.ref_motif, self.eflank_3ss, "ref")
        self.ref_scen1 = Scenario(1, self.eflank_5ss, self.ref_motif, self.eflank_3ss, "ref")
        self.ref_scen2 = Scenario(2, self.eflank_5ss, self.ref_motif, self.eflank_3ss, "ref")

        # Alternate Allele Proteome Effects
        self.alt_scen0 = Scenario(0, self.eflank_5ss, self.alt_motif, self.eflank_3ss, "alt")
        self.alt_scen1 = Scenario(1, self.eflank_5ss, self.alt_motif, self.eflank_3ss, "alt")
        self.alt_scen2 = Scenario(2, self.eflank_5ss, self.alt_motif, self.eflank_3ss, "alt")

    def __str__(self):
        
        c = f"Variant: {Seq.coord2str(chrom=self.chrom, start=self.var_coord, strand=self.strand)}"

        if self.ref_motif and self.ref:
            r = self.ref_motif.lower()[:self.var_pos] + self.ref + self.ref_motif.lower()[self.var_pos+1:]
        else:
            r = "Ref"

        if self.alt_motif and self.alt:
            a = self.alt_motif.lower()[:self.var_pos] + self.alt + self.alt_motif.lower()[self.var_pos+1:]
        else:
            a = "Alt"

        return f"{c}\n{r} --> {a}"

    @staticmethod
    def affected_nagnag(chrom:str, pos:int|str, strand:str, alt:str) -> tuple[int, int, int, str, str]|None|bool:

        ''' Get the NAGNAG affected by a specified variant
        
        Parameters
        ----------
        chrom : str
            chromosome, formatted as chr#
        pos : int
            0-index position of variant
        strand : str {"+", "-"}
            genomic strand
        alt : str {"A", "C", "G", "T"}
            alternate allele base (on the proper strand)

        Returns
        -------
        None
            the variant does not affect a NAGNAG
        gs : int
            0-index chromosomal position of the lower-index end of the motif
        ge : int
            0-index chromosomal position of the higher-index end of the motif
        vp : int {0, 1, 2, 3, 4, 5}
            position of variant in motif
        ref_motif : str
            motif sequence for the reference allele
        alt_motif : str
            motif sequence for the alternate allele
        '''

        try:
            pos = int(pos)
        except:
            return False

        # Surrounding sequences
        ref_region = Seq.gseq(chrom=chrom, start=pos-8, end=pos+9, strand=strand)
        
        # Invalid coordinates
        if not ref_region:
            return False
        
        # Not a variant
        if alt.upper() == ref_region[8].upper():
            return False

        alt_region = ref_region[:8] + alt + ref_region[9:]

        # Search for NAGNAGs
        xn = "((.[^A]G)|(.A[^G])|(.[^A][^G]))"
        nn_re = f"{xn}.AG.AG{xn}"
        nnr = re.search(nn_re, ref_region, re.IGNORECASE)
        nna = re.search(nn_re, alt_region, re.IGNORECASE)

        # Not NAGNAG-affecting
        if (not nnr) and (not nna):
            return None

        # Coordinates
        ms,me = nnr.span() if nnr else nna.span()
        ms += 3
        me -= 3

        # Variant position
        vp = 8 - ms

        # Genomic coordinates
        gs = pos + (ms - 8)
        ge = pos + (me - 8)

        # Sequences
        ref_motif = ref_region[ms:me]
        alt_motif = alt_region[ms:me]

        return gs, ge, vp, ref_motif, alt_motif

    @staticmethod
    def nagnag_effect(ref:str, alt:str) -> str|None:

        ''' Get the effect of a variant on a NAGNAG 
        
        Returns None if invalid

        Parameters
        ----------
        ref : str
            motif sequence for the reference allele
        alt : str
            motif sequence for the alternate allele

        Returns
        -------
        None
            the variant does not affect a NAGNAG
        _ : str {"create", "alter", "destroy"}
            effect of variant on NAGNAG
        '''

        ref = ref.upper()
        alt = alt.upper()
        
        # Validate
        if (len(ref) != 6) or (len(alt) != 6):
            return None
        
        # is a NAGNAG
        nnr = re.fullmatch(".AG.AG", ref)
        nna = re.fullmatch(".AG.AG", alt)

        # NAGNAG-Creating
        if (not nnr) and nna:
            return "create"
        
        # NAGNAG-Altering
        if nnr and nna and (ref != alt):
            return "alter"
        
        # NAGNAG-Destroying
        if nnr and (not nna):
            return "destroy"

        return None
