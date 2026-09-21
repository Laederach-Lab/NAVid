from ..lib import ROOT, pickle, re, BioSeq

class Seq:

    ''' Genomic sequences '''
    
    NAG_STRENGTH = pickle.load(open(f"{ROOT}/data/stats/nag_strength", "rb"))

    @staticmethod
    def hg38(chrom:str=None) -> BioSeq | dict[str, BioSeq]:

        ''' Get genomic sequences for hg38 reference genome
        
        chrom:
        * True: get entire hg38 dictionary {chromosome, sequence}
        * False: get dictionary of non-canonical chromosomes
        * 'chr#': get sequence for chromosome #
        '''

        # format chromosome (key)
        chrom = Seq.format_chrom(chrom)

        # invalid chromosome (return dictionary)
        if chrom == None:
            return pickle.load(open(f"{ROOT}/data/hg38", "rb"))

        # return chromosome
        return pickle.load(open(f"{ROOT}/data/seq/{chrom}", "rb"))
    
    CHR_LEN = {f"chr{c}" : len(pickle.load(open(f"{ROOT}/data/seq/chr{c}", "rb")))
               for c in list(map(str, range(1, 23))) + ["X", "Y", "M"]}
    ''' Lengths of canonical hg38 chromosomes '''

    @staticmethod
    def gseq(chrom:str, start:int, end:int|None=None, strand:str|None=None, upper:bool=True) -> str|None:

        ''' Get genomic sequence

        Parameters
        ----------
        chrom : str
            chromosome, formatted as chr#
        start : int
            0-index position of the first (inclusive) base
        end : int (default = None)
            0-index position of the last (exclusive) base; 
            if None, end = start + 1
        strand : str {"+", "-"} (default = "+")
            genomic strand
        upper : bool (default = True)
            whether to convert the sequence to uppercase, thus ignoring repeat masker indications

        Returns
        -------
        None
            invalid coordinates
        _ : str
            genomic sequence
        '''
        
        # Validate arguments
        chrom = Seq.format_chrom(chrom)
        if not chrom:
            return

        if not Seq.valid_coord(start, chrom):
            return None

        end = end if end else start+1
        if (not Seq.valid_coord(start, chrom)) or (end < start) or (not Seq.valid_coord(end, chrom)):
            return None

        if not Seq.valid_strand(strand):
            return None

        # Get sequence
        gs = Seq.hg38(chrom)[start:end]

        # strand
        if strand == "-":
            gs = gs.reverse_complement()
        
        # repeat masker
        if upper:
            gs = gs.upper()

        return str(gs).upper()

    @staticmethod
    def translate(seq:str) -> str:

        ''' Translate mRNA to amino acids '''

        return str(BioSeq.translate(seq)) if seq else None
    
    @staticmethod
    def rc(seq:str) -> str:

        ''' Return reverse complement '''

        return str(BioSeq(seq).reverse_complement())

    # ===== Validate Values ===== #

    @staticmethod
    def format_chrom(chrom:str) -> str|None:

        ''' Returns chromosome formatted as 'chr#'
        
        if invalid chromosome, return None
        '''

        c = str(chrom).strip().upper().strip("CHR")
        if c in list(map(str, range(1, 23))) + ["X", "Y", "M"]:
            return f"chr{c}"
        
        return None

    @staticmethod
    def valid_coord(coord:int|str, chrom:str=None) -> bool:

        ''' Returns whether coordinate is valid:
            non-negative integer within chromosome (if specified)
        '''

        try:
            coord = int(coord)
        except:
            return False
        
        return (coord >= 0) and ((chrom == None) or ((chrom in Seq.CHR_LEN) and (coord < Seq.CHR_LEN[chrom])))

    @staticmethod
    def valid_range(start:int, end:int, chrom:str=None) -> bool:

        ''' Returns whether a coordinate range is valid:
            ordered (ascending), non-negative integers within chromosome (if specified)
        '''
        
        if (start == None) and (end == None):
            return True
        
        if (end == None):
            return Seq.valid_coord(start, chrom)

        if (start == None):
            return False

        return Seq.valid_coord(start, chrom) and Seq.valid_coord(end, chrom) and (start < end)

    @staticmethod
    def valid_seq(seq:str) -> bool:
        
        ''' Returns whether sequence is valid (string of only A, C, G, T, N) '''
        
        return isinstance(seq, str) and re.fullmatch("[acgtun]+", seq, re.IGNORECASE)
    
    @staticmethod
    def format_base(base:str, dna:bool=True) -> str|None:

        ''' Returns formatted base: (A, C, G, T/U)
        
        if {base} is not a single nucleotide or not A, C, G, T, or U, return None
        '''

        if not Seq.valid_seq(base) or len(base) != 1:
            return None
        
        b = base.strip().upper()

        if dna:
            return b.replace("U", "T")
        
        return b.replace("T", "U")

    @staticmethod
    def valid_strand(strand:str) -> bool:

        ''' Returns whether strand is valid (+, -, None)'''

        return strand in ["+", "-", None]
    
    # ===== Coordinate Format ===== #

    @staticmethod
    def str2coord(coord:str, index:int=0) -> tuple[str, int, str]:

        ''' Parse genomic coordinates

        Parameters
        ----------
        coord : str
            coordinate, formatted as chr#:##[-##][+/-]
        index : int {0, 1} (default = 0)
            indexing used for chromosome position

        Returns
        -------
        chrom : str
            chromosome, formatted as chr#
        start : int
            0-index position of the first (inclusive) base
        end : int
            0-index position of the last (exclusive) base;
            if not a range, end = start + 1
        strand : str {"+", "-"}
            genomic strand \\
            if not specified, use None
        '''

        # chromosome
        chrom = Seq.format_chrom(coord.split(":")[0])

        # position
        if ":" in coord:
            pos = coord.split(":")[1].removesuffix("+").removesuffix("-").split("-")
            start = int("".join(pos[0].split(","))) - index
            end = int("".join(pos[1].split(","))) if len(pos) > 1 else start + 1
        else:
            start = None
            end = None
        
        # strand
        strand = coord[-1] if coord.endswith(("+", "-")) else None

        return chrom, start, end, strand

    @staticmethod
    def coord2str(chrom:str|None, start:int|str|None=None, end:int|None=None, strand:str|None=None) -> str|None:

        ''' Convert chromosomal coordinates to string
        
        Parameters
        ----------
        chrom : str
            chromosome, formatted as chr#
        start : int
            0-index position of the first (inclusive) base
        end : int (default = None)
            0-index position of the last (exclusive) base for a genomic range; 
            if not specified, format as a single nucleotide
        strand : str {"+", "-"} (default = None)
            genomic strand

        Returns
        -------
        _ : str
            coordinates, formatted as [chrom]:[start][-end][strand]
        '''
        
        c = Seq.format_chrom(chrom)
        d = strand if strand in ["+", "-"] else ""

        if start != None:
            end = None if (end == None) or (end == start + 1) else end
        
        s = f"{start:,}" if start else ""
        e = f"-{end:,}" if end else ""
        
        if c and Seq.valid_range(start, end, chrom) and Seq.valid_strand(strand):
            return f"{c}{f":{s}{e}" if s or e else ""}{d}"
        
        return None

    @staticmethod
    def vnt2str(chrom:str|None, pos:int|str|None, strand:str|None, ref:str|None, alt:str|None) -> str|None:

        ''' Convert chromosomal coordinates to variant definition
        
        Parameters
        ----------
        chrom : str
            chromosome, formatted as chr#
        pos : int
            0-index position variant
        strand : str {"+", "-"}
            genomic strand
        ref, alt : str {"A", "C", "G", "T"}
            reference and alternative bases

        Returns
        -------
        _ : str
            variant, formatted as [chrom]:[pos][strand][ref]>[alt]
        '''

        # get location
        coord = Seq.coord2str(chrom=chrom, start=pos, strand=strand)

        if not coord:
            return None
        
        return f"{coord}{ref}>{alt}"

    # ===== Variant Classification ===== #

    @staticmethod
    def sbs_type(ref:str|None, alt:str|None) -> str|None:

        ''' Get single base substitution type
        
        Parameters
        ----------
        ref, alt : str {"A", "C", "G", "T"}
            reference and alternative bases

        Returns
        -------
        _ : str {"Transition", "Transversion"}
            variant type (based on purine/pyrimidine classification)
        '''

        ref = Seq.format_base(ref)
        alt = Seq.format_base(alt)

        if (not ref) or (not alt):
            return None
        
        if (ref == "N") or (alt == "N"):
            return "N/A"
        
        if Seq.nt[ref].type == Seq.nt[alt].type:
            return "Transition"
        
        return "Transversion"

    # ===== Names ===== #

    class Nucleotide:

        def __init__(self, name:str, code:str, pair:str, cat:str):
            self.name = name
            self.code = code
            self.pair = pair
            self.type = cat

    nt = {
        "A" : Nucleotide("Adenine", "A", "T", "Purine"),
        "C" : Nucleotide("Cytosine", "C", "G", "Pyrimidine"),
        "G" : Nucleotide("Guanine", "G", "C", "Purine"),
        "T" : Nucleotide("Thymine", "T", "A", "Pyrimidine"),
        "U" : Nucleotide("Uracil", "U", "A", "Pyrimidine"),
        "N" : Nucleotide("Nucleotide", "N", "N", "None")
    }

    class AminoAcid:
        
        ''' Amino acid names
        
        Properties:
        * `name`: amino acid name {e.g., Alanine}
        * `code3`: 3-letter code {e.g., Ala}
        * `code1`: 1-letter code {e.g., A}
        * `full_name`: `name` (`code3`, `code1`) {e.g., Alanine (Ala, A)}
        * `name_code3`: `name` (`code3`) {e.g., Alanine (Ala)}
        * `name_code1`: `name` (`code1`) {e.g., Alanine (A)}
        * `codes`: `code3` (`code1`) {e.g., Ala (a)}
        '''

        def __init__(self, name:str, code3:str, code1:str):
            self.name = name
            self.code3 = code3
            self.code1 = code1

            self.null = (code1 == ".") or (code1 == "_")

            # e.g., Alanine (Ala, A)
            self.full_name = f"{self.name} ({self.code3}, {self.code1})" if not self.null else self.code1
    
            # e.g., Alanine (Ala)
            self.name_code3 = f"{self.name} ({self.code3})" if not self.null else self.code1

            # e.g., Alanine (A)
            self.name_code1 = f"{self.name} ({self.code1})" if not self.null else self.code1

            # e.g., Ala (A)
            self.codes = f"{self.code3} ({self.code1})" if not self.null else self.code1

    aa = {
        "G" : AminoAcid("Glycine", "Gly", "G"),
        "A" : AminoAcid("Alanine", "Ala", "A"),
        "V" : AminoAcid("Valine", "Val", "V"),
        "L" : AminoAcid("Leucine", "Leu", "L"),
        "I" : AminoAcid("Isoleucine", "Ile", "I"),
        "M" : AminoAcid("Methionine", "Met", "M"),
        "F" : AminoAcid("Phenylalanine", "Phe", "F"),
        "W" : AminoAcid("Tryptophan", "Trp", "W"),
        "P" : AminoAcid("Proline", "Pro", "P"),
        "S" : AminoAcid("Serine", "Ser", "S"),
        "T" : AminoAcid("Threonine", "Thr", "T"),
        "C" : AminoAcid("Cysteine", "Cys", "C"),
        "Y" : AminoAcid("Tyrosine", "Tyr", "Y"),
        "N" : AminoAcid("Asparagine", "Asn", "N"),
        "Q" : AminoAcid("Glutamine", "Gln", "Q"),
        "D" : AminoAcid("Aspartic Acid", "Asp", "D"),
        "E" : AminoAcid("Glutamic Acid", "Glu", "E"),
        "K" : AminoAcid("Lysine", "Lys", "K"),
        "R" : AminoAcid("Arginine", "Arg", "R"),
        "H" : AminoAcid("Histidine", "His", "H"),
        "*" : AminoAcid("STOP", "STOP", "*"),
        "." : AminoAcid("None", ".", "."),
        "_" : AminoAcid("Non-Coding", "_", "_")
    }

    tsn = {
        "E" : aa["E"].codes,
        "Q" : aa["Q"].codes,
        "K" : aa["K"].codes,
        "*" : aa["*"].codes,
        "DID" : "Duplicate Indel (DID)",
        "NID" : "Pre-Indel (NID)",
        "CID" : "Post-Indel (CID)",
        "IDR" : "Indel + Replacement (IDR)",
        "NC" : "No Change (NC)",
        "ET" : "Elongation/Truncation (ET)",
        "_" : "None (non-coding)",
    }
