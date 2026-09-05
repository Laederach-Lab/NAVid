from ..lib import ROOT, pickle, requests, pd
from ..utils import Seq, Scenario, Variant

class State:

    ''' Current app state (input, variant, scenario) 
    
    Instance Variables:
    * chrom (str): current 'Chromosome' input
    * var_pos (int): current 0-index position of variant
    * alt (str): current 'Alternate Allele' input
    * ss5 (int): current 0-index position of 5' splice site
    * strand (str {"+", "-"}): current strand
    * phase (int {0, 1, 2}): current selected phase tab
    * allele (str {"ref", "alt"}): current selected allele

    Class Variables:
    * cv (DataFrame): ClinVar variants (id, chrom, pos, ref, alt)
    '''

    def __init__(self):
        self.reset()

    def reset(self):
        self.chrom = None
        self.var_pos = None
        self.alt = None
        self.ss5 = None
        self.strand = "+"
        self.index = 1
        self.phase = 0
        self.allele = "ref"
        self.example = None

    # ===== Update Values ===== #

    def update_chrom(self, chrom:str) -> str:

        ''' Update chromosome value,
            formatted as 'chr#',
            and return the new value + error message

            Return error message (as string) if invalid input:
            * ERROR: Only canonical chromosomes (1-22, X, Y, M) are supported
        '''

        # no input
        if chrom == "":
            return ""
        
        chrom = Seq.format_chrom(chrom)

        # invalid chromosome
        if not chrom:
            self.chrom = None
            return "ERROR: Only canonical chromosomes (1-22, X, Y, M) are supported"
        
        # valid input
        self.chrom = chrom
        return chrom

    def update_var_pos(self, loc:str|int) -> str:
        
        ''' Update variant location value,
            formatted as a comma-separated integer,
            and return the new value

            Return error message (as string) if invalid input:
            * ERROR: Variant position must be a non-negative integer
        '''

        # no input
        if loc == "":
            return ""
        
        # valid positions
        v = "non-negative" if (self.index == 0) else "positive"

        try:
            l = int("".join(str(loc).split(",")))
            
            if Seq.valid_coord(l, self.chrom):
                self.var_pos = l
                return f"{self.var_pos:,}"
                
        except ValueError:
            self.var_pos = None
            return f"ERROR: Variant position must be a {v} integer"

        self.var_pos = None
        return f"ERROR: Variant position must be a {v} integer"

    def update_alt_allele(self, alt:str) -> str:

        ''' Update alternate allele value,
            formatted as a single uppercase DNA base,
            and return the new value

            Return error message (as string) if invalid input:
            * ERROR: Alternate allele must be a single (A, C, G, or T)
        '''
        
        # no input
        if alt == "":
            return ""
        
        a = str(alt).upper().replace("U", "T")[0]

        # valid input
        if a in ["A", "C", "G", "T"]:
            self.alt = a
            return self.alt

        # invalid alternate allele
        self.alt = None
        return "ERROR: Alternate allele must be a single nucleotide (A, C, G, or T)"
    
    def update_ss5_pos(self, loc:str|int) -> str:
        
        ''' Update 5' splice site location value,
            formatted as a comma-separated integer,
            and return the new value

            Return error message (as string) if invalid input:
            * ERROR: 5' splice site position must be a non-negative integer
        '''

        # no input
        if loc == "":
            return ""

        # valid positions
        v = "non-negative" if (self.index == 0) else "positive"

        try:
            l = int("".join(str(loc).split(",")))
            
            # valid input
            if Seq.valid_coord(l, self.chrom):
                self.ss5 = l
                return f"{self.ss5:,}"

        # invalid position                
        except ValueError:
            self.ss5 = None
            return f"ERROR: 5' splice site position must be a {v} integer"

        self.ss5 = None
        return f"ERROR: 5' splice site position must be a {v} integer"

    def update_strand(self, strand:str) -> str|None:

        ''' Update strand '''

        # invalid input
        if strand not in ["+", "-"]:
            return "+"

        self.strand = strand

        return strand

    def update_index(self, index:int|str) -> str|None:
        
        ''' Update index value '''
        
        # invalid input
        if str(index) not in ["0", "1"]:
            return "1"
        
        self.index = int(index)

        return str(index)

    def update_phase(self, phase:str):

        ''' Update current phase '''

        if phase == "p0":
            self.phase = 0
        elif phase == "p1":
            self.phase = 1
        elif phase == "p2":
            self.phase = 2
        else:
            self.phase = None
    
    def update_allele(self, allele:bool):

        ''' Update current allele '''

        if allele:
            self.allele = "alt" 
        else:
            self.allele = "ref"

    def set_example(self, scen:pd.Series):

        ''' Set example variant '''

        db = {
            "dbSNP_common" : "dbSNP",
            "dbSNP_rare" : "dbSNP",
            "ClinVar" : "ClinVar",
        }[scen["db"]]

        self.example = f"{db} {scen["id"]}"

    def reset_example(self):
        ''' Reset example variant '''
        self.example = None

    # ===== Submit ===== #

    def get_variant(self) -> tuple[list[str], Variant]:
        
        ''' Get variant from input values

            If any values are missing or invalid,
            return empty Variant and error message(s) as a list
            * ERROR: Missing chromosome
            * ERROR: Missing alternate allele
            * ERROR: Missing variant position
            * ERROR: Missing 5' splice site position
            * ERROR: Variant position outside of chromosome
            * ERROR: 5' splice site position outside of chromosome

            Otherwise, return empty list
        '''

        errors = []

        # Check if values are set
        if self.chrom == None:
            errors.append("ERROR: Missing chromosome")
        if self.alt == None:
            errors.append("ERROR: Missing alternate allele")
        if self.var_pos == None:
            errors.append("ERROR: Missing variant position")
        # if self.ss5 == None:
        #     errors.append("ERROR: Missing 5' splice site position")

        if errors:
            return errors, Variant()

        # Check if values are biologically valid
        if self.var_pos >= Seq.CHR_LEN[self.chrom]:
            errors.append("ERROR: Variant position outside of chromosome")
        if (self.ss5 != None) and (self.ss5 >= Seq.CHR_LEN[self.chrom]):
            errors.append("ERROR: 5' splice site position outside of chromosome")

        if errors:
            return errors, Variant()
        
        # update variant
        vnt = Variant(self.chrom, self.var_pos, self.alt, self.ss5, self.strand, self.index)

        return [], vnt

    def get_scenario(self, vnt:Variant) -> Scenario:

        ''' Set current scenario '''

        if self.allele == "ref":
            if self.phase == 0: return vnt.ref_scen0
            if self.phase == 1: return vnt.ref_scen1
            if self.phase == 2: return vnt.ref_scen2

        if self.allele == "alt":
            if self.phase == 0: return vnt.alt_scen0
            if self.phase == 1: return vnt.alt_scen1
            if self.phase == 2: return vnt.alt_scen2

        return Scenario()

    # ===== Search ===== #

    @staticmethod
    def cv() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/dbs/ClinVar_variants", "rb"))
    
    @staticmethod
    def find_variant(db:str, acc_num:str) -> tuple[tuple[str, int, list[str], str]|None, str|None]:

        ''' Search for a pre-existing variant
        
        Parameters
        ----------
        db : str {"dbSNP", "ClinVar"}
            variant database
        acc_num : str
            variant accession number

        Returns
        -------
        _ : tuple | None
            variant information: chromosome, coordinate, alternate allele, and upstream exon position
        _ : str | None
            error message
        '''
        
        if not acc_num:
            return None, None

        acc_num = str(acc_num).strip()

        # dbSNP
        if db == "dbSNP":

            # API Call
            api = requests.get(f"https://clinicaltables.nlm.nih.gov/api/snps/v3/search?terms={acc_num}")
            if api.status_code != 200:
                return None, "Error"
            
            results = api.json()

            # Variant not found
            if acc_num not in results[1]:
                return None, "Variant not found"
            
            info = results[3][results[1].index(acc_num)]
            vnt = {
                "chrom" : f"chr{info[1].replace("MT", "M")}",
                "pos" : int(info[2]),
                "strand" : "+",
                "ref" : info[3].split("/")[0],
                "alt" : ",".join([s.split("/")[1] for s in info[3].split(",")])
            }
        
        # ClinVar
        elif db == "ClinVar":

            # API Call
            api = requests.get(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={acc_num}&retmode=json")
            
            if api.status_code != 200:
                return None, "Error"
            
            results = api.json()

            # Variant not found
            if acc_num not in results["result"]["uids"]:
                return None, "Variant not found"

            for dicts in results["result"][acc_num]["variation_set"][0]["variation_loc"]:
                if dicts["assembly_name"] == "GRCh38":
                    loc = dicts
                    break
            delta = str(results["result"][acc_num]["variation_set"][0]["canonical_spdi"]).split(":")
            vnt = {
                "chrom" : f"chr{loc["chr"].replace("MT", "M")}",
                "pos" : int(loc["start"]) - 1,
                "strand" : "+",
                "ref" : delta[-2],
                "alt" : delta[-1]
            }

        # Not a single-base substitution
        alt = [a.strip().upper() for a in str(vnt["alt"]).split(",") if len(a.strip()) == 1]

        if (len(vnt["ref"]) != 1) or (len(alt) == 0):
            return None, "Variant is not a single-base substitution"
        
        chrom = vnt["chrom"]
        pos = vnt["pos"]
        strand = vnt["strand"]

        return (chrom, pos, alt, strand), None
