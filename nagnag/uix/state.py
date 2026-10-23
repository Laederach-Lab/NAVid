from ..lib import requests, pd
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

        # reset
        if (chrom == "") or (chrom == None):
            self.chrom = None
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

        # reset
        if (loc == "") or (loc == None):
            self.var_pos = None
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
        
        # reset
        if (alt == "") or (alt == None):
            self.alt = None
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

        # reset
        if (loc == "") or (loc == None):
            self.ss5 = None
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
            return self.strand

        self.strand = strand

        return strand

    def update_index(self, index:int|str) -> str|None:
        
        ''' Update index value '''
        
        # invalid input
        if str(index) not in ["0", "1"]:
            return str(self.index)
        
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
            * Missing chromosome
            * Missing alternate allele
            * Missing variant position
            * Variant location is outside of the chromosome
            * 5' splice site location is outside of the chromosome

            Otherwise, return empty list
        '''

        errors = []

        # Check if values are set
        if self.chrom == None:
            errors.append(r"Missing <b>Chromosome</b>")
        if self.alt == None:
            errors.append(r"Missing <b>Alternate Allele</b>")
        if self.var_pos == None:
            errors.append(r"Missing <b>Variant Location</b>")

        if errors:
            return errors, Variant()

        # Check if values are biologically valid
        if self.var_pos >= Seq.CHR_LEN[self.chrom]:
            errors.append(r"<b>Variant Location</b> is outside of the chromosome")
        if (self.ss5 != None) and (self.ss5 >= Seq.CHR_LEN[self.chrom]):
            errors.append(r"<5' Splice Site Location</b> is outside of the chromosome")

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
    def call_api(url:str, params:dict=None, timeout:int=30) -> tuple[bool, dict] | tuple[bool, str]:

        ''' Evaluate API call
        
        Parameters
        ----------
        url : str
            API URL

        Returns
        -------
        success : bool
            was the API call successful (no exceptions)?
        results : str | dict
            error message (success = False) or JSON of API request (success = True) 
        '''

        # request
        try:
            api = requests.get(url=url, params=params, timeout=timeout)
        except:
            return False, "Error"
        
        # exceptions
        try:
            api.raise_for_status()
        except requests.exceptions.Timeout:
            return False, f"Error: connection timeout (status {api.status_code})"
        except requests.exceptions.ConnectionError:
            return False, f"Error: connection failed (status {api.status_code})"
        except:
            return f"Error: status {api.status_code}"

        # JSON
        try:
            results = api.json()
        except:
            return False, "Error: invalid JSON"

        return True, results

    @staticmethod
    def find_variant(db:str, acc_num:str) -> tuple[tuple[str, int, list[str], str], str] | \
                                             tuple[tuple[str, int, list[str], str], None] | \
                                             tuple[None, str] | tuple[None, None]:

        ''' Search for a pre-existing variant
        
        Parameters
        ----------
        db : str {"dbSNP", "ClinVar"}
            variant database
        acc_num : str
            variant accession number

        Returns
        -------
        vnt : tuple | None
            variant information: chromosome, coordinate, alternate allele, and strand
        msg : str | None
            error message
        '''
        
        if not acc_num:
            return None, None

        acc_num = str(acc_num).strip()

        # define ID
        if db == "dbSNP":
            uid = acc_num.removeprefix("rs")
            params = dict(db="snp", retmode="json", id=uid)

        elif db == "ClinVar":
            uid = acc_num.split(".")[0].removeprefix("VCV").lstrip("0")
            params = dict(db="clinvar", retmode="json", id=uid)

        # API Call
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
        success, results = State.call_api(url=url, params=params)

        if not success:
            return None, results

        # Information
        try:

            if db == "dbSNP":
                info = results["result"][uid]
                vnt_type_key = "snp_class"
                vnt_type_value = "snv"
                spdi_key = "spdi"

                # coordinates
                chrom, pos = str(info["chrpos"]).split(":")

            elif db == "ClinVar":
                info = results["result"][uid]["variation_set"][0]
                vnt_type_key = "variant_type"
                vnt_type_value = "single nucleotide variant"
                spdi_key = "canonical_spdi"

                # GRCh38/hg38 coordinates
                for vdict in info["variation_loc"]:
                    if vdict["assembly_name"] == "GRCh38":
                        chrom = vdict["chr"]
                        pos = vdict["start"]
                        break

            # non-SNP
            if info[vnt_type_key] != vnt_type_value:
                return None, "Variant is not a single-base substitution"

            # alleles
            spdi = str(info[spdi_key]).split(",")

            vnt = {
                "chrom" : f"chr{chrom.replace("MT", "M")}",
                "pos" : int(pos) - 1,
                "strand" : "+",
                "ref" : ",".join(set(a.split(":")[2] for a in spdi)),
                "alt" : ",".join(set(a.split(":")[3] for a in spdi)),
            }

        # variant not found
        except:
            return None, "Variant not found"

        # Non-canonical chromosome
        if not Seq.format_chrom(vnt["chrom"]):
            return None, "Variant not on a canonical chromosome"
            
        # Not a single-base substitution
        alt = [a.strip().upper() for a in str(vnt["alt"]).split(",") if len(a.strip()) == 1]

        if (len(vnt["ref"]) != 1) or (len(alt) == 0):
            return None, "Variant is not a single-base substitution"

        # Return        
        chrom = vnt["chrom"]
        pos = vnt["pos"]
        strand = vnt["strand"]

        return (chrom, pos, alt, strand), None
