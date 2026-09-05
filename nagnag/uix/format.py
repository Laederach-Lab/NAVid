class Fmt:

    ''' App label formatting '''

    def __init__(self):
        pass
    
    @staticmethod
    def effs(key:str="", default:str="", abbr:bool=False) -> dict[str, str]|list[str]|str:
        
        ''' Variant effects (create, alter, destroy) '''

        d = {
            False : {
                "create" : "NAGNAG-Creating",
                "alter" : "NAGNAG-Altering",
                "destroy" : "NAGNAG-Destroying"
            },
            True : {
                "create" : "Create", 
                "alter" : "Alter",
                "destroy" : "Destroy"
            }
        }[abbr]

        if not key:
            return list(d.items())
        
        if key.lower() == "d":
            return d

        if key == "k":
            return list(d.keys())

        if key == "v":
            return list(d.values())
        
        if key in d:
            return d[key]
        
        return default

    @staticmethod
    def scen_effs(key:str="", default:str="", abbr:bool=False) -> dict[str, str]|list[str]|str:
        
        ''' Scenario effects (create_alt, alter_ref, alter_alt, destroy_ref) '''
        
        d = {
            False : {
                "create_alt" : "NAGNAG-Creating (Alternate Allele)", 
                "alter_ref" : "NAGNAG-Altering (Reference Allele)",
                "alter_alt" : "NAGNAG-Altering (Alternate Allele)", 
                "destroy_ref" : "NAGNAG-Destroying (Reference Allele)"
            },
            True : {
                "create_alt" : "Create", 
                "alter_ref" : "Alter (Ref)",
                "alter_alt" : "Alter (Alt)", 
                "destroy_ref" : "Destroy"
            }
        }[abbr]

        if not key:
            return list(d.items())
        
        if key.lower() == "d":
            return d
        
        if key.lower() == "k":
            return list(d.keys())

        if key.lower() == "v":
            return list(d.values())
        
        if key in d:
            return d[key]
        
        return default

    @staticmethod
    def dbs(key:str="", default:str="", vnt:bool=False, prop:bool=True) -> dict[str, str]|list[str]|str:
        
        ''' Databases (exp_stoch, exp_splice, canon, dbSNP_common, dbSNP_rare, CLinVar) '''

        d = {
            (False, True) : {
                "exp_stoch" : "Expected (Stochastic)",
                "exp_splice" : "Expected (Splice Selection)",
                "canon" : "hg38 Reference",
                "dbSNP_common" : "dbSNP Common",
                "dbSNP_rare" : "dbSNP Rare",
                "ClinVar" : "ClinVar",
            },
            (False, False) : {
                "canon" : "hg38 Reference",
                "dbSNP_common" : "dbSNP Common",
                "dbSNP_rare" : "dbSNP Rare",
                "ClinVar" : "ClinVar",
            },
            (True, True) : {
                "dbSNP_common" : "dbSNP Common",
                "dbSNP_rare" : "dbSNP Rare",
                "ClinVar" : "ClinVar",
            },
            (True, False) : {
                "dbSNP_common" : "dbSNP Common",
                "dbSNP_rare" : "dbSNP Rare",
                "ClinVar" : "ClinVar",
            }
        }[(vnt, prop)]

        if not key:
            return list(d.items())
        
        if key.lower() == "d":
            return d
        
        if key.lower() == "k":
            return list(d.keys())

        if key.lower() == "v":
            return list(d.values())
        
        if key in d:
            return d[key]
        
        return default

    @staticmethod
    def metrics(key:str="", default:str="", phases:bool=False, inds:bool=False) -> dict[str, str|list[str]]|list[str|list[str]]|str|list[str]:
        
        ''' Amino acid transition metrics (vc, aat, aatt, ins, del) '''
        
        d = {
            (False, False) : {
                "vc" : "Variable Codons",
                "aat" : "Variable Amino Acids",
                "aatt" : "Transition Type",
                "ins" : "Insertions",
                "del" : "Deletions"
            },
            (True, False) : {
                "vc" : ["p0", "p1", "p2"],
                "aat" : ["p0", "p1", "p2"],
                "aatt" : ["p0", "p1", "p2"],
                "ins" : ["p0", "p1", "p2", "all"],
                "del" : ["p1", "p2", "all"]
            },
            (True, True) : {
                "vc" : [0, 1, 2],
                "aat" : [0, 1, 2],
                "aatt" : [0, 1, 2],
                "ins" : [0, 1, 2, 3],
                "del" : [1, 2, 3],
            },
            (False, True) : {}
        }[(phases, inds)]

        if not key:
            return list(d.items())
        
        if key.lower() == "d":
            return d
        
        if key.lower() == "k":
            return list(d.keys())

        if key.lower() == "v":
            return list(d.values())
        
        if key in d:
            return d[key]
        
        return default

    @staticmethod
    def phases(key:str="", default:str="", ind:bool=False) -> dict[str, str|int]|list[str|int]|str|int:
        
        ''' Phases (nc, p0, p1, p2, all) '''
        
        d = {
            False : {
                "nc" : "Non-Coding",
                "p0" : "Phase 0",
                "p1" : "Phase 1",
                "p2" : "Phase 2",
                "all" : "All Phases"
            },
            True : {
                "nc" : -1,
                "p0" : 0,
                "p1" : 1,
                "p2" : 2,
                "all" : 3,
            }
        }[ind]

        if not key:
            return list(d.items())
        
        if key.lower() == "d":
            return d
        
        if key.lower() == "k":
            return list(d.keys())

        if key.lower() == "v":
            return list(d.values())
        
        if key in d:
            return d[key]
        
        return default

    col_names = {
        "ref" : "Reference",
        "exp_stoch" : "Expected (Stochastic)",
        "exp_splice" : "Expected (Splice Selection)",
        "canon" : "hg38 Reference",
        "exp" : "Expected",
        "ps" : "Proximally-Spliced",
        "ds" : "Distally-Spliced",
        "as" : "Alternatively-Spliced",
        "create" : "NAGNAG-Creating (Alternate Allele)",
        "create_alt" : "NAGNAG-Creating (Alternate Allele)",
        "alter_ref" : "NAGNAG-Altering (Reference Allele)",
        "alter_alt" : "NAGNAG-Altering (Alternate Allele)",
        "destroy" : "NAGNAG-Destroying (Reference Allele)",
        "destroy_ref" : "NAGNAG-Destroying (Reference Allele)",
        "dbSNP_common" : "dbSNP Common",
        "dbSNP_rare" : "dbSNP Rare",
        "ClinVar" : "ClinVar",
        "phase" : "Phase",
        "vc_ps" : "Proximal Variable Codons",
        "vc_ds" : "Distal Variable Codons",
        "aa_ps" : "Proximal Variable AAs",
        "aa_ds" : "Distal Variable AAs",
        "aatype" : "Amino Acid Transition Type",
        "val" : "Amino Acid",
        "count" : "Number of Scenarios",
        "prop" : "Proportion of Scenarios"
    }

    prop_colors = {
        "exp_stoch" : "#3B3838",
        "exp_splice" : "#01295F",
        "canon" : "#D25A46",
        "ps" : "#006034",
        "ds" : "#84007A",
        "as" : "#046198",
        "dbSNP_common" :  "#265E19",
        "dbSNP_rare" :  "#2E6366",
        "ClinVar" : "#FFB521",
        "create" : "#ED790C",
        "alter" : "#BD5B5B",
        "destroy" : "#2D5F91",
        "alter_ref" : "#A04040",
        "alter_alt" : "#D08A8A"
    }
