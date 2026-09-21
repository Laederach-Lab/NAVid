import os, pathlib, argparse, pickle, tarfile, itertools, pandas as pd, numpy as np

# Arguments
parser = argparse.ArgumentParser(prog="update-data.py", description="Update data based on analyses from nagnag-3ss (https://github.com/hallep/nagnag-3ss)")
parser.add_argument("in-dir", type=str, metavar="INPUT/DIR", help="path to nagnag-3ss directory")
parser.add_argument("out-dir", type=str, metavar="OUTPUT/DIR", help="path to NAVid directory")
args = parser.parse_args()

NAGNAG_DIR = pathlib.Path(args.in_dir)
if not NAGNAG_DIR.is_dir():
    print("Invalid input directory!")
    parser.error("Invalid input directory!")

NAVID_DIR = pathlib.Path(args.out_dir)
if not NAVID_DIR.is_dir():
    parser.error("Invalid output directory!")

# Directories
TMP_DIR = f"{NAGNAG_DIR}/app_data"
OUT_DIR = f"{NAVID_DIR}/nagnag/data"

os.makedirs(f"{NAGNAG_DIR}/app_data")
for subdir in ["dbs/ClinVar", "export", "obs/ref", "obs/vnt", "seq", "stats/tsns", "stats/vnts"]:
    os.makedirs(f"{OUT_DIR}/{subdir}", exist_ok=True)

# Utilities
chroms = [f"chr{c}" for c in list(map(str, range(1, 23))) + ["X", "Y", "M"]]

dtypes = {

    # 1off_splice_sites.txt
    "csite_starts" : "str",
    "csite_pos" : "str",
    "csite_inds" : "str",

    # vcf columns
    "chrom" : "str",
    "pos" : "int",
    "id" : "str",
    "ref" : "str",
    "alt" : "str",
    "qual" : "str",
    "filter" : "str",
    "info" : "str",
}

vnt_dbs = vnt_dbs_all = ["dbSNP_common", "dbSNP_rare", "ClinVar"]
if pathlib.Path(f"{NAGNAG_DIR}/variants/bed/HGMD_splice.bed").exists:
    vnt_dbs_all = ["dbSNP_common", "dbSNP_rare", "ClinVar", "HGMD_splice"]

canon_fmt = {c : x for x,c in enumerate(chroms)}

def chr_fmt(c:str) -> int:
    if c in canon_fmt:
        return canon_fmt[c]
    return 0

# NAG Strength (stats/nag_strength)
def nag_strength():
    strength = pd.read_csv(f"{NAGNAG_DIR}/stats/1nag_motif_freq.txt", sep="\t", index_col=0
                           )["prop_ss"].loc[["AAG", "CAG", "GAG", "TAG"]].to_dict()
    pickle.dump(strength, open(f"{OUT_DIR}/stats/nag_strength", "wb"))

# hg38 Reference Genome Sequences (hg38, seq/)
def sort_hg38():

    # Canonical Chromosomes (dictionary)
    hg38:dict = pickle.load(open(f"{NAGNAG_DIR}/src/hg38", "rb"))
    keys = list(hg38)
    [hg38.pop(c) for c in keys if c not in canon_fmt]
    pickle.dump(hg38, open(f"{OUT_DIR}/hg38", "wb"))

    # By Chromosome
    hg38:dict = pickle.load(open(f"{NAGNAG_DIR}/src/hg38", "rb"))

    # canonical chromosomes
    for chrom in list(canon_fmt):
        pickle.dump(hg38.pop(chrom), open(f"{OUT_DIR}/seq/{chrom}", "wb"))

    # non-canonical chromosomes
    pickle.dump(hg38, open(f"{OUT_DIR}/seq/non-canonical", "wb"))

# NAGNAG Instances (obs/ref/)
def canon_sites():
    
    ''' Canonical splice sites
    
    Columns: "chrom", "motif_start", "motif_end", "strand", "motif", "splice_type"
    '''

    cols = ["chrom", "ssite_start", "ssite_end", "strand", "ssite_seq", "csite_pos"]

    # DataFrame
    sites = pd.read_csv(f"{NAGNAG_DIR}/sites/nagnag_3ss.txt", sep="\t", usecols=cols)[cols]
    sites["csite_pos"] = sites["csite_pos"].apply(lambda x: "PS" if (x == "0") else "DS" if (x == "1") else "AS")
    sites.columns = ["chrom", "motif_start", "motif_end", "strand", "motif", "splice_type"]

    # Sort
    sites["x"] = sites["chrom"].apply(lambda x: chr_fmt(x)).astype(int)
    sites.sort_values(by=["x", "motif_start", "strand"], inplace=True)
    sites.drop(columns="x", inplace=True)

    # Save
    sites.to_csv(f"{TMP_DIR}/reference_sites.txt", sep="\t", index=False)
    pickle.dump(sites, open(f"{OUT_DIR}/obs/ref/reference_sites", "wb"))

def canon_scenarios():
    
    ''' Canonical splice scenarios

    Columns: "chrom", "motif_start", "motif_end", "strand", "phase", "motif", "splice_type", "up_end",
    "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aatt", "ins_aa", "del_aa"
    '''

    cols = ["chrom", "ssite_start", "ssite_end", "strand", "phase", "ssite_seq", "splice_type",
            "up_end", "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aattype", "ins_aa", "del_aa"]
    
    # DataFrame
    scens = pd.read_csv(f"{NAGNAG_DIR}/sites/nagnag_scens.txt", sep="\t", usecols=cols)[cols]
    scens.columns = ["chrom", "motif_start", "motif_end", "strand", "phase", "motif", "splice_type",
                     "up_end", "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aatt", "ins_aa", "del_aa"]

    # Sort
    scens["x"] = scens["chrom"].apply(lambda x: chr_fmt(x)).astype(int)
    scens.sort_values(by=["x", "motif_start", "strand", "phase"], inplace=True)
    scens.drop(columns="x", inplace=True)
    
    # Save
    scens.to_csv(f"{TMP_DIR}/reference_scenarios.txt", sep="\t", index=False)
    pickle.dump(scens, open(f"{OUT_DIR}/obs/ref/reference_scenarios", "wb"))

# Variant Instances (obs/vnt/)
def variant_events():

    ''' Variant Events IDs
    
    Columns: "id", "chrom", "pos", "strand", "db", "effect", "ref", "alt", "motif_start", "motif_end"
    '''

    cols = ["vnt_id", "chrom", "vnt_coord", "strand", "alt", "db", "vnt_effect", "ref_seq", "alt_seq", "motif_start", "motif_end"]

    def process_db(db:str):
        scen = pd.read_csv(f"{NAGNAG_DIR}/variants/found/{db}_nagnag_vnt_events.txt",
                           sep="\t", index_col=0, dtype=dtypes)
        scen["vnt_effect"] = scen["vnt_effect"].str.lower()
        scen["db"] = db

        # sort
        scen["x"] = scen["chrom"].apply(lambda x: chr_fmt(x)).astype(int)
        scen.sort_values(by=["x", "vnt_coord", "strand"], inplace=True)
        scen.drop(columns="x", inplace=True)

        return scen[cols]

    # Stack
    df = pd.concat((process_db(db) for db in vnt_dbs), axis=0, ignore_index=True)
    df.columns = ["id", "chrom", "pos", "strand", "allele", "db", "effect", "ref", "alt", "motif_start", "motif_end"]

    # Save
    df.to_csv(f"{TMP_DIR}/variant_events.txt", sep="\t", index=False)    
    pickle.dump(df, open(f"{OUT_DIR}/obs/vnt/variant_events", "wb"))

def variant_scenarios():

    ''' Variant scenario IDs 
    
    Columns: "id", "chrom", "pos", "strand", "phase", "db", "effect", "ref", "alt",
    "motif_start", "motif_end", "up_end", "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aatt", "ins_aa", "del_aa"
    '''

    # Columns
    info = ["vnt_id", "chrom", "vnt_coord", "strand", "alt", "phase", "db", "vnt_effect", "ref_seq", "alt_seq", "motif_start", "motif_end", "up_end"]
    prop = ["vc_ps", "vc_ds", "aa_ps", "aa_ds", "aattype", "ins_aa", "del_aa"]

    names = ["id", "chrom", "pos", "strand", "allele", "phase", "db", "effect", "ref", "alt",
             "motif_start", "motif_end", "upstream_end",
             "vc_ps", "vc_ds", "aa_ps", "aa_ds", "aatt", "ins_aa", "del_aa"]

    def process_db(db:str):

        scen = pd.read_csv(f"{NAGNAG_DIR}/variants/found/{db}_nagnag_vnt_scenarios.txt",
                           sep="\t", index_col=0, dtype=dtypes)
        scen["db"] = db
        
        def process_effect(effect:str, allele:str):

            e = scen[scen["vnt_effect"].str.lower() == effect][info + [f"{p}_{allele}" for p in prop]]
            e["vnt_effect"] = f"{effect}_{allele}"
            e.columns = names

            # sort
            e["x"] = e["chrom"].apply(lambda x: chr_fmt(x)).astype(int)
            e.sort_values(by=["x", "pos", "strand", "phase"], inplace=True)
            e.drop(columns="x", inplace=True)

            return e
        
        return pd.concat((process_effect("create", "alt"),
                          process_effect("alter", "ref"),
                          process_effect("alter", "alt"),
                          process_effect("destroy", "ref")
                          ), axis=0, ignore_index=True)

    # Stack
    df = pd.concat((process_db(db) for db in vnt_dbs), axis=0, ignore_index=True)

    # Save
    df.to_csv(f"{TMP_DIR}/variant_scenarios.txt", sep="\t", index=False)
    pickle.dump(df, open(f"{OUT_DIR}/obs/vnt/variant_scenarios", "wb"))

# Statistics (stats/tsns/, stats/vnts/)
def variant_frequencies():

    ''' Variant Frequencies

    Effect Frequencies
    * index: effect
    * columns: (freq, db)
        
    Variant Type Frequencies
    * index: (effect, ref, alt)
    * columns: (freq, db)
    '''

    effects = ["create", "alter", "destroy"]

    # Variant Events
    obs = [pd.read_csv(f"{NAGNAG_DIR}/variants/found/{db}_nagnag_vnt_events.txt", sep="\t", index_col=0, dtype=dtypes) for db in vnt_dbs]
    events = [list(zip(df["ref_seq"].str.upper(), df["alt_seq"].str.upper(), df["vnt_effect"].str.lower())) for df in obs]

    n_events = np.array([len(df) for df in events]).reshape(1, -1)

    # ===== Effect Frequencies ===== #

    # Counts
    count_types = np.array([[(df["vnt_effect"].str.lower() == e).sum() for df in obs] for e in effects])
    count_df = pd.DataFrame(count_types, index=effects, columns=vnt_dbs)

    # Proportions
    prop_df = count_df.copy()
    for n,c in zip(n_events[0], vnt_dbs):
        prop_df[c] = count_df[c] / n

    # Combine & Save
    type_df = pd.concat([count_df, prop_df], axis=1, keys=["count", "prop"])
    type_df.index.names = ["effect"]
    pickle.dump(type_df, open(f"{OUT_DIR}/stats/vnts/effect_frequencies", "wb"))

    # ===== Variant Event Frequencies ===== #

    # Motif Possibilities
    N = ["A", "C", "G", "T"]
    B = ["C", "G", "T"]
    H = ["A", "C", "T"]

    # canonical NAGNAGs
    canon = [f"{n1}AG{n2}AG" for n1,n2 in itertools.product(N, N)]

    # 1-off NAGNAGs
    off1 = [f"{n1}{b}G{n2}AG" for n1,b,n2 in itertools.product(N, B, N)]
    off2 = [f"{n1}A{h}{n2}AG" for n1,h,n2 in itertools.product(N, H, N)]
    off4 = [f"{n1}AG{n2}{b}G" for n1,n2,b in itertools.product(N, N, B)]
    off5 = [f"{n1}AG{n2}A{h}" for n1,n2,h in itertools.product(N, N, H)]
    off = off1 + off2 + off4 + off5

    # Variant Possibilities
    create_poss = [(r, a, "create") for r,a in itertools.product(off, canon)]
    alter_poss = [(r, a, "alter") for r,a in itertools.product(canon, canon)]
    destroy_poss = [(r, a, "destroy") for r,a in itertools.product(canon, off)]
    
    # Counts
    count_create = np.array([[df.count(p) for p in create_poss] for df in events]).T
    count_alter = np.array([[df.count(p) for p in alter_poss] for df in events]).T
    count_destroy = np.array([[df.count(p) for p in destroy_poss] for df in events]).T

    # Observed Variant Events
    create = {(r,a) : c for (r,a,_),c in zip(create_poss, count_create) if sum(c) > 0}
    alter = {(r,a) : c for (r,a,_),c in zip(alter_poss, count_alter) if sum(c) > 0}
    destroy = {(r,a) : c for (r,a,_),c in zip(destroy_poss, count_destroy) if sum(c) > 0}

    # DataFrames
    create_df = pd.DataFrame(data=create.values(), index=create.keys(), columns=vnt_dbs)
    alter_df = pd.DataFrame(data=alter.values(), index=alter.keys(), columns=vnt_dbs)
    destroy_df = pd.DataFrame(data=destroy.values(), index=destroy.keys(), columns=vnt_dbs)

    # Combine SNP Databases
    count = pd.concat([create_df, alter_df, destroy_df], keys=["create", "alter", "destroy"])
    count.index.names = ["effect", "ref", "alt"]

    # Proportions
    prop = count.copy()
    for n,c in zip(n_events[0], vnt_dbs):
        prop[c] = count[c] / n

    # Combine & Save
    event_df = pd.concat([count, prop], axis=1, keys=["count", "prop"])
    pickle.dump(event_df, open(f"{OUT_DIR}/stats/vnts/event_frequencies", "wb"))

def transition_frequencies():

    ''' Get transition frequencies 
    
    Variable Codons
    * index: (phase, vc_ps, vc_ds)
    * columns: (freq, effect, database)
    
    Variable Amino Acids
    * index: (phase, aa_ps, aa_ds)
    * columns: (freq, effect, database)
    
    Transition Types
    * index: (phase, aattype)
    * columns: (freq, effect, database)
    
    Outcomes
    * index: (cat, phase, val)
    * columns: (freq, effect, database)
    '''

    # variant effects
    effects = {
        "create" : "create_alt",
        "alter_ref" : "alter_ref",
        "alter_alt" : "alter_alt",
        "destroy" : "destroy_ref"
    }

    # columns to rename
    rename = {
        "ref" : ("ref", "canon"),
        "ref_ps" : ("ref", "ps"),
        "ref_ds" : ("ref", "ds"),
        "ref_as" : ("ref", "as"),
        "exp_stoch" : ("ref", "exp_stoch"),
        "exp_splice" : ("ref", "exp_splice")
    }
    rename.update({f"{db}_{ef}" : (effects[ef], db) for db,ef in itertools.product(vnt_dbs_all, list(effects))})

    def fmt(i:str, j:str) -> tuple[str, str, str]:
        x,y = rename[j]
        return (i, x, y)

    # columns to drop
    to_drop = list(itertools.product(["count", "prop"], vnt_dbs_all))
    to_drop += list(itertools.product(["prop"], ["exp_ps", "exp_ds", "exp_as"]))

    # column order
    order = list(itertools.product(["count"], ["ref"], ["canon", "ps", "ds", "as"]))
    order += list(itertools.product(["count"], effects.values(), vnt_dbs_all))
    order += list(itertools.product(["prop"], ["ref"], ["exp_stoch", "exp_splice", "canon", "ps", "ds", "as"]))
    order += list(itertools.product(["prop"], effects.values(), vnt_dbs_all))

    # ===== Variable Codons ===== #

    # Frequencies
    vcC = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/variable_codons_freq_count.txt", sep="\t", index_col=[0, 1, 2], header=0)
    vcP = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/variable_codons_freq.txt", sep="\t", index_col=[0, 1, 2], header=0)
    vc_df = pd.concat([vcC, vcP], axis=1, keys=["count", "prop"])

    vc_df.drop(columns=to_drop, inplace=True)
    vc_df.drop(columns=list(itertools.product(["count", "prop"], ["aa_ps", "aa_ds", "aattype", "ins", "del"])), inplace=True)
    vc_df.drop(labels="HGMD_splice", level=2, axis=1, inplace=True, errors="ignore")

    vc_df.columns = pd.MultiIndex.from_tuples([fmt(i, j) for i,j in vc_df.columns.values])
    vc_df = vc_df[order]
        
    # Save
    pickle.dump(vc_df, open(f"{OUT_DIR}/stats/tsns/variable_codon_frequencies", "wb"))

    # ===== Variable Amino Acids ===== #

    # Frequencies
    aaC = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/variable_amino_acids_freq_count.txt", sep="\t", index_col=[0, 1, 2], header=0)
    aaP = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/variable_amino_acids_freq.txt", sep="\t", index_col=[0, 1, 2], header=0)
    aa_df = pd.concat([aaC, aaP], axis=1, keys=["count", "prop"])

    aa_df.drop(columns=to_drop, inplace=True)
    aa_df.drop(columns=list(itertools.product(["count", "prop"], ["num_vc", "prop_vc", "aattype", "ins", "del"])), inplace=True)
    aa_df.drop(labels="HGMD_splice", level=2, axis=1, inplace=True, errors="ignore")

    aa_df.columns = pd.MultiIndex.from_tuples([fmt(i, j) for i,j in aa_df.columns.values])
    aa_df = aa_df[order]

    # Save
    pickle.dump(aa_df, open(f"{OUT_DIR}/stats/tsns/variable_amino_acid_frequencies", "wb"))

    # ===== Amino Acid Transition Types ===== #

    # Frequencies
    ttC = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/amino_acid_transitions_freq_count.txt", sep="\t", index_col=[0, 1], header=0)
    ttP = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/amino_acid_transitions_freq.txt", sep="\t", index_col=[0, 1], header=0)
    tt_df = pd.concat([ttC, ttP], axis=1, keys=["count", "prop"])

    tt_df.drop(columns=to_drop, inplace=True)
    tt_df.drop(columns=list(itertools.product(["count", "prop"], ["num_vc", "prop_vc", "num_aat", "prop_aat"])), inplace=True)
    tt_df.drop(labels="HGMD_splice", level=2, axis=1, inplace=True, errors="ignore")

    tt_df.columns = pd.MultiIndex.from_tuples([fmt(i, j) for i,j in tt_df.columns.values])
    tt_df = tt_df[order]

    # Save
    pickle.dump(tt_df, open(f"{OUT_DIR}/stats/tsns/transition_type_frequencies", "wb"))

    # ===== Amino Acid Outcomes ===== #

    # Frequencies
    outC = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/amino_acid_outcomes_freq_count.txt", sep="\t", index_col=[0, 1, 2], header=0)
    outP = pd.read_csv(f"{NAGNAG_DIR}/variants/proteome/amino_acid_outcomes_freq.txt", sep="\t", index_col=[0, 1, 2], header=0)
    out_df = pd.concat([outC, outP], axis=1, keys=["count", "prop"]).reorder_levels(order=["cat", "phase", "aa"], axis="index")

    out_df.drop(columns=to_drop, inplace=True)
    out_df.drop(labels="HGMD_splice", level=2, axis=1, inplace=True, errors="ignore")

    out_df.columns = pd.MultiIndex.from_tuples([fmt(i, j) for i,j in out_df.columns.values])
    out_df = out_df[order]

    # Save
    pickle.dump(out_df, open(f"{OUT_DIR}/stats/tsns/aa_outcome_frequencies", "wb"))

# Compressed Data (export/)
def compress_data():

    to_compress = {
        "variant_effects" : ["event_frequencies", "effect_frequencies"],
        "proteome_effects" : ["variable_codons", "variable_amino_acids", "transition_types", "aa_outcomes"]
    }

    # compress frequencies
    for setname, filenames in to_compress.items():
        with tarfile.open(f"{OUT_DIR}/export/{setname}.tar.gz", "w:gz") as tar:
            for fn in filenames:
                tar.add(name=f"{TMP_DIR}/{fn}.txt", arcname=f"{fn}.txt")

    # compress observed
    with tarfile.open(f"{OUT_DIR}/export/observed_nagnags.tar.gz", "w:gz") as tar:
        for fn in ["reference_sites", "reference_scenarios", "variant_events", "variant_scenarios"]:
            tar.add(name=f"{TMP_DIR}/{fn}.txt", arcname=f"{fn}.txt")

# === RUN === #

nag_strength()

sort_hg38()

canon_sites()
canon_scenarios()

variant_events()
variant_scenarios()

variant_frequencies()
transition_frequencies()

compress_data()
