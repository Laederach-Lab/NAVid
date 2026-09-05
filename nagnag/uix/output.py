from .format import Fmt
from .interface import Interface
from ..lib import ROOT, pickle, product, pd, np, plt, patches, ui, render
from ..utils import Seq, Scenario, Variant

class Output:
    
    ''' Formatted app output (tables, plots, etc.)  '''

    def __init__(self):
        pass
    
    # ===== Uploads ===== #

    @staticmethod
    def processed_upload(filepath:str) -> pd.DataFrame:
        data = pd.read_csv(filepath, sep="\t", header=None)

        data["effect"] = ["none"] * len(data)
        data["ref_motif"] = ["."] * len(data)
        data["alt_motif"] = ["."] * len(data)
        
        # for each variant:
        for i,vnt in data.iterrows():

            # validate data
            chrom = Seq.format_chrom(vnt[0])
            data.at[i, 0] = chrom
            
            try:
                if Seq.valid_coord(int(vnt[1])):
                    pos = int(vnt[1])
            except:
                data.at[i, "effect"] = "invalid"
                continue
            else:
                data.at[i, 1] = pos

            strand = vnt[2] if Seq.valid_strand(vnt[2]) else "+"
            data.at[i, 2] = strand
            
            alt = Seq.format_base(vnt[3])
            data.at[i, 3] = alt

            result = Variant.affected_nagnag(chrom, pos, strand, alt)

            if result == False:
                data.at[i, "effect"] = "invalid"
            
            elif result != None:
                data.at[i, "alt_motif"] = result[3]
                data.at[i, "ref_motif"] = result[4]
                data.at[i, "effect"] = Variant.nagnag_effect(result[3], result[4])

        return data

    # ===== Data ===== #

    # ----- Variants ----- #

    @staticmethod
    def vnt_events() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/obs/vnt/variant_events", "rb"))

    @staticmethod
    def vnt_scens() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/obs/vnt/variant_scenarios", "rb"))
    
    @staticmethod
    def ref_sites() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/obs/ref/reference_sites", "rb"))
    
    @staticmethod
    def ref_scens() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/obs/ref/reference_scenarios", "rb"))
    
    @staticmethod
    def event_freqs() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/stats/vnts/event_frequencies", "rb"))
    
    @staticmethod
    def effect_freqs() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/stats/vnts/effect_frequencies", "rb"))

    # ----- Transitions ----- #
    
    @staticmethod
    def vc() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/stats/tsns/variable_codon_frequencies", "rb"))
    
    @staticmethod
    def aat() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/stats/tsns/variable_amino_acid_frequencies", "rb"))
    
    @staticmethod
    def aatt() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/stats/tsns/transition_type_frequencies", "rb"))
    
    @staticmethod
    def out() -> pd.DataFrame:
        return pickle.load(open(f"{ROOT}/data/stats/tsns/aa_outcome_frequencies", "rb"))

    # ===== Tables ===== #

    @staticmethod
    def get_df(metric:str, phase:int=None, freq:str=None) -> pd.DataFrame:
        
        # Metric
        match metric.split("_")[0]:
            case "vc":
                df = Output.vc()
            case "aat":
                df = Output.aat()
            case "aatt":
                df = Output.aatt()
            case "ins":
                df = Output.out().loc["ins"]
            case "del":
                df = Output.out().loc["del"]
            case _:
                return pd.DataFrame()

        # Phase
        if phase in Fmt.metrics(metric, phases=True, inds=True):
            df = df.loc[phase]

        # Frequency
        if freq in ["prop", "count"]:
            df = df[freq]

        return df
        
    @staticmethod
    def scen_freq_table(vnt:Variant, scen:Scenario, metric:str, fmt:bool) -> pd.DataFrame:

        # Empty
        if vnt.empty or vnt.not_nagnag or (not scen.canon):
            return pd.DataFrame()
        
        # Get source DataFrame
        df = Output.get_df(metric.split("_")[0], scen.phase)
        if df.empty:
            return pd.DataFrame()

        # Get this scenario
        if metric == "vc":
            df = df.loc[(scen.vc_ps, scen.vc_ds)]
        elif metric == "aat":
            df = df.loc[(scen.aa_ps, scen.aa_ds)]
        elif metric == "aatt":
            df = df.loc[(scen.aatt)]
        elif metric.startswith("ins_") or metric.startswith("del_"):
            a = Output.cat_aa(metric, scen)
            if not a:
                return pd.DataFrame()
            df = df.loc[a]
        
        d = [(f"{vnt.effect}_{scen.allele}", db) for db in ["dbSNP_common", "dbSNP_rare", "ClinVar"]]
        n = [("ref", "canon")]
        p = [("ref", "exp_stoch"), ("ref", "exp_splice"), ("ref", "canon")]
        
        # Formatted DataFrame
        return pd.DataFrame({
            "Database" : Fmt.dbs("v" if fmt else "k"),
            "Number of Scenarios" : ["-", "-"] + df["count"][n+d].apply(lambda x: f"{int(x):,}").to_list(),
            "Proportion of Scenarios" : df["prop"][p+d].apply(lambda x: f"{x:.3%}").to_list()
        })

    # ===== Frequency Dropdown ===== #
    
    @staticmethod
    def vnt_ids_table(id:str, db:str, vnt:Variant, scen:Scenario, fmt:bool=False) -> pd.DataFrame:

        # Affected NAGNAG
        if id in ["nagnag_effect", "motif"]:
            data = Output.vnt_events()
            data = data[(data["db"] == db) & (data["effect"] == vnt.effect)]

            if id == "motif":
                data = data[(data["ref"] == vnt.ref_motif) & (data["alt"] == vnt.alt_motif)]

        # Protome Effects
        elif id.startswith(("vc_", "aat_", "aatt_", "ins_", "del_")):
            data = Output.vnt_scens()
            data = data[(data["db"] == db) & (data["effect"] == f"{vnt.effect}_{scen.allele}") &
                        (data["phase"] == scen.phase)]

            if id.startswith(("ins_", "del_")):
                a = Output.cat_aa(id, scen)
                if not a:
                    return pd.DataFrame()

            if id.startswith("vc_"):
                data = data[(data["vc_ps"] == scen.vc_ps) & (data["vc_ds"] == scen.vc_ds)]
            elif id.startswith("aat_"):
                data = data[(data["aa_ps"] == scen.aa_ps) & (data["vc_ds"] == scen.vc_ds)]
            elif id.startswith("aatt_"):
                data = data[data["aatt"] == scen.aatt]
            elif id.startswith("ins_"):
                data = data[data["ins_aa"].str.contains(a)]
            elif id.startswith("del_"):
                data = data[data["del_aa"].str.contains(a)]

        else:
            return pd.DataFrame()

        # Format
        if fmt:
            data = data[["id", "chrom", "pos", "strand", "ref", "alt"]]
            data = data.drop_duplicates(inplace=False)
            data.columns = ["Variant ID", "Chromosome", "Position", "Strand", "Reference", "Alternate"]

        return data

    # ----- Affected NAGNAG ----- #

    @staticmethod
    def nagnag_effect_plot(vnt:Variant) -> plt.Figure:
        
        effect = str(vnt.effect).lower()

        # Figure
        fig = plt.figure()
        ax = fig.add_subplot()
        ax.spines[["top", "right"]].set_visible(False)
        
        # Bar width/positioning
        width = 0.85
        p = 0.0125
        w = (width - (3 * p)) / 4

        x = np.arange(3)
        offset = [-((w+p)*3/2), -((w+p)/2), ((w+p)/2), ((w+p)*3/2)]

        vnt_dbs = {k:v for k,v in Fmt.dbs(vnt=True)}

        # Series
        dbs = {k : (v, Fmt.prop_colors[k]) for k,v in vnt_dbs.items()}

        # Bar Design
        alpha = 0.5
        match str(effect).lower():
            case "create":
                c, a, d = 1, alpha, alpha
            case "alter":
                c, a, d = alpha, 1, alpha
            case "destroy":
                c, a, d = alpha, alpha, 1
            case _:
                c, a, d = 1, 1, 1

        data = Output.effect_freqs()["prop"]

        # Bars
        for o,(k,(_,h)) in zip(offset, dbs.items()):
            ax.bar(x+o, data[k], width=w, color=[(h,c), (h,a), (h,d)])

        # Labels
        ax.set_xticks(range(3), [t.replace("-","-\n") for t in Fmt.effs("v")])
        if effect in data.index.values:
            ax.get_xticklabels()[data.index.to_list().index(effect)].set_fontweight("bold")
        ax.set_ylabel("Proportion of Variants", weight="bold")

        # Legend
        b = [ax.bar([0], [0], color=h, label=l) for l,h in dbs.values()]
        ax.legend(handles=b)

        return fig

    @staticmethod
    def motif_table(vnt:Variant) -> pd.DataFrame:

        # Empty
        if vnt.empty or vnt.not_nagnag:
            return None
        
        eff = str(vnt.effect).lower()
        ref = vnt.ref_motif
        alt = vnt.alt_motif

        df = Output.event_freqs().loc[(eff, ref, alt)]
        
        # No other occurances
        if sum(df.values) == 0:
            return None
        
        c = df.loc["count"].apply(lambda x: f"{int(x):,}").to_frame()
        p = df.loc["prop"].apply(lambda x: f"{x:.3%}").to_frame()

        df = pd.concat((c, p), axis=1).reset_index(names="Database")
        df["Database"] = [Fmt.dbs(d) for d in df["Database"]]
        df.columns = ["Database", "Number of Variants", "Proportion of Variants"]
        df.drop_duplicates(inplace=True)
        df.sort_values(by="Number of Variants", ascending=False, inplace=True)

        return df

    # ----- Proteome Effects ----- #

    @staticmethod
    def splice_junction_plot(vnt:Variant) -> plt.Figure:
        
        # Figure
        fig, ax = plt.subplots(figsize=(12.5, 4), layout="constrained")

        if vnt.empty or vnt.not_nagnag or vnt.not_scenario:
            return fig
        
        # Annotations
        cs5 = 0.45
        pcs = 1.4
        dcs = 1.55

        cU = "#445F44"
        cD = "#9B5B08"
        cN = "#A5240D"

        def plot_junction(y:int, nagnag:str, allele:str):
            
            # Splice Type
            if (nagnag[1:3] == "AG") and (nagnag[4:6] == "AG"):
                stype = "AS"
                nagnag_txt = nagnag[:3].lower() + nagnag[3:].upper()
            elif nagnag[1:3] == "AG":
                stype = "PS"
                nagnag_txt = nagnag[:3].lower() + nagnag[3:].upper()
            elif nagnag[4:6] == "AG":
                stype = "DS"
                nagnag_txt = nagnag.lower()

            # Exons
            eMin, eMax = y - 0.065, y + 0.075

            ax.vlines(cs5, eMin, eMax, colors=cU, lw=2, zorder=2)
            ax.add_patch(patches.Rectangle((-0.025, eMin), width=0.475, height=0.14, fc=(cU, 0.7), ec=cU, lw=0, zorder=1))

            if stype == "AS":
                ax.add_patch(patches.Rectangle((pcs, eMin), width=0.15, height=0.14, fc=(cN, 0.6), ec=cN, lw=0, zorder=1))
                ax.vlines([pcs, dcs], eMin, eMax, colors=[cN, cD], lw=2, zorder=2)
            else:
                ax.vlines(pcs if (stype == "PS") else dcs, eMin, eMax, colors=cD, lw=2, zorder=2)

            if stype == "PS":
                ax.add_patch(patches.Rectangle((pcs, eMin), width=0.625, height=0.14, fc=(cD, 0.65), ec=cD, lw=0, zorder=1))
            else:
                ax.add_patch(patches.Rectangle((dcs, eMin), width=0.475, height=0.14, fc=(cD, 0.65), ec=cD, lw=0, zorder=1))

            # Intron
            iMin, iMax = y - 0.02, y + 0.02

            ax.hlines(y, xmin=0.625, xmax=1.075, color="black", lw=3)
            ax.plot([0.825, 0.85], [iMin, iMax], color="black")
            ax.plot([0.85, 0.875], [iMin, iMax], color="black")

            # Splice Lines
            s = y / 0.15
            def draw_splice(c3:float, dy:float=0.1):
                y0 = (y + 0.085) if (s == 1) else (y - 0.075)
                y1 = y + ((dy + 0.085) * s)
                ax.plot([cs5, (cs5 + c3) / 2, c3], [y0, y1, y0], color="black")
              
            if stype == "PS":
                draw_splice(pcs)
            elif stype == "DS":
                draw_splice(dcs)
            elif stype == "AS":
                draw_splice(pcs)
                draw_splice(dcs, dy=0.175)

            # Description
            ax.text(-0.025, y+0.085, f"{allele.capitalize()} Allele", ha="left", va="bottom", size=15)

            # Bases
            # upstream exon
            up = "..." + vnt.eflank_5ss[-6:].upper() + vnt.iflank_5ss[:3].lower()
            for x,b in zip(np.linspace(0.025, 0.725, 15), up):
                ax.text(x, y, b, ha="center", va="center", family="monospace", size=23)

            # upstream intron
            intron = vnt.iflank_3ss[-9:-6].lower()
            for x,b in zip(np.linspace(1.125, 1.225, 3), intron):
                ax.text(x, y, b, ha="center", va="center", family="monospace", size=23)

            # nagnag
            for x,b in zip(np.linspace(1.275, 1.525, 6), nagnag_txt):
                ax.text(x, y, b, ha="center", va="center", family="monospace", size=23, weight="bold")

            # downstream exon
            down = vnt.eflank_3ss[:6].upper() + "..."
            for x,b in zip(np.linspace(1.575, 1.975, 9), down):
                ax.text(x, y, b, ha="center", va="center", family="monospace", size=23)

        # reference + alternative alleles
        plot_junction(0.15, vnt.ref_motif.upper(), "Reference")
        plot_junction(-0.15, vnt.alt_motif.upper(), "Alternate")

        # Variant
        xv = np.linspace(1.275, 1.525, 6)[vnt.var_pos]
        ax.hlines([-0.185, 0.115], xmin=xv-0.0165, xmax=xv+0.02, colors="black")
        sd = patches.ArrowStyle("->", head_length=0.7, head_width=0.35)
        ax.annotate("", xy=(xv, -0.075), xytext=(xv, 0.075), arrowprops=dict(arrowstyle=sd, color="black", lw=2))
        
        # Axes
        ax.margins(0)
        ax.set_axis_off()
        ax.set_xlim(-0.05, 2.05)
        ax.set_ylim(-0.45, 0.45)

        return fig

    @staticmethod
    def cat_aa(metric:str, scen:Scenario) -> str:
        
        if (metric[:3] != "ins") and (metric[:3] != "del"):
            return
        
        aa = scen.ins_aa if (metric[:3] == "ins") else scen.del_aa
        if (aa == None) or (aa == ["_"]) or (aa == ["."]):
            return

        try:
            i = int(metric.split("_")[1][1:])
            return aa[i]
        except ValueError:
            return
        except IndexError:
            return

    @staticmethod
    def tsn_freq_desc(vnt:Variant, scen:Scenario, metric:str):
        
        # Empty
        if vnt.empty or vnt.not_nagnag:
            return

        # Insertion/Deletion
        if metric in ["ins", "del"]:
            return None
        
        # Description
        if metric == "vc":
            this = f"({scen.vc_ps}, {scen.vc_ds})"
        elif metric == "aat":
            this = f"({scen.aa_ps}, {scen.aa_ds})"
        elif metric == "aatt":
            this = Seq.aa[scen.aatt].codes if (scen.phase == 0) else Seq.tsn[scen.aatt]
        elif metric.startswith("ins"):
            a = Output.cat_aa(metric, scen)
            if not a:
                return
            this = f"{Seq.aa[a].name} insertion"
        elif metric.startswith("del"):
            a = Output.cat_aa(metric, scen)
            if not a:
                return
            this = f"{Seq.aa[a].name} deletion"
        
        data = Output.scen_freq_table(vnt, scen, metric, False)
        if data.empty:
            return None
        
        data = data.set_index("Database")["Proportion of Scenarios"]
        df = data.iloc[3:].apply(lambda x: float(x[:-1]) / 100).sort_values(ascending=False, inplace=False)
        
        c = float(data.loc["canon"][:-1]) / 100
        
        # No other occurances
        if sum(df.values) == 0:
            return Interface.freq(
                f"{this} is not listed in dbSNP or ClinVar",
                f"Canonically, it occurs in {c:.2%} of phase {scen.phase} splice scenarios"
            )
        
        v = v = df.values[0]

        # comparison to canonical
        n = c * 0.05
        if (v > (c-n)) and (v < (c+n)):
            compare = "about as frequent as"
        elif v < c:
            compare = "less frequent than"
        elif v > c:
            compare = "more frequent than"
        
        return Interface.freq(
            f"{this} is most common in {Fmt.col_names[df.index.values[0]]}, representing {v:.2%} of splice scenarios",
            f"This is {compare} what occurs canonically ({c:.2%})"
        )

    @staticmethod
    def rendered_tsn_freq_table(vnt:Variant, scen:Scenario, metric:str):
        
        # Empty
        if vnt.empty or vnt.not_nagnag:
            return None

        # Navset for ins/del
        if metric in ["ins", "del"]:

            # Get amino acids
            aa = scen.ins_aa if (metric == "ins") else scen.del_aa

            # No inserted/deleted amino acid(s)
            if (aa == None) or (aa == ["_"]) or (aa == ["."]):
                return None
            
            return ui.navset_pill(
                *[
                    ui.nav_panel(
                        Seq.aa[a].codes,
                        Interface.freq_div(
                            desc_id = f"{metric}_a{i}_p{scen.phase}_freq_desc",
                            display_id = f"{metric}_a{i}_p{scen.phase}_freq_table",
                            button_id = f"{metric}_a{i}_p{scen.phase}_obs_vnt_ids",
                            display = ui.output_data_frame
                        ),
                    ) for i,a in enumerate(aa)
                ],
            )
        
        data = Output.scen_freq_table(vnt, scen, metric, True)
        if data.empty:
            return None

        return render.DataTable(
            data = data,
            height = None,
            styles = [
                {"class" : "table-values numeric"},
                {"cols" : 0, "class" : "table-index"}
            ],
        )

    # ===== Data Tabs ===== #

    @staticmethod
    def variant_effects_table(eff:str, freq:str, fmt:bool) -> pd.DataFrame:

        # Data
        if eff == "effects_summary":
            df = Output.effect_freqs()[freq]
            df.index = Fmt.effs("v")
        else:
            df = Output.event_freqs().loc[eff.split("_")[0]][freq]

        df.columns = Fmt.dbs("v", vnt=True)

        # Format Numbers
        if fmt:
            for c in df.columns:
                if freq == "prop":
                    df[c] = df[c].apply(lambda x: f"{x:.3%}")
                else:
                    df[c] = df[c].apply(lambda x: f"{int(x):,d}")

        # Indices
        if eff == "effects_summary":
            df.insert(0, "Effect", df.index)
        else:
            df.insert(0, "Reference", df.index.get_level_values(0))
            df.insert(1, "Alternate", df.index.get_level_values(1))

        return df

    @staticmethod
    def tsn_table(metric:str, phase:int=None, effects:str|list[str]=None, dbs:list=None, freq:str="prop", fmt:bool=True) -> pd.DataFrame:
        
        # Metric, Phase, Frequency
        df = Output.get_df(metric, phase, freq)

        # Effects
        if effects == None:
            effects = Fmt.scen_effs("k")
        
        effects = [e for e in effects if (e in Fmt.scen_effs("k"))]

        # Databases
        if dbs == None:
            dbs = Fmt.dbs("k")
        valid = Fmt.dbs("k", prop=(freq == "prop"))
        dbs = [d for d in dbs if d in valid]

        # Columns
        ref = ["exp_stoch", "exp_splice", "canon"]
        canon = [d for d in dbs if d in ref]
        variants = [d for d in dbs if d not in ref]
        cols = [("ref", d) for d in canon] + list(product(effects, variants))
        
        if len(cols) == 0:
            return pd.DataFrame()
        
        # Sort
        if metric == "aatt":
            if phase == 0:
                df = df.loc[["E", "Q", "K", "*"]]
            else:
                df = df.loc[["DID", "NID", "CID", "IDR", "NC", "ET"]]
        elif freq == "prop":
            df = df.sort_values(by=[("ref", "exp_splice"), ("ref", "canon")], axis=0, ascending=False)
        elif freq == "prop":
            df = df.sort_values(by=("ref", "canon"), axis=0, ascending=False)

        df = df[cols]

        # Remove null values
        if metric in ["ins", "del"]:
            df = df[df[cols[0]].astype(str).apply(lambda x: float(x.removesuffix("%"))) > 0]

        # Formatting
        if fmt:

            # Data values
            for c in df.columns:
                if freq == "prop":
                    df[c] = df[c].apply(lambda x: f"{x:.3%}")
                else:
                    df[c] = df[c].apply(lambda x: f"{int(x):,d}")

            # Column names
            col_names = [Fmt.col_names[d] for d in canon]
            if (len(effects) == 1) and (len(variants) == 1):
                col_names += [Fmt.col_names[freq]]
            elif len(effects) == 1:
                col_names += [Fmt.col_names[v] for v in variants]
            elif len(variants) == 1:
                col_names += [Fmt.col_names[e] for e in effects]
            elif len(effects) > 1:
                col_names += [f"{Fmt.col_names[e]}\n{Fmt.col_names[v]}" for e,v in product(effects, variants)]
            df.columns = col_names        

            # Index names
            ind_names = {
                "vc": ["Proximal Codons", "Distal Codons"],
                "aat" : ["Proximal Residues", "Distal Residues"],
                "aatt" : ["Transition Type"],
                "ins" : ["Amino Acid"],
                "del" : ["Amino Acid"],
            }[metric]
            [df.insert(i, n, df.index.get_level_values(i)) for i,n in enumerate(ind_names)]

            # Index Values
            if metric == "aatt":
                df["Transition Type"] = df["Transition Type"].apply(lambda x: Seq.tsn[x])
            if metric in ["ins", "del"]:
                df["Amino Acid"] = df["Amino Acid"].apply(lambda x: Seq.aa[x].full_name)
            
        return df

    @staticmethod
    def obs_sites_table(db:str, effects:str|list[str]=None, fmt:bool=True) -> pd.DataFrame:

        # Dataframe
        if db == "canon":
            df = Output.ref_sites()

            # Format values
            if fmt:
                st = {
                    "PS" : "Proximal",
                    "DS" : "Distal",
                    "AS" : "Alternative",
                }
                df["splice_type"] = df["splice_type"].apply(lambda x: st[x])

                df.columns = ["Chromosome", "Site Start", "Site End", "Strand", "Motif", "Splice Type"]

        else:
            df = Output.vnt_events()

            # Isolate database
            df = df[df["db"] == db]
            df.drop(columns="db", inplace=True)

            # Filter effects
            if effects != None:
                if isinstance(effects, str):
                    effects = [effects]
                df = df[df["effect"].apply(lambda x: x in effects)]

            # Format values
            if fmt:
                df["effect"] = df["effect"].apply(lambda x: Fmt.effs(x, abbr=True))
                df.columns = ["Variant ID", "Chromosome", "Position", "Strand", "Variant",
                              "Effect", "Reference", "Alternate", "Motif Start", "Motif End"]

        return df

    @staticmethod
    def obs_scens_table(db:str, phase:int=None, effects:str|list[str]=None, fmt:bool=True) -> pd.DataFrame:

        # Dataframe
        if db == "canon":
            df = Output.ref_scens()
            cols = ["Chromosome", "Site Start", "Site End", "Strand", "Motif", "Splice Type"]

            # Format values
            if fmt:
                st = {
                    "PS" : "Proximal",
                    "DS" : "Distal",
                    "AS" : "Alternative",
                }
                df["splice_type"] = df["splice_type"].apply(lambda x: st[x])

        else:
            df = Output.vnt_scens()
            cols = ["Variant ID", "Chromosome", "Position", "Strand", "Variant", "Effect",
                    "Reference", "Alternate", "Motif Start", "Motif End"]

            # Isolate database
            df = df[df["db"] == db]
            df.drop(columns="db", inplace=True)

            # Filter effects
            if effects != None:
                if isinstance(effects, str):
                    effects = [effects]
                df = df[df["effect"].apply(lambda x: x in effects)]
                
            # Format values
            if fmt:
                df["effect"] = df["effect"].apply(lambda x: Fmt.scen_effs(x, abbr=True))

        # Filter by phase
        if phase == None:
            cols.insert(4, "Phase")
        else:
            df = df[df["phase"] == phase]
            df.drop(columns="phase", inplace=True)

        # Columns
        if fmt:
            if phase == -1:
                cols.extend(["5' Splice Site"])
                df.drop(columns=["vc_ps", "vc_ds", "aa_ps", "aa_ds", "aatt", "ins_aa", "del_aa"], inplace=True)
            else:
                cols.extend(["5' Splice Site", "Proximal Codons", "Distal Codons", "Proximal Residues", "Distal Residues",
                             "Transition", "Insertions", "Deletions"])
            df.columns = cols

        return df
