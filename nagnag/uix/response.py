from .format import Fmt

class Response:
    
    ''' App response functions

    Can be executed using the exec() function
    '''

    def __init__(self):
        pass

    # ========== SEARCH TAB ========== #

    # === Input === #

    @staticmethod
    def reactive_input_text(var:str) -> str:
        return f"""
@reactive.effect
@reactive.event(input.{var})
def _():
    v = state.get().update_{var}(input.{var}())
    
    # Invalid input
    if v.startswith("ERROR"):
        ui.update_text(id="{var}", value="")
        
        # error message
        ui.modal_show(Interface.error_modal(v))

    # Valid input
    else:
        ui.update_text(id="{var}", value=v)
        """

    @staticmethod
    def reactive_input_radio(var:str) -> str:
        return f"""
@reactive.effect
@reactive.event(input.{var})
def _():
    state.get().update_{var}(input.{var}())
        """
    
    @staticmethod
    def reactive_select_allele(alt:str) -> str:
        return f"""
@reactive.effect
@reactive.event(input.select_alt{alt})
def _():
    state.get().update_alt_allele("{alt}")
    ui.update_text(id="alt_allele", value="{alt}")
    ui.modal_remove()
        """

    # === Results === #

    @staticmethod
    def reactive_render_vnt_ids(id:str) -> str:

        def fn() -> str:
            if id == "nagnag_effect": return "nagnag_effect"
            if id == "motif": return "motif_effect"
            if id.startswith("vc_"): return "variable_codons"
            if id.startswith("aat_"): return "amino_acids"
            if id.startswith("aatt_"): return "transition_type"
            if id.startswith("ins_"): return "same_inserted_aa"
            if id.startswith("del_"): return "same_deleted_aa"
            return ""

        def table(db:str) -> str:
            return f"""
@render.data_frame
def {id}_{db}_vnts_table():
    data = Output.vnt_ids_table(id="{id}", db="{db}", vnt=vnt.get(), scen=scen.get(), fmt=True)
    if data.empty:
        return None
    
    return render.DataTable(
        data = data,
        width = "100%",
        height = "45dvh",
        styles = [
            {{"class" : "table-values"}},
            {{"cols" : 0, "class" : "table-index"}}
        ],
    )

@render.download(filename="{db}_variants_same_{fn()}.txt")
def {id}_{db}_vnts_download():
    data = Output.vnt_ids_table(id="{id}", db="{db}", vnt=vnt.get(), scen=scen.get(), fmt=False)
    if data.empty:
        return None
    
    yield data.to_csv(sep="\t", header=True, index=False)
            """
        
        return f"""
@reactive.effect
@reactive.event(input.{id}_obs_vnt_ids)
def _():
    return ui.modal_show(Interface.vnt_ids_modal("{id}"))

{"\n".join([table(db) for db in Fmt.dbs("k", vnt=True)])}
        """
    
    # ----- Proteome Effects ----- #

    @staticmethod
    def render_tsn_result(phase:str, metric:str) -> str:
        def val() -> str:
            if metric in ["ins", "del"]:
                return f"', '.join([Seq.aa[a].codes for a in scen.get().{metric}_aa])"
            
            if metric == "aatt":
                return "f'{Seq.tsn[scen.get().aatt]}'"
            
            return f"scen.get().{metric}"
        
        return f"""
@render.ui
def {metric}_{phase}():

    # No Variant
    if vnt.get().empty:
        return
    
    # Not a NAGNAG
    if vnt.get().not_nagnag or (not scen.get().canon):
        return Interface.result("-")
    
    return Interface.result({val()})
        """

    @staticmethod
    def render_results_metric_table(phase:str, metric:str) -> str:
        i = Fmt.phases(phase, ind=True)
        m = "aa" if (metric == "aat") else metric

        return f"""
@render.data_frame
def {metric}_{phase}_dataframe():
    
    # Invalid
    if vnt.get().empty or vnt.get().not_nagnag or (not scen.get().canon):
        return
    
    # Data
    if vnt.get().effect == None:
        return
    
    e = f"{{vnt.get().effect}}_{{scen.get().allele}}"
    data = Output.tsn_table(metric="{metric}", phase={i}, effects=[e], dbs=Fmt.dbs("k"), fmt=True)
    if data.empty:
        return None

    # row to highlight
    {f"""
    rows = list(data.iloc[:,0])
    i = [rows.index(Seq.aa[a].full_name) for a in scen.get().{metric}_aa if (a != ".") and (a != "_")]
    
    """ if metric in ["ins", "del"] else
    
    f"""

    # Rows + This Scenario
        {f"""
    rows = list(data.index)
    this = scen.get().aatt
        """ if (metric == "aatt") else
        
        f"""
    rows = list(data.iloc[:,0])
    this = scen.get().{m}_ps
        """ if (phase == "p0") else
        
        f"""
    rows = list(zip(data.iloc[:,0], data.iloc[:,1]))
    this = (scen.get().{m}_ps, scen.get().{m}_ds)
        """}

    i = rows.index(this)
    """}
    

    return render.DataTable(
        data = data,
        width = "100%",
        height = "35dvh",
        styles = [
            {{"class" : "table-values numeric"}},
            {{
                "cols" : {[0, 1] if metric in ["vc", "aat"] else 0},
                "class" : "table-index"
            }},
            {{
                "rows" : i,
                "class" : "table-highlight"
            }}
        ]
    )
        """

    @staticmethod
    def render_tsn_freq(phase:str, metric:str) -> str:
        if metric in ["ins", "del"]:
            display = f"Output.rendered_tsn_freq_table(vnt.get(), scen.get(), '{metric}')"
        else:
            display = f"ui.output_data_frame('{metric}_{phase}_freq_table')"

        return f"""
@render.ui
def {metric}_{phase}_freq_desc():
    return Output.tsn_freq_desc(vnt.get(), scen.get(), "{metric}")

@render.ui
def {metric}_{phase}_freq_display():
    return {display}

@render.data_frame
def {metric}_{phase}_freq_table():
    return Output.rendered_tsn_freq_table(vnt.get(), scen.get(), "{metric}")
        """

    # ========== UPLOAD TAB ========== #

    @staticmethod
    def render_download_link(file:dict[str, str]) -> str:
        return f"""
@render.download(filename="{file["id"]}_nagnags.txt")
def download_{file["id"]}():
    data = Output.processed_upload("{file["datapath"]}")
    yield data.to_csv(sep="\t", header=False, index=False)
        """

    # ========== DATA TABS ========== #

    @staticmethod
    def render_download_existing(id:str, file:str) -> str:
        return f"""
@render.download(filename="{file}")
def {id}():
    return f"{{ROOT}}/data/export/{file}"
        """    

    # ===== Variant Effects ===== #
    
    @staticmethod
    def render_variant_effects_tables(eff:str, freq:str) -> str:
        h = "100%" if (eff == "effects_summary") else "62.5dvh"

        return f"""
@render.data_frame
def {eff}_{freq}_dataframe():
    return render.DataTable(
        data = Output.variant_effects_table(eff="{eff}", freq="{freq}", fmt=True),
        width = "100%",
        height = "{h}",
        styles = [
            {{"class" : "table-values numeric"}},
            {{
                "cols" : {0 if (eff == "effects_summary") else [0,1]},
                "class" : "table-index"
            }},
        ]
    )
        """

    # ===== Proteome NAGNAGs ===== #

    @staticmethod
    def render_proteome_effects_tables(metric:str, phase:str="", eff:str="") -> str:
        i = Fmt.phases(phase, ind=True)
        effect = ["exp_stoch", "exp_splice", "canon"] if (eff == "canon") else [eff]        
        db = Fmt.dbs("k", vnt=(eff != "canon"))
        h = "100%" if (phase == "p0") or (metric == "aatt") else "50dvh"

        return f"""
@render.data_frame
def {metric}_{phase}_{eff}_dataframe():
    data = Output.tsn_table(metric="{metric}", phase={i}, effects={effect}, dbs={db}, freq=input.tsn_data_freq(), fmt=True)
    if data.empty:
        return None

    return render.DataTable(
        data = data,
        width = "100%",
        height = "{h}",
        styles = [
            {{"class" : "table-values numeric"}},
            {{
                "cols" : {[0, 1] if metric in ["vc", "aat"] else 0},
                "class" : "table-index"
            }},
        ]
    )
        """

    # ===== Observed NAGNAGs ===== #

    @staticmethod
    def render_observed_sites_tables(db:str, phase:str=None) -> str:
        fn = "obs_sites_table" if (phase == None) else "obs_scens_table"
        p = "sites" if (phase == None) else phase
        i = {
            "nc" : -1,
            "p0" : 0,
            "p1" : 1,
            "p2" : 2,
            None : None
        }[phase]

        return f"""
@output
@render.data_frame
def {db}_{p}_dataframe():
    return render.DataTable(
        data = Output.{fn}(db="{db}", {f"phase={i}, " if (phase != None) else ""}effects=list(input.effects_to_show()), fmt=True),
        width = "100%",
        height = "57.5dvh",
        styles = [
            {{"class" : "table-values"}},
            {f"""
            {{"cols" : 0, "class" : "table-index"}}
            """ if (db != "canon") else ""}
        ]
    )
        """
