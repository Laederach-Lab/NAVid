from ..lib import ROOT, pd, np, Inputs, Outputs, Session, reactive, render, ui
from ..utils import Seq, Scenario, Variant
from ..uix import State, Interface, Output, Response, Fmt

def server(input:Inputs, output:Outputs, session:Session):
    
    # ----- Variables ----- #

    state = reactive.value(State())
    vnt = reactive.value(Variant())
    scen = reactive.value(Scenario())

    uploads = reactive.value([])
    curr_toasts = reactive.value([])

    env = globals().copy()
    env.update(locals())

    # ===== Upload ===== #
    
    # Upload
    @reactive.effect
    @reactive.event(input.file_upload)
    def _():

        names = [f["id"] for f in uploads.get()]
        files = uploads.get()
        new_files = []

        toasts = curr_toasts.get()

        for f in input.file_upload():

            # sanitize filename
            f["name"] = "".join(filter(lambda x: str(x).isalnum() or (x == "_"), str(f["name"]).split(".")[0]))

            # get instance (files with the same name)
            f["inst"] = names.count(f["name"])
            f["id"] = f"{f["name"]}{"" if (f["inst"] == 0) else f"_{f["inst"]}"}"

            # empty
            if f["size"] == 0:
                ui.show_toast(
                    Interface.toast(
                        "File Size Error", 
                        fr"<b>{f["name"]}</b> is empty!", 
                        id=f"empty_{f["id"]}"
                    )
                )
                toasts.append(f"empty_{f["id"]}")
                continue
            
            # too big
            if f["size"] > 5e+7:
                bytes = {1e+9 : "GB", 1e+6 : "MB", 1e+3 : "KB", 1e+0 : "B"}
                for b, unit in bytes.items():
                    if f["size"] > b: d = b
                    else: break
                
                ui.show_toast(
                    Interface.toast(
                        "File Size Error",
                        [fr"<b>{f["name"]}</b> ({round(f["size"] / d, 1)} {unit}) is too big!",
                         "Maximum file size: 50MB"],
                        id=f"gt50mb{f["id"]}"
                    )
                )
                toasts.append(f"gt50mb{f["id"]}")
                continue

            # too few columns
            data = pd.read_csv(f["datapath"], sep="\t", header=None, comment="#", skip_blank_lines=True)
            if len(data.columns) < 4:
                ui.show_toast(
                    Interface.toast(
                        "File Format Error",
                        [fr"<b>{f["name"]}</b> only has {len(data.columns)} columns!",
                         r"4 columns (<i>chromosome</i>, <i>position</i>, <i>strand</i>, <i>alternate allele</i>) are expected"],
                        id=f"lt4_columns_{f["id"]}"
                    )
                )
                toasts.append(f"lt4_columns_{f["id"]}")
                continue

            # invalid variants
            ni = 0
            for c, p, _, a in data.values[:, :4]:
                if ((Seq.format_chrom(c) == None) or 
                    (not Seq.valid_coord(p)) or 
                    (Seq.format_base(a) == None)):
                    ni += 1

            if ni > 0:
                pl = "s" if (ni > 1) else ""
                verb = "have" if (ni > 1) else "has"
                ui.show_toast(
                    Interface.toast(
                        "Warning",
                        [fr"{ni} variant{pl} in <b>{f["name"]}</b> {verb} invalid data type(s)",
                        "Their effect is marked as <b>invalid</b>"],
                        toast_type="warning", id=f"invalid_data_{f["id"]}"
                    )
                )
                toasts.append(f"invalid_data_{f["id"]}")

            # add file to lists
            new_files.append(f)
            names.append(f["name"])

        # update variables
        uploads.set(files + new_files)
        curr_toasts.set(toasts)

        for file in new_files:
            exec(Response.render_download_link(file), env)
    
    # Clear
    @reactive.effect
    @reactive.event(input.clear_uploads)
    def _():
        uploads.set([])

        # close toasts
        for tid in curr_toasts.get():
            ui.hide_toast(tid)
        curr_toasts.set([])

    # Submit
    @render.ui
    @reactive.event(input.file_upload, input.clear_uploads)
    def download_processed_files():
        return ui.div(*[Interface.download_link(f["id"]) for f in uploads.get()])
    
    # Examples
    @render.text
    def input_columns():
        return "chrom	pos	strand	alt_allele"

    @render.text
    def example_input():
        return """
chr1	179889308	+	A
chr1	182602483	-	G
chr1	183258105	+	G
chr1	197122304	-	G
chr1	202434651	+	G
chr1	205775089	-	A
        """

    @render.text
    def example_output():
        return """
chr1	179889308	+	A	destroy	TAACAG	TAGCAG
chr1	182602483	-	G	alter	TAGGAG	TAGTAG
chr1	183258105	+	G	alter	GAGGAG	CAGGAG
chr1	197122304	-	G	create	TAGTAG	TAATAG
chr1	202434651	+	G	create	TAGCAG	TAACAG
chr1	205775089	-	A	destroy	GAACAG	GAGCAG
        """
    
    # ===== Input ===== #

    exec(Response.reactive_input_text("chrom"), env)
    exec(Response.reactive_input_text("var_pos"), env)
    exec(Response.reactive_input_text("alt_allele"), env)
    exec(Response.reactive_input_text("ss5_pos"), env)
    exec(Response.reactive_input_radio("strand"), env)
    exec(Response.reactive_input_radio("index"), env)

    # Switch Strand
    @reactive.effect
    @reactive.event(input.switch_strand)
    def _():
        s = {"+":"-", "-":"+"}[input.strand()]
        state.get().update_strand(s)
        ui.update_radio_buttons(id="strand", selected=s)

        a = Seq.rc(input.alt_allele())
        state.get().update_alt_allele(a)
        ui.update_text(id="alt_allele", value=a)

    # Switch Index
    @reactive.effect
    @reactive.event(input.switch_index)
    def _():
        i = {"0":"1", "1":"0"}[str(input.index())]
        d = {"0":1, "1":-1}[str(input.index())]

        state.get().update_index(i)
        ui.update_radio_buttons(id="index", selected=i)

        if input.ss5_pos():
            p = state.get().update_ss5_pos(int("".join(input.ss5_pos().split(","))) + d)
            ui.update_text(id="ss5_pos", value=p)
        
        if input.var_pos():
            p = state.get().update_var_pos(int("".join(input.var_pos().split(","))) + d)
            ui.update_text(id="var_pos", value=p)

    def submit():

        # close toasts
        for tid in curr_toasts.get():
            ui.hide_toast(tid)

        toasts = []

        errors, v = state.get().get_variant()
        
        # error messages
        if errors:
            for c, e in enumerate(errors):
                ui.show_toast(Interface.toast("Input Error", e, id=f"input_error_{c}"))
                toasts.append(f"input_error_{c}")
            curr_toasts.set(toasts)

        vnt.set(v)
        scen.set(state.get().get_scenario(v))

        a = (v.effect == "create")
        l = "Alternate Allele" if a else "Reference Allele"
        ui.update_switch(id="allele", label=l)
        state.get().update_allele(a)

        # no 5'ss
        if not errors and vnt.get().not_scenario:
            ui.show_toast(
                Interface.toast(
                    "Missing 5' splice site", 
                    fr"Specify the <b>5' Splice Site Location</b> in order to view proteome effects",
                    toast_type="dark", id="missing_5ss_toast"
                )
            )
            toasts.append("missing_5ss_toast")

        curr_toasts.set(toasts)

    # Submit
    @reactive.effect
    @reactive.event(input.submit)
    def _():
        state.get().reset_example()
        submit()
    
    # Reset
    @reactive.effect
    @reactive.event(input.reset)
    def _():

        state.set(State())
        vnt.set(Variant())
        scen.set(Scenario())

        # Inputs
        ui.update_text(id="chrom", value="")
        ui.update_text(id="var_pos", value="")
        ui.update_text(id="alt_allele", value="")
        ui.update_text(id="ss5_pos", value="")
        
        ui.update_radio_buttons(id="strand", selected="+")
        ui.update_radio_buttons(id="index", selected="1")
    
    # Get Example
    @reactive.effect
    @reactive.event(input.example)
    def _():

        # random scenario
        scen = Output.vnt_scens().iloc[np.random.randint(0, 76745)]

        # update state
        c = state.get().update_chrom(scen["chrom"])
        ui.update_text(id="chrom", value=c)

        p = state.get().update_var_pos(scen["pos"])
        ui.update_text(id="var_pos", value=p)

        a = state.get().update_alt_allele(scen["allele"])
        ui.update_text(id="alt_allele", value=a)

        s = state.get().update_strand(scen["strand"])
        ui.update_radio_buttons(id="strand", selected=s)

        u = state.get().update_ss5_pos(scen["upstream_end"])
        ui.update_text(id="ss5_pos", value=u)

        i = state.get().update_index(0)
        ui.update_radio_buttons(id="index", selected=i)
        
        state.get().set_example(scen)

        submit()
        
    # ----- Scenario ----- #

    # Phase
    @reactive.effect
    @reactive.event(input.tsns_phase_navset)
    def _():
        state.get().update_phase(input.tsns_phase_navset())
        scen.set(state.get().get_scenario(vnt.get()))

    # Allele
    @reactive.effect
    @reactive.event(input.allele)
    def _():

        # alter
        if vnt.get().effect == "alter":
            a = "Alternate Allele" if input.allele() else "Reference Allele"
            ui.update_switch(id="allele", label=a)
            state.get().update_allele(input.allele())
            scen.set(state.get().get_scenario(vnt.get()))
        
        # create
        elif vnt.get().effect == "create":
            if not input.allele():
                ui.show_toast(
                    ui.toast(
                        "Information only available for the alternate allele of NAGNAG-creating variants!",
                        type = "primary",
                        duration_s = 3,
                        position = "bottom right",
                        closable = True
                    )
                )
            ui.update_switch(id="allele", label="Alternate Allele", value=True)
        
        # destroy
        elif vnt.get().effect == "destroy":
            if input.allele():
                ui.show_toast(
                    ui.toast(
                        "Information only available for the reference allele of NAGNAG-destroying variants!",
                        type = "primary",
                        duration_s = 3,
                        position = "bottom right",
                        closable = True
                    )
                )
            ui.update_switch(id="allele", label="Reference Allele", value=False)

    # ----- Variant Look-Up ----- #
    
    # Search
    @reactive.effect
    @reactive.event(input.search)
    def _():
        ui.modal_show(Interface.search_modal())

    # Submit Search
    @reactive.effect
    @reactive.event(input.submit_search)
    def _():
        v,_ = State.find_variant(input.db_select(), input.acc_num())

        # error
        if v == None:
            return

        # use variant
        c,p,_,s = v

        c = state.get().update_chrom(c)
        ui.update_text(id="chrom", value=c)

        p = state.get().update_var_pos(p)
        ui.update_text(id="var_pos", value=p)

        s = state.get().update_strand(s)
        ui.update_radio_buttons(id="strand", selected=s)

        state.get().update_index(0)
        ui.update_radio_buttons(id="index", selected="0")

        state.get().update_ss5_pos(None)
        ui.update_text(id="ss5_pos", value="")

    # Error Message
    @render.ui
    @reactive.event(input.submit_search)
    def search_error():
        _, err = State.find_variant(input.db_select(), input.acc_num())
        return ui.p(err, style="color:red; line-height:1; font-weight:550") if err else None
    
    # Alternate Alleles Select
    @render.ui
    @reactive.event(input.submit_search)
    def alt_select():
        v, _ = State.find_variant(input.db_select(), input.acc_num())
        return Interface.alt_allele_select(v[2]) if v else None

    for alt in ["A", "C", "G", "T", "U", "N"]:
        exec(Response.reactive_select_allele(alt), env)

    # ===== Search Results ===== #

    @render.ui
    def search_results_output():
        return Interface.search_results_output(vnt.get())
        
    # ----- Variant ----- #

    @render.ui
    def vnt_id():
        
        # No Variant
        if vnt.get().empty:
            return

        loc = Seq.vnt2str(vnt.get().chrom, vnt.get().var_coord, vnt.get().strand, vnt.get().ref, vnt.get().alt)
        acc_num = state.get().example

        return Interface.result(loc, acc_num)

    @render.ui
    def vnt_type():
        
        # No Variant
        if vnt.get().empty:
            return
        
        t = Seq.sbs_type(vnt.get().ref, vnt.get().alt)
        return Interface.result(f"{vnt.get().ref} > {vnt.get().alt}: {t}")

    @render.ui
    def vnt_affect_status():
        
        # No Variant
        if vnt.get().empty:
            return

        # Not NAGNAG-affecting
        if vnt.get().not_nagnag:
            return Interface.result("Variant does not affect a NAGNAG")

        return Interface.result(f"Variant {vnt.get().effect}s a NAGNAG")

    # ----- Affected NAGNAG ----- #

    @render.ui
    def nagnag_location():

        # No Variant
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")

        return Interface.result(Seq.coord2str(vnt.get().chrom, vnt.get().motif_start, vnt.get().motif_end, vnt.get().strand))

    # Variant Effect
    @render.ui
    def nagnag_effect():

        # No Variant
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")

        if vnt.get().effect not in Fmt.effs("k"):
            return Interface.result("-")

        return Interface.result(Fmt.effs(vnt.get().effect))

    @render.ui
    def nagnag_effect_freq_desc():

        # No Variant        
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")

        e = str(vnt.get().effect).lower()

        df = Output.effect_freqs()["prop"].T[e].sort_values(ascending=False, inplace=False)
        d = Fmt.dbs(df.index.values[0])
        p = df.values[0]

        return Interface.freq(
            f"{Fmt.effs(e)} variants are most common in {d}",
            f"They represent {p:.2%} of observed NAGNAG-affecting {d} variants"
        )

    @render.ui
    def nagnag_effect_freq_display():
        return ui.output_plot("nagnag_effect_freq_plot")

    @render.plot
    def nagnag_effect_freq_plot():
        return Output.nagnag_effect_plot(vnt.get())

    exec(Response.reactive_render_vnt_ids("nagnag_effect"), env)

    @render.ui
    def nag_strength():

        # No Variant
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")

        eff = vnt.get().effect
        pos = vnt.get().var_pos
        ref = vnt.get().ref_motif
        alt = vnt.get().alt_motif
        
        line1, line2 = "", ""

        if eff == "create":

            # original vs. created NAG
            o = alt[3:] if (pos < 3) else alt[:3]
            c = alt[:3] if (pos < 3) else alt[3:]
            
            if o == c:
                compare, s = "as strong as", "="
            else:
                compare, s = {
                    True : ("stronger than", ">"),
                    False : ("weaker than", "<")
                }[(Seq.NAG_STRENGTH[c] - Seq.NAG_STRENGTH[o]) > 0]

            line1 = f"The created NAG is {compare} the original NAG"
            line2 = f"{c} {s} {o}"

        elif eff == "alter":

            # original vs. created NAG
            o = ref[3:] if (pos < 3) else ref[:3]
            r = ref[:3] if (pos < 3) else ref[3:]
            a = alt[:3] if (pos < 3) else alt[3:]

            if o == r:
                compare = "proximal" if (pos < 3) else "distal"
            else:
                compare, sOA = {
                    True : ("stronger", ">"),
                    False : ("weaker", "<")
                }[(Seq.NAG_STRENGTH[r] - Seq.NAG_STRENGTH[o]) > 0]

            change, sRA = {
                True : ("increases", ">"),
                False : ("decreases", "<")
            }[(Seq.NAG_STRENGTH[a] - Seq.NAG_STRENGTH[r]) > 0]

            line1 = f"The variant {change} the strength of the {compare} NAG",
            line2 = f"{f"{r} {sOA} {o}; " if (o != r) else ""}{a} {sRA} {r}"

        elif eff == "destroy":

            # original vs. destroyed NAG
            o = ref[3:] if (pos < 3) else ref[:3]
            d = ref[:3] if (pos < 3) else ref[3:]
            
            if o == d:
                line1 = f"The {"proximal" if (pos < 3) else "distal"} NAG is destroyed"
                line2 = f"Both NAGs are of equal strength: ({o} = {d})"
            else:
                compare, s = {
                    True : ("stronger", ">"),
                    False : ("weaker", "<")
                }[(Seq.NAG_STRENGTH[d] - Seq.NAG_STRENGTH[o]) > 0]

                line1 = f"The {compare} NAG is destroyed"
                line2 = f"{d} {s} {o}"

        return Interface.result(line1, line2)

    # NAGNAG Motif
    @render.ui
    def motif_ref():

        # No Variant        
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")
        
        return Interface.result(vnt.get().ref_motif)

    @render.ui
    def motif_alt():

        # No Variant        
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")
        
        return Interface.result(vnt.get().alt_motif)

    @render.ui
    def motif_freq_desc():
        
        # No Variant
        if vnt.get().empty:
            return
        
        # Not a NAGNAG
        if vnt.get().not_nagnag:
            return Interface.result("-")
        
        e = str(vnt.get().effect).lower()
        r = vnt.get().ref_motif
        a = vnt.get().alt_motif

        df = Output.event_freqs().loc[(e, r, a)].sort_values(inplace=False, ascending=False)
        
        # No other occurances
        if sum(df.values) == 0:
            return Interface.freq(
                "A new variant!",
                f"{r} > {a} is not listed in dbSNP or ClinVar"
            )

        # Other occurances
        counts = df["count"]
        counts = counts[counts > 0]

        db = counts.index.values[0]
        n = counts.values[0]
        one = (n == 1)

        # define other values
        match len(counts) - 1:
            case 0:
                other = ""
            case 1:
                other = Fmt.dbs(counts.index.values[1])
            case 2:
                other = " and ".join([Fmt.dbs(d) for d in counts.index.values[1:]])
            case _:
                index = counts.index.values[1:]
                other = ", ".join([Fmt.dbs(d) for d in index[:-1]]) + ", and " + Fmt.dbs(index[-1])

        return Interface.freq(
            f"{n:.0f} case{"" if one else "s"} of {r} > {a} {"has" if one else "have"} been recorded in {Fmt.dbs(db)}",
            f"This type of motif variant also occurs in {other}" if other else ""
        )

    @render.ui
    def motif_freq_display():
        return ui.div(
            ui.output_data_frame("motif_freq_table"),
            style = "width:100%"
        )

    @render.data_frame
    def motif_freq_table():
        return render.DataTable(
            data = Output.motif_table(vnt.get()),
            height = None,
            styles = [
                {"class" : "table-values numeric"},
                {"cols" : 0, "class" : "table-index"}
            ]
        )

    exec(Response.reactive_render_vnt_ids("motif"), env)

    # ----- Proteome Effects ----- #

    @render.image
    def splice_junction_image():
        fp = f"{ROOT}/utils/splice_junction.svg"
        Output.splice_junction_plot(vnt.get()).savefig(fp)
        return dict(src=fp, style="max-width:100%; max-height:350px")

    for p in ["p0", "p1", "p2"]:
        
        # Values
        for m in ["vc_ps", "vc_ds", "aa_ps", "aa_ds", "aatt", "ins", "del"]:
            exec(Response.render_tsn_result(p, m), env)
    
        # Frequency + Output Tables
        for m in ["vc", "aat", "aatt", "ins", "del", "ins_a0", "ins_a1", "del_a0", "del_a1"]:
            exec(Response.render_tsn_freq(p, m), env)
            exec(Response.reactive_render_vnt_ids(f"{m}_{p}"), env)
            exec(Response.render_results_metric_table(p, m), env)
        
    # ===== Variant Effects ===== #

    for e in ["effects_summary", "create_motifs", "alter_motifs", "destroy_motifs"]:
        for f in ["count", "prop"]:
            exec(Response.render_variant_effects_tables(e, f), env)

    exec(Response.render_download_existing("download_variants_tar_gz", "variant_effects.tar.gz"), env)

    # ====== Proteome Effects ===== #
    
    for m,phases in Fmt.metrics(phases=True):
        for p in phases:
            for e in ["canon"] + Fmt.scen_effs("k"):
                exec(Response.render_proteome_effects_tables(metric=m, phase=p, eff=e), env)

    exec(Response.render_download_existing("download_proteome_tar_gz", "proteome_effects.tar.gz"), env)

    # ===== Observed NAGNAGs ===== #

    @reactive.effect
    @reactive.event(input.site_type_to_show)
    def _():
        if input.site_type_to_show() == "sites":
            ui.update_checkbox_group(id="effects_to_show", choices=Fmt.effs("d"), selected=Fmt.effs("k"))
        elif input.site_type_to_show() == "scenarios":
            ui.update_checkbox_group(id="effects_to_show", choices=Fmt.scen_effs("d"), selected=Fmt.scen_effs("k"))

    @render.ui
    def observed_nagnags_panel():
        return ui.navset_pill(
            Interface.database_data_panel("canon", input.site_type_to_show()),
            Interface.database_data_panel("dbSNP_common", input.site_type_to_show()),
            Interface.database_data_panel("dbSNP_rare", input.site_type_to_show()),
            Interface.database_data_panel("ClinVar", input.site_type_to_show()),
            id = "observed_nagnags_navset"
        )

    for d in ["canon", "dbSNP_common", "dbSNP_rare", "ClinVar"]:
        exec(Response.render_observed_sites_tables(d), env)
        for p in ["nc", "p0", "p1", "p2"]:
            exec(Response.render_observed_sites_tables(d, p), env)
    
    exec(Response.render_download_existing("download_observed_tar_gz", "observed_nagnags.tar.gz"), env)
