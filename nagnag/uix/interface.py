from .format import Fmt
from ..lib import faicons, ui
from ..utils import Variant

class Interface:

    ''' User interface components (input, reactive, modules, pages, etc.) '''

    def __init__(self):
        pass

    # ===== Input Components ===== #
    
    @staticmethod
    def input_text(id:str, label:str, tooltip:str="", initial:str="", example:str=""):
        l = ui.p(label, style="font-weight:600; padding:0; margin:0"),
        return ui.input_text(
            id = id,
            label = ui.tooltip(l, tooltip) if tooltip else l,
            value = initial,
            placeholder = example,
            update_on = "blur"
        )
    
    @staticmethod
    def input_radio(id:str, label:str, choices:list|dict, tooltip:str="", initial:str="", inline:bool=True):
        l = ui.p(label, style=f"font-weight:600; line-height:{None if inline else 2.5}; padding:0; margin:0")
        return ui.input_radio_buttons(
            id = id,
            label = ui.tooltip(l, tooltip) if tooltip else l,
            choices = choices,
            selected = initial,
            inline = inline
        )
    
    @staticmethod
    def input_check(id:str, label:str, choices:list|dict, tooltip:str="", initial:str|list[str]=""):
        l = ui.p(label, style="font-weight:600; line-height:2.5; padding:0; margin:0")
        return ui.input_checkbox_group(
            id = id,
            label = ui.tooltip(l, tooltip) if tooltip else l,
            choices = choices,
            selected = list(initial),
        )
    
    @staticmethod
    def input_button(id:str, label:str, tooltip:str="", class_:str="button", style:str=""):
        return ui.input_action_button(
            id = id,
            label = ui.tooltip(label, tooltip) if tooltip else label,
            class_ = class_,
            style = style
        )
    
    @staticmethod
    def download_button(id:str, label:str, tooltip:str="", class_:str="", div_class:str="download-div"):
        return ui.div(
            ui.download_button(
                id = id,
                label = ui.tooltip(label, tooltip) if tooltip else label,
                style = "padding: 12px 24px",
                class_ = f"button {class_}"
            ),
            class_ = div_class
        )

    # ===== Output Components ===== #

    @staticmethod
    def card_w_header(header:ui.TagChild, *content:ui.TagChild):
        return ui.card(
            ui.card_header(
                ui.div(header, style="margin:16px"),
                class_ = "card-header"
            ),
            ui.div(*content, style="margin:0; padding:0"),
        )

    @staticmethod
    def result(value:str="", sub:str=None, color:str="black", subcolor:str="black"):
        if sub:
            return ui.div(
                ui.div(
                    ui.h1(value, class_="results-value"),
                    style = f"display:flex; width:100%; justify-content:center; padding:15px 0px 7px; color:{color}"
                ),
                ui.div(
                    ui.h4(sub, class_="results-value sub"),
                    style = f"display:flex; width:100%; justify-content:center; color:{subcolor}"
                ),
            )

        return ui.h1(value, class_="results-value", style=f"color:{color}")

    @staticmethod
    def freq_div(desc_id:str, display_id:str, button_id:str, display):
        return ui.div(
            ui.div(
                ui.output_ui(desc_id),
                style = "margin:0px 0px 12px"
            ),
            ui.div(
                display(display_id),
                style = "margin:0px 0px 12px"
            ),
            ui.div(
                Interface.input_button(
                    id = button_id,
                    label = "Observed Variant IDs",
                    class_ = "button",
                ),
                class_ = "ids-div"
            )
        )

    @staticmethod
    def freq(line1:str="", line2:str=""):

        if (not line1) and (not line2):
            return
        
        return ui.div(
            ui.h4(line1, class_="freq-value-1"),
            ui.h5(line2, class_="freq-value-2")
        )
    
    @staticmethod
    def download_link(name:str):
        return ui.div(
            ui.download_link(id=f"download_{name}_nagnags", label=name),
            style = "padding: 6px 0px 0px"
        )

    # ===== Modules ===== #

    @staticmethod
    def toast(title:str, message:str|list[str], toast_type:str="error", duration:int=None, id:str=None):
        icon, color = {
            "error" : ("circle-xmark", "#C10000"), 
            "warning" : ("triangle-exclamation", "#F9B939"),
            "dark" : ("circle-info", "#000000"),
        }[toast_type]

        if isinstance(message, str):
            message = [message]

        return ui.toast(
            *[ui.div(ui.HTML(s)) for s in list(message)],
            header = ui.toast_header(
                title = ui.div(title, class_="toast-header-text"),
                icon = faicons.icon_svg(icon, fill=color)
            ),
            type=toast_type,
            duration_s=duration,
            id=id
        )

    @staticmethod
    def search_modal():
        return ui.modal(

            ui.div(
                # Accession Number ('acc_num')
                ui.div(
                    Interface.input_text(id="acc_num", label="Accession Number", example="rs100000",
                                        tooltip="dbSNP of ClinVar variant accession number"),
                ),

                # Database ('db_select')
                ui.div(
                    Interface.input_radio(id="db_select", label="Variant Database",
                                        tooltip="Database in which variant has been identified",
                                        choices=["dbSNP", "ClinVar"], initial="dbSNP"),
                    style = "margin:15px 0px"
                ),

                # Submit ('submit_search'), Error Message ('search_error')
                ui.div(

                    # Submit ('submit_search')
                    ui.div(
                        Interface.input_button(id="submit_search", label="Search", style="padding:6px 8px"),
                        style = "float:left; width:fit-content; height:fit-content; margin-right:12px"
                    ),

                    # Search Error ('search_error')
                    ui.div(
                        ui.output_ui("search_error"),
                        style = "display:flex; width:100%; height:20px; padding-top:2px"
                    ),
                    style = "display:flex; width:100%; height:fit-content; align-items:center"
                ),

                # Alternate Alleles
                ui.div(
                    ui.output_ui("alt_select"),
                )
            ),
            
            size = "m",
            footer = None,
            easy_close=True,
        )
    
    @staticmethod
    def alt_allele_select(alt:list[str]):

        def button(a:str):
            return ui.div(
                Interface.input_button(id=f"select_alt{a}", label=a, style="padding:4px 16px"),
                style = "float:left; width:fit-content; height:fit-content; margin-right:16px"
            )
        
        return ui.div(
            ui.div(
                ui.tooltip(
                    ui.p("Alternate Allele", style="line-height:1; font-weight:bold"),
                    "Alternate alleles, based on the '+' strand"
                ),
                style = "height:20px; margin-top:16px; margin-bottom:6px"
            ),

            ui.div(
                *[button(a) for a in alt],
                style = "display:flex; width:100%; height:32px"
            ),

            style = "margin-top:12px"
        )

    @staticmethod
    def download_table(id:str, label:str, style:str="width:fit-content", button_class:str="download_button", download:bool=True):
        return ui.div(
            ui.div(
                ui.output_data_frame(f"{id}_dataframe"),
                style = style
            ),
            Interface.download_button(
                id=f"{id}_download",
                label=f"Download {label}",
                class_=button_class
            ) if download else None,
        )

    @staticmethod
    def vnt_ids_modal(id:str):
        def panel(db:str):
            return ui.nav_panel(
                Fmt.dbs(db),
                ui.div(
                    ui.div(
                        ui.output_data_frame(f"{id}_{db}_vnts_table"),
                        style = "height:47.5dvh"
                    ),
                    Interface.download_button(
                        id = f"{id}_{db}_vnts_download",
                        label = "Download Variants"
                    ),
                    style = "margin-top:12px"
                )
            )
        
        return ui.modal(
            ui.navset_pill(
                panel("dbSNP_common"),
                panel("dbSNP_rare"),
                panel("ClinVar"),
            ),
            title = None,
            size = "l",
            easy_close = True,
            footer = None
        )

    @staticmethod
    def navset_pill_divider():
        return ui.nav_control(
            ui.div(
                ui.div(
                    style = "padding-top:10%; margin:0px 15px; height:80%; border-left:1.5px solid lightgray"
                ),
                style = "display:flex; align-items:center; height:100%"
            )
        )

    # ===== Layouts ===== #

    @staticmethod
    def search_results_output(vnt:Variant=Variant()):

        def variant_tab():
            return ui.nav_panel(
                "Variant",

                ui.div(
                    ui.layout_columns(
                        Interface.results_card(id="vnt_id", title="Variant ID", tall=True),
                        Interface.results_card(id="vnt_type", title="Variant Type", tall=True,
                                               tooltip="The effect of the variant on the DNA sequence"),
                        Interface.results_card(id="vnt_affect_status") if (not vnt.empty) else None,
                        col_widths = (7, 5, 12) if (not vnt.empty) else (7, 5),
                        gap = "18px"
                    ),
                    class_ = "page-div"
                ),

                value = "variant"
            )
        
        def nagnag_tab():
            return ui.nav_panel(
                "Affected NAGNAG",

                ui.div(
                    ui.layout_columns(
                        Interface.results_card(id="nagnag_location", title="Location",
                                               tooltip="0-index chromosomal position of the affected NAGNAG"),
                        Interface.results_card(id="nagnag_effect", title="Variant Effect",
                                               tooltip="How the variant affects the NAGNAG",
                                               freq=True, id_button=True, tall=True,),
                        Interface.results_card(id="motif<SUB>", title="Motif",
                                               tooltip="The NAGNAG motif in the reference and alternate alleles",
                                               labels={"ref":"Reference", "alt":"Alternate"}, freq=True, id_button=True),
                        Interface.results_card(id="nag_strength", title="NAG Strength",
                                               tooltip="Effect of the variant on the NAG strengths", tall=True),
                        col_widths = (12, 6, 6, 12),
                        gap = "18px"
                    ),
                    class_ = "page-div"
                ),

                value = "affected_nagnag"
            )

        def proteome_tab():

            def phase_tab(p:int):

                tabs = [
                    ("vc", "Variable Codons"),
                    ("aat", "Variable Amino Acids"),
                    ("aatt", "Transition Type"),
                    ("ins", "Amino Acid Insertions")
                ]
                if p != 0:
                    tabs.append(("del", "Amino Acid Deletions"))

                def data_tab(metric:str, title:str):
                    return ui.nav_panel(
                        title,
                        ui.div(
                            ui.output_data_frame(f"{metric}_p{p}_dataframe"),
                            style = "margin-top:12px; "
                        ),
                        value = metric
                    )
                
                return ui.nav_panel(
                    f"Phase {p}",

                    ui.div(
                        ui.layout_columns(

                            # Results
                            Interface.results_card(id=f"vc<SUB>_p{p}", title="Variable Codons",
                                                   tooltip="The mRNA codons that differ between the proximally- and distally-spliced transcript isoforms",
                                                   labels={"ps":"Proximal", "ds":"Distal"},
                                                   freq=True, id_button=True),
                            Interface.results_card(id=f"aat<SUB>_p{p}", title="Variable Amino Acids",
                                                   tooltip="The amino acids that differ between the protein isoforms that result from proximally and distally splicing",
                                                   labels={"ps":"Proximal", "ds":"Distal"},
                                                   freq=True, id_button=True),
                            Interface.results_card(id=f"aatt_p{p}", title="Transition Type",
                                                   tooltip="The transition type describing the relationship between the amino acid sequences of the two protein isoforms",
                                                   freq=True, id_button=True, tall=True),
                            Interface.results_card(id=f"ins_p{p}", title="Inserted Amino Acids",
                                                   tooltip="The amino acid(s) absent in the distal protein isoform but present in the proximal protein isoform",
                                                   freq=True, id_button=True),
                            Interface.results_card(id=f"del_p{p}", title="Deleted Amino Acids",
                                                   tooltip="The amino acid(s) present in the distal protein isoform but absent in the proximal protein isoform",
                                                   freq=True, id_button=True),

                            # Data tables
                            Interface.card_w_header(
                                ui.h4("Statistics", class_="results-card-title"),
                                ui.navset_pill(
                                    *[data_tab(v, t) for v,t in tabs],
                                    id = f"phase{p}_data_navset"
                                ),
                            ),

                            col_widths = (4, 4, 4, 6, 6, 12),
                            gap = "18px"
                        ),

                        class_ = "page-div",
                    ),

                    value = f"p{p}"
                )

            value = vnt.effect == "create"
            label = {
                False : "Reference Allele",
                True : "Alternate Allele",
            }[value]

            return ui.nav_panel(
                "Proteome Effects",
                
                ui.div(
                    Interface.card_w_header(
                        ui.h5("Splice Junction", class_="results-card-title"),
                        ui.div(
                            ui.output_image("splice_junction_image", width="fit-content", height="fit-content", fill=True),
                            style = "display:flex; justify-content:center; width:100%; height:100%"
                        ),
                    ),
                    class_ = "page-div",
                ),
                
                ui.div(
                    ui.navset_pill(
                        phase_tab(0),
                        phase_tab(1),
                        phase_tab(2),

                        ui.nav_spacer(),
                        ui.nav_control(
                            ui.div(
                                ui.input_switch(
                                    id = "allele",
                                    label = label,
                                    value = value
                                ),
                                class_ = "switch",
                            )
                        ),
                        id = "tsns_phase_navset"
                    ),

                    class_ = "page-div"
                ),
                
                value = "proteome_effects"
            )

        if vnt.empty or vnt.not_nagnag:
            tabs = [variant_tab()]
        elif vnt.not_scenario:
            tabs = [variant_tab(), nagnag_tab()]
        else:
            tabs = [variant_tab(), nagnag_tab(), proteome_tab()]
        
        return ui.navset_pill(
            *tabs,
            id = "search_results_navset"
        )

    @staticmethod
    def results_card(id:str, title:str="", tooltip:str="", labels:dict[str]=None,
                     freq:bool=False, id_button:bool=False, tall:bool=None):

        h5 = ui.h5(title, class_="results-card-title")
        header = ui.tooltip(h5, tooltip) if tooltip else h5
        
        card_id = id.replace("<SUB>", "")
        id = id.replace("aat", "aa") if (id.startswith("aat<SUB>")) else id

        def content():
            if labels:
                sub = {i : (k,v) for i,(k,v) in enumerate(labels.items())}
                return ui.div(
                    
                    # Values
                    ui.div(
                        ui.div(
                            ui.output_ui(id.replace("<SUB>", f"_{sub[0][0]}")),
                            class_ = "results-compare value"
                        ),
                        ui.div(
                            ui.output_ui(id.replace("<SUB>", f"_{sub[1][0]}")),
                            class_ = "results-compare value"
                        ),
                        style = "display:table-row"
                    ),

                    # Labels
                    ui.div(
                        ui.div(ui.p(sub[0][1]), class_ = "results-compare label"),
                        ui.div(ui.p(sub[1][1]), class_ = "results-compare label"),
                        style = "display:table-row"
                    ),

                    style = "display:table; width:100%"
                )

            return ui.div(
                ui.div(
                    ui.output_ui(card_id),
                    class_ = f"results-single {"normal" if (not tall) else "tall"}",
                ),
                class_ = "results-content"
            )
        
        def freq_accordion():
            return ui.accordion(
                ui.accordion_panel(
                    ui.p("Frequency", class_="freq-accordion-title"),
                    ui.card(
                        Interface.freq_div(
                            desc_id = f"{card_id}_freq_desc",
                            display_id = f"{card_id}_freq_display",
                            button_id = f"{card_id}_obs_vnt_ids",
                            display = ui.output_ui
                        ),
                        full_screen = True
                    ) if (id[:3] not in ["ins", "del"]) else 
                    ui.card(
                        ui.output_ui(f"{card_id}_freq_display"),
                        full_screen = True
                    ),
                    value = "freq"
                ),
                open = False,
                width = "100%"
            
            ) if freq else None
        
        if title:
            return Interface.card_w_header(header, content(), freq_accordion())
        
        return ui.card(
            ui.div(
                content(),
                freq_accordion(),
                style = "display:flex; flex-direction:column; width:100%; height:100%"
            )
        )
    
    @staticmethod
    def effects_data_panel(eff:str, name:str):

        def card(freq:str, title:str, label:str):
            return Interface.card_w_header(
                ui.h5(title, class_="results-card-title"),
                Interface.download_table(id=f"{eff}_{freq}", label=label, download=False),
            )

        return ui.nav_panel(
            name,
            ui.layout_columns(
                card("count", "Number of Variants", "Counts"),
                card("prop", "Proportion of Variants", "Proportions"),
                col_widths = (6, 6),
                gap = "18px",
                style = "padding:0px; margin:0px"
            ),
            value = eff
        )

    @staticmethod
    def tsn_metric_data_panel(metric:str):

        def phase_tab(p:int):
        
            def effect_tab(e:str):
                return ui.nav_panel(
                    Fmt.scen_effs(e, "hg38 Reference"),
                    ui.card(
                        Interface.download_table(id=f"{metric}_{p}_{e}", label="Statistics", download=False),
                        style = "margin:0px"
                    ),
                    value = e
                )

            return ui.nav_panel(
                Fmt.phases(p),
                ui.div(
                    ui.navset_pill_list(
                        effect_tab("canon"),
                        effect_tab("create_alt"),
                        effect_tab("alter_ref"),
                        effect_tab("alter_alt"),
                        effect_tab("destroy_ref"),

                        widths = (3, 9),
                        id = f"{metric}_{p}_effect_navset",
                    ),
                    class_ = "page-div list"
                ),
                value = p
            )

        return ui.nav_panel(
            Fmt.metrics(metric),
            ui.div(
                ui.navset_pill(
                    *[phase_tab(p) for p in Fmt.metrics(metric, phases=True)],
                    id = f"{metric}_phase_navset"
                ),
                class_ = "page-div"
            ),
            value = metric
        )

    @staticmethod
    def database_data_panel(db:str, site:str):

        if site == "sites":
            return ui.nav_panel(
                Fmt.dbs(db),

                ui.div(
                    ui.card(
                        Interface.download_table(id=f"{db}_sites", label="Sites", download=False),
                        style = "margin:0px"
                    ),

                    class_ = "page-div"
                ),
            )
        
        def phase_tab(p:str):
            return ui.nav_panel(
                Fmt.phases(p),
                ui.card(
                    Interface.download_table(id=f"{db}_{p}", label="Scenarios", download=False),
                    style = "margin:0px"
                ),
                value = p
            )
        
        return ui.nav_panel(
            Fmt.dbs(db),

            ui.div(
                ui.navset_pill_list(
                    phase_tab("nc"),
                    phase_tab("p0"),
                    phase_tab("p1"),
                    phase_tab("p2"),
                    
                    widths = (2, 10),
                    id = f"{db}_phase_data_navset"
                ),

                class_ = "page-div list"
            ),

            value = db
        )

    @staticmethod
    def about_card(title:str, body:str):
        return Interface.card_w_header(
            ui.h4(title, class_="results-card-title"),
            ui.HTML(f"<p>{body}</p>"),
        )

    @staticmethod
    def table_styles() -> str:

        def style(name:str, ncols:int=1) -> str:
            
            # header
            s = f"""
#{name}_dataframe th {{
        text-align: right;
}}
            """

            # index headers
            for i in range(ncols):
                s += f"""
#{name}_dataframe th:nth-child({i+1}) {{
        text-align: left !important;
}}
                """

            return s

        styles = ""

        # search results
        for p in ["p0", "p1", "p2"]:
            for m in ["vc", "aat", "aatt", "ins", "del"]:
                styles += style(f"{m}_{p}", 2 if m in ["vc", "aat"] else 1) 
        
        # variant effects
        for e in ["effects_summary", "create_motifs", "alter_motifs", "destroy_motifs"]:
            for f in ["count", "prop"]:
                styles += style(f"{e}_{f}", 1 if e == "effects_summary" else 2)
        
        # proteome effects
        for m,phases in Fmt.metrics(phases=True):
            for p in phases:
                for e in ["canon"] + Fmt.scen_effs("k"):
                    styles += style(f"{m}_{p}_{e}", 2 if m in ["vc", "aat"] else 1)

        return styles
