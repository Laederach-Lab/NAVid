from ..lib import ui
from ..uix import Interface, Fmt

app_ui = ui.page_fluid(

    # ----- STYLE (CSS) ----- #

    # Results Card
    ui.tags.style(
        """
        .results-card-title {
            font-weight: 600;
            margin: 0px;
            padding: 0px;
        }

        .results-content {
            display: flex;
            width: 100%;
            margin: 0px;
            justify-content: center;
            align-items: center;
        }

        .results-single {
            display: flex;
            width: 100%;
            align-items: center;
            justify-content: center;

            &.normal {
                height: 100px;
            };

            &.tall {
                height: 135px;
            }
        }

        .results-compare {
            display: table-cell;
            width: 50%;
            margin: 0px;
            text-align:center;

            &.value {
                height: 90px;
                padding-bottom: 10px;
                vertical-align: bottom;
            };

            &.label {
                height: 45px;
                vertical-align: top;
                font-weight: 425;
                font-style: italic;
            }
        }

        .results-value {
            font-weight: bold;
            text-align: center;
            
            &.sub {
                font-weight: 500;
            }
        }

        .ids-div {
            display: flex;
            justify-content: center;
        }
        """
    ),

    # Frequency Dropdown
    ui.tags.style(
        """
        .freq-accordion-title {
            font-weight: 550;
            color: #00314e;
        }

        .freq-value-1 {
            font-weight: 550;
        }

        .freq-value-2 {
            line-height: 1.35;
        }
        """
    ),

    # NAG Strength
    ui.tags.style(
        """
        .strength-value-1 {
            font-weight: bold;
            justify-content: center;
        }

        .strength-value-2 {
            font-weight: 500;
            line-height: 1.35;
            justify-content: center;
        }
        """
    ),

    # Tables
    ui.tags.style(
        """
        .table th {
            font-weight: 600;
            text-align: left;
        }

        .table-values {
            border-spacing: 0px;
            padding: 0px 0px;

            &.numeric {
                text-align: right;
            }
        }

        .table-index {
            font-weight: 600;
            text-align: left !important;
            # white-space: pre-wrap;
        }

        .table-highlight {
            font-weight: bold;
            background-color: #ffff40 !important;
        }
        """
    ),

    ui.tags.style(Interface.table_styles()),

    # Text
    ui.tags.style(
        """
        .error-message {
            text-align: center;
            color: #EE0000;
            font-weight: 550
        }

        .about-title {
            font-weight: 500;
        }

        .toast.text-bg-primary {
            background-color: #f3cccc !important;
            color: #c10000 !important;
        }
        """
    ),

    # Input
    ui.tags.style(
        """
        .switch {
            display: flex;
            align-items: center;
            justify-content: right;
            padding-right: 8px;
            margin: 0;
            height: 100%;
            width: 100%;
        }

        .switch label {
            padding-left: 8px;
            padding-right: 0;
            margin-right: 0;
            width: 140px;
            font-weight: 550;
            font-size: 12;
        }

        .checkbox label {
            display: inline-block;
            padding-left: 30px;
            text-indent: -30px;
        }
        """
    ),

    # Layout
    ui.tags.style(
        """
        div {
            padding: 0px;
            margin: 0px
        }

        .page-div {
            margin-top: 12px;
        }

        .input-div {
            margin-top: 12px;
        }

        .download-div {
            padding: 12px 0px
        }

        .card-header {
            background-color: #cce5f3 !important;
            color: #00314e !important;
            padding: 0px;
        }

        .shiny-card-container {
            min-height: 500px !important
        }

        .shiny-input-container {
            margin: 4px 0px 0px 0px !important;
        }

        .shiny-input-container:not(.shiny-input-container-inline) {
            width: fit-content;
        }

        .input-group {
            width: 300px;
        }
        """
    ),

    # Styles
    ui.tags.style(
        """
        .button, .btn-file {
            --bs-btn-color: #00314e !important;
            --bs-btn-border-color: #00314e !important;
            --bs-btn-hover-bg: #00314e !important;
            --bs-btn-hover-border-color: #00314e !important;
        }

        input:checked {
            background-color: #00314e !important;
            border-color: #00314e !important;
        }

        # .progress-bar {
        #     background-color: #00314e !important;
        # }
        """
    ),

    # ----- CONTENT ----- #

    # Title
    ui.div(
        ui.h1("NAVid", style="color:#00314e; font-weight:650; line-height:1; margin:20px 0px 12px"),
        ui.h4("Identify and assess NAGNAG-affecting variants",
              style="color:#00314e; font-weight:500; line-height:1; margin-bottom: 16px"),
        style = "margin-bottom:16px; border-bottom:2px solid #00314e",
    ),

    # Tab Bar
    ui.navset_pill(
        
        # Variant Search
        ui.nav_panel(
            "Search",
            
            ui.div(
                ui.layout_sidebar(

                    # Input Sidebar
                    ui.sidebar(

                        # Examples ('example')
                        ui.div(
                            Interface.input_button(id="example", label="Get Example", style="padding:12px; font-weight:600; width:100%"),
                            class_ = "input-div"
                        ),

                        # Chromosome ('chrom')
                        ui.div(
                            Interface.input_text(id="chrom", label="Chromosome", example="chr1",
                                                 tooltip="Chromosome on which variant is located\nformat: chr#"),
                            class_ = "input-div"
                        ),

                        # Variant Position ('var_pos')
                        ui.div(
                            Interface.input_text(id="var_pos", label="Variant Location", example="100,000",
                                                 tooltip="Chromosome position of the variant"),
                            class_ = "input-div"
                        ),

                        # Alternate Allele ('alt_allele')
                        ui.div(
                            Interface.input_text(id="alt_allele", label="Alternate Allele", example="A, C, G, T",
                                                 tooltip="Base specific to the variant"),
                            class_ = "input-div"
                        ),
                        
                        # Search Button ('search')
                        ui.div(
                            Interface.input_button(id="search", label="Find Existing Variant", style="padding:8px; width:100%",
                                                   tooltip="Search dbSNP and ClinVar accession numbers"),
                            class_ = "input-div"
                        ),

                        # 5' Splice Site Position ('ss5_pos')
                        ui.div(
                            Interface.input_text(id="ss5_pos", label="5' Splice Site Location", example="100,000",
                                                 tooltip="Chromosome position of the first base at the 5' end of the intron"),
                            
                            class_ = "input-div"
                        ),

                        # Strand ('strand'), Switch ('switch_strand')
                        ui.div(
                            ui.div(
                                {"style":"display:table; padding:0px; margin:0px; height:fit-content; width:100%"},

                                # Strand ('strand')
                                ui.div(
                                    Interface.input_radio(id="strand", label="Strand", choices=["+", "-"], initial="+",
                                                          tooltip="Genomic strand on which the variant/transcript is located"),
                                    style = "display:table-cell; padding:0px; margin:0px; width:fit-content"
                                ),

                                # Switch ('switch_strand')
                                ui.div(
                                    Interface.input_button(id="switch_strand", label="Switch", style="padding:6px 8px", tooltip="Use opposite strand"),
                                    style = "display:table-cell; padding:0px; margin:0px; text-align:right; vertical-align:bottom"
                                )
                            ),
                            class_ = "input-div"
                        ),

                        # Index ('index'), Switch ('switch_index')
                        ui.div(
                            ui.div(
                                {"style":"display:table; padding:0px; margin:0px; height:fit-content; width:100%"},

                                # Index ('index')
                                ui.div(
                                    Interface.input_radio(id="index", label="Index", choices=["0", "1"], initial="1",
                                                          tooltip="Coordinate of the first base in the chromosome"),
                                    style = "display:table-cell; padding:0px; margin:0px; width:fit-content"
                                ),

                                # Switch ('switch_index')
                                ui.div(
                                    Interface.input_button(id="switch_index", label="Switch", style="padding:6px 8px", tooltip="Use other index"),
                                    style = "display:table-cell; padding:0px; margin:0px; text-align:right; vertical-align:bottom"
                                )
                            ),

                            class_ = "input-div"
                        ),

                        # Submit Button ('submit')
                        ui.div(
                            Interface.input_button(id="submit", label="Submit", style="padding:12px; font-weight:650; width:100%"),
                            class_ = "input-div"
                        ),

                        # Reset Button ('reset')
                        ui.div(
                            Interface.input_button(id="reset", label="Reset", style="padding:8px 16px; width:fit-content"),
                            class_ = "input-div", style = "display: flex; width: 100%; justify-content: left"
                        ),
                        
                        position = "left",
                        gap = "0px",
                        padding = "42px 25px 25px",
                    ),

                    # Results
                    ui.output_ui("search_results_output"),
                    gap = "0px",
                    padding = "12px"
                ),

                class_ = "page-div"
            ),

            value = "search"
        ),

        # Variant Search Upload
        ui.nav_panel(
            "Upload",

            ui.layout_column_wrap(
                
                # Upload
                ui.div(
                    Interface.card_w_header(
                        ui.h4("Upload", class_="results-card-title"),

                        ui.div(
                            
                            # instructions
                            ui.p("Files should be tab-delimited with the first 4 columns as follows:", style="line-height:1; margin-bottom:12px"),
                            ui.output_text_verbatim("input_columns"),
                            ui.p("(see right for examples)", style="line-height:1; margin-top:12px; margin-bottom:6px"),
                            ui.p("Other columns are permissible but will be ignored.", style="line-height:1"),

                            # upload
                            ui.div(
                                ui.div(
                                    ui.input_file(
                                        id = "file_upload",
                                        label = "",
                                        multiple = True,
                                        accept = ["text/plain", ".tsv", ".csv", ".txt"]
                                    ),
                                    style = "display:table-cell; padding-right:16px; width:fit-content"
                                ),

                                # clear uploads
                                ui.div(
                                    Interface.input_button(id="clear_uploads", label="Clear", style="padding:6px 12px; font-weight:650"),
                                    style = "display:table-cell; width:100%"
                                ),

                                style = "display:table; margin:0px; padding:0px"
                            ),

                            ui.div(ui.output_ui("download_processed_files")),

                            # style = "height:fit-content; min-height:400px"
                        ),
                    ),
                    class_ = "page-div"
                ),

                # Example
                ui.div(
                    Interface.card_w_header(
                        ui.h4("Example", class_="results-card-title"),

                        ui.div(
                            # example input
                            ui.h5("Input", style="font-weight:600; padding-bottom:3px"),
                            ui.output_text_verbatim("example_input"),
                            
                            # example output
                            ui.h5("Output", style="font-weight:600; padding:12px 0px 3px"),
                            ui.output_text_verbatim("example_output"),
                        )
                    ),

                    class_ = "page-div"
                ),

                width = 0.5
            )
        ),

        # Spacer
        Interface.navset_pill_divider(),

        # Variant Effects
        ui.nav_panel(
            "Variant Effects",
            
            ui.div(
                ui.navset_pill_list(
                    Interface.effects_data_panel("effects_summary", "Summary"),
                    
                    # divider
                    ui.nav_control(
                        ui.div({"style":"margin:15px 0px; border-top:1px solid lightgray"})
                    ),

                    Interface.effects_data_panel("create_motifs", "NAGNAG-Creating Motifs"),
                    Interface.effects_data_panel("alter_motifs", "NAGNAG-Altering Motifs"),
                    Interface.effects_data_panel("destroy_motifs", "NAGNAG-Destroying Motifs"),

                    # download button
                    ui.nav_control(
                        Interface.download_button(id="download_variants_tar_gz", label="Download All"),
                    ),

                    widths = (2, 10),
                    id = "variant_effects_navset"
                ),

                class_ = "page-div list"
            ),

            value = "variant_effects"
        ),

        # Proteome Effects
        ui.nav_panel(
            "Proteome Effects",

            ui.div(
                ui.layout_sidebar(

                    # Sidebar
                    ui.sidebar(
                        Interface.input_radio(id="tsn_data_freq", label="Frequency",
                                            choices={"prop":"Proportion", "count":"Counts"},
                                            initial="prop", inline=False),
                        Interface.download_button(id="download_proteome_tar_gz", label="Download All"),

                        width = "225px",
                        gap = "12px",
                        padding = "25px"
                    ),

                    # Metric Tabs
                    ui.navset_pill(
                        Interface.tsn_metric_data_panel("vc"),
                        Interface.tsn_metric_data_panel("aat"),
                        Interface.tsn_metric_data_panel("aatt"),
                        Interface.tsn_metric_data_panel("ins"),
                        Interface.tsn_metric_data_panel("del"),
                        id = "proteome_effects_navset"
                    ),

                    gap = "0px",
                    padding = "12px",
                ),

                class_ = "page-div"
            ),

            value = "proteome_effects"
        ),

        # Observed NAGNAGs
        ui.nav_panel(
            "Observed NAGNAGs",
            
            ui.div(
                ui.layout_sidebar(

                    # Sidebar
                    ui.sidebar(
                        
                        # Type
                        Interface.input_radio(id="site_type_to_show", label="Observed", choices={"sites":"Sites", "scenarios":"Scenarios"},
                                            initial="sites", inline=False),

                        # Effects
                        Interface.input_check(id="effects_to_show", label="Variant Effects",
                                            choices=Fmt.scen_effs("d"), initial=Fmt.scen_effs("d")),
                        
                        # Download
                        Interface.download_button(id="download_observed_tar_gz", label="Download All"),

                        gap = "12px",
                        padding = "25px"
                    ),

                    # Database Tabs
                    ui.output_ui("observed_nagnags_panel"),

                    gap = "0px",
                    padding = "12px",
                ),

                class_ = "page-div"
            ),
            
            value = "observed_nagnags"
        ),

        # Spacer
        Interface.navset_pill_divider(),

        # About
        ui.nav_panel(
            "About",

            ui.div(
                ui.layout_columns(
                    Interface.about_card(
                        "Usage",
                        """<strong>NAVid</strong> allows users to identify genomic variants that create, alter, or destroy a NAGNAG 3' splice site and assess the proteomic consequences of that NAGNAG.<br>
                        All coordinates used are based on human reference genome GRCh38/hg38. Other genomes and assemblies are not currently supported."""
                    ),

                    Interface.about_card(
                        "Authors",
                        """<strong>NAVid</strong> was developed using <a href="https://shiny.posit.co/py/">Shiny for Python</a> by Halle Pearce as part of the <a href="https://ribosnitch.bio.unc.edu/">Laederach Lab</a> at the University of North Carolina at Chapel Hill, Department of Biology."""
                    ),
                    Interface.about_card(
                        "Contact",
                        """If you have questions or run into issues, please contact us at <a href='mailto:hpearce@unc.edu'>hpearce@unc.edu</a> or <a href='mailto:alain@unc.edu'>alain@unc.edu</a>."""
                    ),
                    
                    col_widths = (12, 6, 6),
                    gap = "18px",
                ),

                class_ = "page-div"
            ),

            value = "about"
        ),

        id = "app_navset"
    ),

    title = "NAVid",
    lang = "en",
    zoom = "100%"
)
