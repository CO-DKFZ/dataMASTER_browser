
COMPONENT$region = list()

COMPONENT$region$ui = div(
	h3("Query by regions"),
	div(
		textInput("region_search", "Gene/Genomic region", placeholder = "TP53 or chr1:10000-20000"),
		textAreaInput("region_sample_list", "Patients/Samples (optional)"),
		actionButton("region_search_submit", "Submit")
	),
	br(),
	tabsetPanel(
		type = "tabs",
		tabPanel(qq("IGV browser"),
			uiOutput(qq("region_igv_ui")),
			style = "padding:16px 0px"
		),
		tabPanel(qq("Genes"),
			uiOutput(qq("region_genes_ui")),
			style = "padding:16px 0px"
		),
		tabPanel(qq("Table"),
			uiOutput(qq("region_table_ui")),
			style = "padding:16px 0px"
		)
	)
)

COMPONENT$region$server = function(input, output, session) {
	observeEvent(input$region_search_submit, {
		
		updateTabItems(session, "sidebar_menu", "region")

		region_search = input$region_search
		updateTextInput(session, "region_search", value = region_search)

		if(!grepl("^(chr)?[0-9XYxy]+:", region_search)) {
			gr_g = as.data.frame(GENCODE[GENCODE$gene_name == region_search])
			region_search = paste0("chr", gr_g[, 1], ":", gr_g[, 2], "-", gr_g[, 3])
			if(nrow(gr_g) == 0) {
				region_search = ""
			}
		}
		
		gr_search = parse_gr_str(region_search)

		if(length(gr_search) == 0 || nrow(gr_g) == 0) {
			output[["region_igv_ui"]] = renderUI({
				p("Region string should be a gene symbol or in a format of '1:1000-2000'.")
			})
			return(NULL)
		}

		samples = NULL
		if(!grepl("^\\s*$", input$region_sample_list)) {
			samples = parse_samples(input$region_sample_list)
			if(length(samples) == 0) {
				output[["region_igv_ui"]] = renderUI({
					p("No sample is found.")
				})
				return(NULL)
			}
		}

		showNotification("Loading IGV JavaScript browser", duration = 4, type = "message")

		gr_search = GRanges(seqnames = gr_search[[1]], ranges = IRanges(gr_search[[2]], gr_search[[3]]))

		experiments = c("snv", "snv_germline", "indel", "indel_germline", "sv", "sv_germline", "fusion")
		names(experiments) = experiments
		
		output[[qq("region_igv_ui")]] = renderUI({
				
			div(
				div(
					column(8, selectInput("region_igv_track", label = "Tracks", choices = experiments, selected = "snv")),
					column(2, actionButton("region_add_track_submit", "Add new track"), style="margin-top:24px"),
					div(style = "clear:both"),
					style = "width:500px;margin-left:-15px;"
				),
				igvShinyOutput("region_igv"),
			)
		})

		output[["region_igv"]] = renderIgvShiny({

			options = parseAndValidateGenomeSpec(genomeName=GENOME,  initialLocus=region_search,
	                                      stockGenome=TRUE, dataMode="stock",
	                                      fasta=NA, fastaIndex=NA, 
	                                      genomeAnnotation=NA)
			igvShiny(options)
			
		})

		ov = findOverlaps(GENCODE, gr_search)
		ghits = GENCODE[unique(queryHits(ov))]
		ghits = ghits[ghits$gene_type == "protein_coding"]

		output[[qq("region_genes_ui")]] = renderUI({
			tableOutput("region_genes")
		})

		output[["region_genes"]] = renderTable({
			
			columns = c("GENENAME", "ENTREZID", "ENSEMBL")
			tb = NULL
			for(g in ghits$gene_name) {
				tbb = select(org.Hs.eg.db, keys = g, keytype = "SYMBOL", columns = columns)
				tbb = lapply(tbb, function(x) paste(unique(x), collapse = ", "))
				tbb = as.data.frame(tbb)
				tb = rbind(tb, tbb)
			}

			tb$ENTREZID = qq("<a href='https://www.ncbi.nlm.nih.gov/gene/@{tb$ENTREZID}' target='_black'>@{tb$ENTREZID}</a>", collapse = FALSE)
			tb$ENSEMBL = qq("<a href='http://www.ensembl.org/id/@{tb$ENSEMBL}' target='_blank'>@{tb$ENSEMBL}</a>", collapse = FALSE)
			
			colnames(tb) = c("Symbol", "Name", "EntrezID", "EnsemblID")
			tb
		}, sanitize.text.function = function(x) x)


		glhits = lapply(experiments, function(e) {
			grl = DB[[e]]@assays
			if(!is.null(samples)) {
				grl = grl[intersect(names(grl), samples)]
			}
			grr = unlist(grl)
			ov = findOverlaps(grr, gr_search)
			grr[unique(queryHits(ov))]
		})

		output[[qq("region_table_ui")]] = renderUI({
			ll = sapply(glhits, length) > 0
			al = lapply(experiments[ll], function(e) {
				actionLink(qq("region_table_@{e}"), class="region_table_link", label = qq("@{length(glhits[[e]])} @{e}"), style="padding:4px 10px;color:#aaaaaa", onclick="$('.region_table_link').css('color', '#aaaaaa');$(this).css('color', '#3c8dbc');")
			})
			names(al) = NULL
			div(
				do.call(div, al),
				div(
					DTOutput("region_table_output"),
					style = "margin:10px 15px; padding: 10px; background-color:white;"
				)
			)
		})

		for(e in experiments) {
			local({
				e = e
				observeEvent(input[[qq("region_table_@{e}")]], {
					output[["region_table_output"]] = renderDT({
						gr = glhits[[e]]
						nm = names(gr)
						names(gr) = NULL
						df = as.data.frame(gr)
						df = cbind(Sample = nm, df)
						df
					}, server = TRUE, options = list(scrollX = TRUE))
				})
			})
		}


		observeEvent(input$region_add_track_submit, {
			e = input$region_igv_track

			if(e == "gencode v19") {
				loadGFF3TrackFromLocalData(session, id = "Gencode v19", trackName = e, tbl.gff3 = GENCODE_GFF)
			} else {

				if(is.null(samples)) {
					gr = unlist(DB[[e]]@assays)
				} else {
					gr = unlist(DB[[e]]@assays[samples])
				}
				# ov = findOverlaps(gr, gr_search)
				# gr = gr[unique(queryHits(ov))]
				gr_nm = names(gr)
				names(gr) = NULL
				tbl = as.data.frame(gr)[, 1:3]
				colnames(tbl) = c("chr", "start", "end")
				tbl$name = ""
				tbl$score = gr_nm
				tbl[, 1] = as.character(tbl[, 1])
				loadBedTrack(session, id = "region_igv", trackName = e, tbl = tbl)
			}	
		})
	})
}
