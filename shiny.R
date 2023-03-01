
## if you have some packages missing, use the following command

# setRepositories(ind = 1:4)
# install.packages(c("MultiAssayExperiment", "RaggedExperiment", "org.Hs.eg.db", "shiny", "shinydashboard",
# 	"GetoptLong", "DT", "devtools" "igvShiny", "gtrellis", "ComplexHeatmap", "InteractiveComplexHeatmap", "circlize"))
# devtools::install_github("paul-shannon/igvShiny")
# 

setRepositories(ind = 1:4)
for(pkg in c("MultiAssayExperiment", "RaggedExperiment", "org.Hs.eg.db", "shiny", "shinydashboard",
 		"GetoptLong", "DT", "devtools", "gtrellis", "ComplexHeatmap", "InteractiveComplexHeatmap", 
 		"circlize", "latex2exp")) {
	if(!require(pkg, character.only = TRUE)) {
		install.packages(pkg)
		library(pkg)
	}
}
if(!require(igvShiny)) {
	devtools::install_github("paul-shannon/igvShiny")
}
library(igvShiny)

library(MultiAssayExperiment)
library(RaggedExperiment)
library(org.Hs.eg.db)

library(GetoptLong)
library(DT)
library(igvShiny)
library(gtrellis)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(latex2exp)
library(RColorBrewer)
library(jsonlite)
library(cola)
library(shiny)
library(shinydashboard)
library(dashboardthemes)

load("dataMASTER.RData")
load("gencode19_gns_lite.RData")
GENCODE = gencode19_gns_lite

source("src/config.R")
source("src/utils.R")


DB = dataMASTER
META = colData(DB)

EXPERIMENTS = experiments(DB)
EXPERIMENT_NAME = names(EXPERIMENTS)
EXPERIMENTS = lapply(EXPERIMENTS, function(obj) {
	if(inherits(obj, c("RangedSummarizedExperiment", "RaggedExperiment"))) {
		obj[as.vector(seqnames(rowRanges(obj)) %in% CHROMOSOME)]
	} else {
		obj
	}
})
CD_LIST = lapply(EXPERIMENT_NAME, function(x) colData(DB[, , x]))
names(CD_LIST) = EXPERIMENT_NAME

CNV_NORMALIZED = readRDS("normalized_CNV.rds")
CNV_GERMLINE_NORMALIZED = readRDS("normalized_CNV_germline.rds")

COMPONENT = list()

foo = function() {
for(f in list.files(path = "src", pattern = "comp_", full.name = TRUE)) {
	source(f)
}


ui = dashboardPage(
	dashboardHeader(title = "dataMASTER browser"),
	dashboardSidebar(
		sidebarMenu(
			id = "sidebar_menu",
			menuItem("Overview", tabName = "overview"),
			menuItem("Experiments", startExpanded = TRUE,
				.list = lapply(EXPERIMENT_NAME, function(e) {
					menuSubItem(e, tabName = e)
				})
			),
			menuItem("Query regions", tabName = "region"),
			menuItem("Patient", tabName = "patient"),
			menuItem("Cohort", tabName = "cohort")
		)
    ),
	dashboardBody(
		includeCSS("www/master.css"),
		includeScript("www/master.js"),
		shinyDashboardThemes(
	    	theme = "onenote"
	    ),
		tabItems(
			.list = lapply(names(COMPONENT), function(e) {
				tabItem(tabName = e,
					COMPONENT[[e]]$ui
				)
			})
		),
		style = "min-width:",
		tags$script(HTML('
$("#sidebarItemExpanded").append("<hr><p style=\'color:black;padding:4px 15px;\'>Developed by <a style=\'color:#3c8dbc;\' href=\'https://github.com/jokergoo/\' target=\'_blank\'>Zuguang Gu</a>. Source code is available at <a style=\'color:#3c8dbc;\' href=\'https://github.com/CO-DKFZ/dataMASTER_browser\' target=\'_blank\'>GitHub</a>.</p>");
	'))
	),
	tags$head(tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/bootstrap-3-typeahead/4.0.2/bootstrap3-typeahead.min.js"))
)

server = function(input, output, session) {

	lapply(names(COMPONENT), function(e) {
		COMPONENT[[e]]$server(input, output, session)
	})
}

print(shinyApp(ui, server))
}

foo()
