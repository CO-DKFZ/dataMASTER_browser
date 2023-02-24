
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

library(shiny)
library(shinydashboard)
library(GetoptLong)
library(DT)
library(igvShiny)
library(gtrellis)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(latex2exp)


load("dataMASTER.RData")
load("gencode19_gns_lite.RData")
GENCODE = gencode19_gns_lite


DB = dataMASTER
META = colData(DB)
EXPR = names(experiments(DB))

COMPONENT = list()

foo = function() {
source("src/config.R")
source("src/component.R")
source("src/comp_snv.R")
source("src/comp_snv_germline.R")
source("src/comp_indel.R")
source("src/comp_indel_germline.R")
source("src/utils.R")

ui = dashboardPage(
	dashboardHeader(title = "dataMASTER browser"),
	dashboardSidebar(
		sidebarMenu(
			id = "sidebar_menu",
			menuItem("Overview", tabName = "overview", icon = icon("landmark")),
			menuItem("Experiments", icon = icon("flask"),
				.list = c(list(menuSubItem("overview", tabName = "experiment")), lapply(EXPR, function(e) {
					menuSubItem(e, tabName = e)
				}))
			),
			menuItem("Gene", tabName = "gene", icon = icon("leaf")),
			menuItem("Patient", tabName = "patient", icon = icon("bed")),
			menuItem("Cohort", tabName = "cohort", icon = icon("people-group")),
			menuItem("Signatures", tabName = "signature", icon = icon("tower-broadcast"))

		)
    ),
	dashboardBody(
		includeCSS("www/master.css"),
		includeScript("www/master.js"),
		tabItems(
			.list = lapply(names(COMPONENT), function(e) {
				tabItem(tabName = e,
					COMPONENT[[e]]$ui
				)
			})
		)
	),
)

server = function(input, output, session) {

	lapply(names(COMPONENT), function(e) {
		COMPONENT[[e]]$server(input, output, session)
	})
}

print(shinyApp(ui, server))
}