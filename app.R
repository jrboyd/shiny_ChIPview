library(shiny)
library(magrittr)
library(data.table)
library(GenomicRanges)
library(seqsetvis)
library(BiocFileCache)
source("module_scatter_plot.R")
bfc = BiocFileCache("~/.cache_shiny_MM_diffbind")
bfcif = function(bfc, rname, FUN, force_overwrite = FALSE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname")) == 0){
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)
        
    }else{
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        message("loading from cache...")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        message("running function and caching...")
        res = FUN()
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}

# path to diffbind result
db_file = "example_diffbind.txt"
data_dirs = dir(file.path("/slipstream/galaxy/uploads/working", 
                          "/qc_framework/output_AF_MCF10_CTCF/"), 
                pattern = "CTCF.+ed$", full.names = TRUE)
# load db_file to create db_gr
db_dt = fread(db_file)
colnames(db_dt)[1] = "id"
k = grepl("chr[0-9]{1,2}$", db_dt$Chr)
table(db_dt[k,]$Chr)
db_dt = db_dt[k,]
db_gr = GRanges(db_dt)

var_names = colnames(mcols(db_gr))
var_names = var_names[sapply(mcols(db_gr), is.numeric)]
var_names = var_names[var_names != "id"]

bw_files = sapply(data_dirs, function(d)dir(d, pattern = "FE.bw$", full.names = TRUE))
bw_files_dt = data.table(files = bw_files)
bw_files_dt[, c("cell", "mark") := tstrsplit(basename(files), "_", keep = 1:2)]
bam_files = sapply(data_dirs, function(d)dir(d, pattern = ".bam$", full.names = TRUE))
bam_files_dt = data.table(files = bam_files)
bam_files_dt[, c("cell", "mark") := tstrsplit(basename(files), "_", keep = 1:2)]

# prof_dt = bfcif(bfc, digest::digest(list(bw_files_dt, db_gr)),
#                 function(){
#                     ssvFetchBigwig(bw_files_dt, db_gr,
#                                 return_data.table = TRUE,
#                                 n_cores = 8)
#                 })

prof_dt = bfcif(bfc, digest::digest(list(bam_files_dt, db_gr)),
                function(){
                    ssvFetchBam(bam_files_dt, db_gr, 
                                return_data.table = TRUE, 
                                n_cores = 8, 
                                target_strand = "both")
                })
setkey(prof_dt, id)

#db_gr and prof_dt must be set

# Define UI for miles per gallon app ----
ui <- pageWithSidebar(
    
    # App title ----
    headerPanel("Interactive ChIP-seq"),
    
    # Sidebar panel for inputs ----
    sidebarPanel(
        selectInput("x_variable", label = "X", choices = var_names, selected = var_names[1]),
        selectInput("y_variable", label = "Y", choices = var_names, selected = var_names[2])#,
        # uiOutput("GroupingsAvailable")
    ),
    
    
    # Main panel for displaying outputs ----
    mainPanel(                   
        fixedRow(
            column(width = 3,
                   radioButtons("RadioScatterplotType", label = "Plot Type", choices = c("volcano", "standard"), selected = "volcano"),
                   # shinyjs::hidden(checkboxInput("CheckShowHelpers", label = "Show Plot Help", value = FALSE)),
                   checkboxInput("CheckFixCoords", label = "Fixed Coordinates", value = TRUE)
            ),
            column(width = 9,
                   plotOutput("xy_values", width = 640, height = 640,
                              click = "xy_click",
                              brush = brushOpts(
                                  id = "xy_brush"
                              )),
                   plotOutput("summaryPlot"),
                   plotOutput("detailPlot")
            )
        )
    )
)

# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
    # shinyjs::useShinyjs()
    features_gr = reactiveVal(db_gr)
    xy_dt = reactiveVal(db_dt)
    selected_dt = reactiveVal(data.table())
    
    set_selected = function(new_selected_dt){
        selected_dt(new_selected_dt)
    }
    
    
    server_scatter_plot(input, output, session, 
                        get_features = features_gr, 
                        get_profiles = xy_dt, 
                        set_selected = set_selected)
    
    observe({
        sel_dt = selected_dt()
        if(is.null(sel_dt)) return(NULL)
        if(nrow(sel_dt) == 0) return(NULL)
        showNotification(paste(as.character(nrow(sel_dt)), "selected"))    
    })
    
    output$summaryPlot = renderPlot({
        sel_dt = selected_dt()
        if(is.null(sel_dt)) return(NULL)
        if(nrow(sel_dt) == 0) return(NULL)
        num_detail = 5
        summary_dt = prof_dt[.(as.character(sel_dt$id))]
        agg_dt = summary_dt[, .(y = mean(y)), by = .(x, cell, mark, strand)]
        sp_agg_dt = applySpline(agg_dt, by_ = c("strand", "cell"), n = 10)
        ggplot(sp_agg_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_grid(".~cell")
    })
    
    output$detailPlot = renderPlot({
        sel_dt = selected_dt()
        if(is.null(sel_dt)) return(NULL)
        if(nrow(sel_dt) == 0) return(NULL)
        num_detail = 5
        summary_dt = prof_dt[.(as.character(sel_dt$id))]
        ids = unique(summary_dt$id)
        ids = sample(ids, min(num_detail, length(ids)))
        detail_dt = summary_dt[id %in% ids]
        sp_detail_dt = applySpline(detail_dt, by_ = c("strand", "cell", "id"), n = 10)
        ggplot(sp_detail_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_grid("id~cell")
    })
        
        
        
        
}

shinyApp(ui, server)