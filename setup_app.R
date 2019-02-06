library(shiny)
library(magrittr)
library(data.table)
library(GenomicRanges)
library(seqsetvis)
library(BiocFileCache)
source("module_scatter_plot.R")
bfc = BiocFileCache("~/.cache_shiny_ChIPview")
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