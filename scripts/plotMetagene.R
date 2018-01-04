library(tidyverse)
library(forcats)
#library(viridis)
library(psych)

import = function(path){
    read_tsv(path,
    	 col_names=c("factor", "assay", "sample", "index", "position","cpm"),
    	 col_types=cols(factor=col_character(), assay=col_character(),
    	                sample=col_character(), index=col_integer(),
    	                position=col_double(), cpm=col_double())) %>%
    filter(cpm != "NA") %>%
    return()
}

format_xaxis_kb = function(refptlabel){
    function(x) if_else(x==0, refptlabel, as.character(x))
}
format_xaxis_nt = function(refptlabel){
    function(x) if_else(x==0, refptlabel, as.character(x*1000))
}
label_xaxis = function(ggp, refptlabel, upstream, downstream){
    if(upstream>500 | downstream>500){
        ggp = ggp +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                               labels=format_xaxis_kb(refptlabel=refptlabel),
                               name=paste("distance from", refptlabel, "(kb)"),
                               limits = c(-upstream/1000, downstream/1000),
                               expand=c(0,0))
    } else{
        ggp = ggp +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                               labels=format_xaxis_nt(refptlabel=refptlabel),
                               name=paste("distance from", refptlabel, "(kb)"),
                               limits = c(-upstream/1000, downstream/1000),
                               expand=c(0,0))
    }
}

plotmeta = function(intable.plus, intable.minus, trim_pct, upstream, downstream,
                    refptlabel, factors, ylabel, outpath){
    raw = intable.plus %>% import() %>%
            right_join(intable.minus %>% import,
                       by=c("factor","assay","sample","index","position")) %>% 
            filter(factor==factors)
    raw$cpm.y = -raw$cpm.y
    raw$sample = fct_inorder(raw$sample, ordered = TRUE)
    raw$factor = fct_inorder(raw$factor, ordered = TRUE)
    nindices = max(raw$index, na.rm=TRUE)
    
    df = raw %>% group_by(factor, assay, position) %>%
            summarise(pos.mean = winsor.mean(cpm.x, trim=trim_pct, na.rm=TRUE),
                      neg.mean = winsor.mean(cpm.y, trim=trim_pct, na.rm=TRUE))
    
    metagene = ggplot(data = df, aes(x=position))
    
    if(refptlabel=="TATA"){
        metagene = metagene +
                    annotate(geom="rect", xmin=0,xmax=8e-3,
                             ymin=min(df$neg.mean)*1.05,
                             ymax=max(df$pos.mean)*1.05,
                             fill="grey65")
    }
    metagene = metagene +
                geom_col(aes(y=pos.mean), fill="#08306b", alpha=.90) +
                geom_col(aes(y=neg.mean), fill="#2171b5", alpha=.90) +
                theme_bw() +
                scale_y_continuous(expand=c(0,0),
                                   name="normalized counts") +
                ggtitle(paste("mean", factors, "coverage"),
                        subtitle = paste(nindices, ylabel)) +
                theme(plot.title = element_text(size=12, face="bold"),
                      plot.subtitle = element_text(size=12),
                      strip.text = element_text(size=12, face="bold"),
                      strip.background = element_blank(),
                      axis.text.x = element_text(size=12, face="bold", color="black"),
                      axis.text.y = element_text(size=10, color="black"),
                      axis.title = element_text(size=12, face="bold"),
                      panel.grid.major.x=element_line(color="grey60"),
                      panel.grid.minor.x=element_line(color="grey80"),
                      panel.grid.major.y=element_line(color="grey80"),
                      panel.grid.minor.y=element_line(color="grey80"),
                      panel.spacing.x= unit(.5, "cm")) +
                facet_grid(.~assay)
    metagene = metagene %>% label_xaxis(upstream=upstream, downstream=downstream, refptlabel=refptlabel)
    ggsave(outpath, metagene, width=11, height=8, units="cm")
}

plotmeta(intable.plus = snakemake@input[["plus"]],
         intable.minus = snakemake@input[["minus"]],
         trim_pct = snakemake@params[["trim_pct"]],
         upstream = snakemake@params[["upstream"]],
         downstream= snakemake@params[["downstream"]],
         refptlabel = snakemake@params[["refptlabel"]],
         factors = snakemake@params[["factor"]],
         ylabel = snakemake@params[["ylabel"]],
         outpath = snakemake@output[[1]])