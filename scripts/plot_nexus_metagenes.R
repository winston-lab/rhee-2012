#!/usr/bin/env Rscript
library(argparse)
library(psych)
library(tidyverse)
library(forcats)
library(viridis)
library(ggthemes)

parser = ArgumentParser()
parser$add_argument('--inplus', dest='intable_plus', type='character')
parser$add_argument('--inminus', dest='intable_minus', type='character')
parser$add_argument('--inprotection', dest='intable_protection', type='character')
parser$add_argument('-s', dest='samplelist', type='character', nargs='+')
parser$add_argument('-t', dest='type', type='character')
parser$add_argument('-u', dest='upstream', type='integer')
parser$add_argument('-d', dest='downstream', type='integer')
parser$add_argument('-p', dest='trim_pct', type='double')
parser$add_argument('-r', dest='refptlabel', type='character', nargs='+')
parser$add_argument('-f', dest='factor', type='character', nargs='+')
parser$add_argument('-l', dest='scaled_length', type='integer')
parser$add_argument('-e', dest='endlabel', type='character', nargs='+')
parser$add_argument('-y', dest='ylabel', type='character', nargs='+')
parser$add_argument('--out1', dest='smeta_group_out', type='character')
parser$add_argument('--out2', dest='smeta_sample_out', type='character')
parser$add_argument('--out3', dest='pmeta_group_out', type='character')
parser$add_argument('--out4', dest='pmeta_sample_out', type='character')
parser$add_argument('--out5', dest='pmeta_goverlay_out', type='character')
parser$add_argument('--out6', dest='pmeta_soverlay_out', type='character')
parser$add_argument('--out7', dest='pmeta_soverlay_bygroup_out', type='character')
parser$add_argument('--out8', dest='meta_heatmap_group_out', type='character')
parser$add_argument('--out9', dest='meta_heatmap_sample_out', type='character')

args = parser$parse_args()

import = function(path){
    read_tsv(path,
    	 col_names=c("group", "sample", "index", "position","cpm"),
    	 col_types=cols(group=col_character(), sample=col_character(), index=col_integer(), position=col_double(), cpm=col_double())) %>%
    filter(cpm != "NA") %>%
    return()
}

format_xaxis = function(refptlabel, upstream, dnstream){
    function(x){
        if (first(upstream)>500 | first(dnstream)>500){
            return(if_else(x==0, refptlabel, as.character(x)))    
        }    
        else {
            return(if_else(x==0, refptlabel, as.character(x*1000)))
        }
    }
}

theme_default = theme_light() +
    theme(strip.text = element_text(size=12, face="bold", color="black"),
          strip.text.y = element_text(angle=0),
          axis.text.x = element_text(size=10, face="bold", color="black"),
          axis.text.y = element_text(size=10, color="black"),
          axis.title = element_text(size=10, face="plain"),
          panel.grid.major.x = element_line(color="grey60"),
          panel.grid.minor.x = element_line(color="grey80"),
          panel.grid.major.y = element_line(color="grey80"),
          panel.grid.minor.y = element_line(color="grey80"),
          panel.spacing.x = unit(0.5, "cm"),
          strip.background=element_blank(),
          plot.title = element_text(size=12, face="bold"),
          plot.subtitle = element_text(size=10))

main = function(intable_plus, intable_minus, intable_protection, samplelist, type, upstream, dnstream,
                trim_pct, refptlabel, factor, scaled_length, endlabel, ylabel, smeta_group_out, smeta_sample_out,
                pmeta_group_out, pmeta_sample_out, pmeta_goverlay_out,
                pmeta_soverlay_out, pmeta_soverlay_bygroup_out, meta_heatmap_group_out, meta_heatmap_sample_out){
    x_label = function(ggp){
        if(type=="absolute"){
            ggp = ggp +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   labels=format_xaxis(refptlabel=refptlabel,
                                                       upstream=upstream,
                                                       dnstream=dnstream),
                                   name=paste("distance from", refptlabel,
                                              if_else(upstream>500 | dnstream>500, "(kb)", "(nt)")),
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
        }
        else {
            ggp = ggp +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
            
        }
        return(ggp)
    }

    stranded_meta = function(df){
        metagene_base = ggplot(data=df, aes(x=position)) +
            geom_vline(xintercept=0, size=1, color="grey65")
        if(type=="scaled"){
            metagene_base = metagene_base +
                geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
        }
        metagene_base = metagene_base +
            geom_col(aes(y=pos.mean), fill="#114477", color="#114477", alpha=.9, size=0.1) +
            geom_col(aes(y=neg.mean), fill="#4477AA", color="#4477AA", alpha=.9, size=0.1) +
            ggtitle(paste("mean", factor, "ChIP-nexus coverage"),
                    subtitle=paste(nindices, ylabel)) +
            ylab("normalized counts") +
            theme_default
        return(x_label(metagene_base))
    }

    protection_meta = function(df){
        metagene_base = ggplot(data = df, aes(x=position)) +
            geom_vline(xintercept = 0, size=1, color="grey65")
        if(type=="scaled"){
            metagene_base = metagene_base +
                geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
        }
        metagene_base = metagene_base +
            geom_ribbon(aes(ymax=mean+1.96*sem, ymin=mean-1.96*sem),
                        fill="#114477", alpha=0.4, size=0) +
            geom_line(aes(y=mean), color="#114477", alpha=0.9) +
            scale_y_continuous(limits=c(0, NA), name="normalized counts") +
            ggtitle(paste("mean", factor, "protection"),
                    subtitle = paste(nindices, ylabel)) +
            theme_default
        return(x_label(metagene_base))
    } 

    meta_heatmap = function(df){
        metagene_base = ggplot(data = df, aes(x=position, y=0, fill=mean)) +
            geom_raster() +
            scale_fill_viridis(option="inferno", name="normalized\ncounts",
                               guide=guide_colorbar(barheight=8, barwidth=1)) +
            scale_y_continuous(breaks=0,expand=c(0,0), name=NULL) +
            ggtitle(paste("mean", factor, "protection"),
                    subtitle = paste(nindices, ylabel)) +
            theme_default +
            theme(axis.text.y = element_blank(),
                  strip.text.y = element_text(hjust=1, angle=-180),
                  legend.title = element_text(size=10, face="plain"),
                  axis.ticks.x = element_line(color="black", size=1),
                  axis.ticks.y = element_blank(),
                  panel.border = element_blank())
        return(x_label(metagene_base))
    }

    stranded = intable_plus %>% import() %>%
        right_join(intable_minus %>% import,
                   by=c("group","sample","index","position")) %>% 
        filter(sample %in% samplelist) %>% 
        mutate(cpm.y = -cpm.y) %>%
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
    
    repl_df = stranded %>% select(group, sample) %>% distinct() %>%
        group_by(group) %>% mutate(replicate=row_number()) %>% ungroup() %>% 
        select(-group)
    
    nindices = max(stranded$index, na.rm=TRUE)
    nsamples = length(fct_unique(stranded$sample))
    ngroups = length(fct_unique(stranded$group))
    
    sdf_group = stranded %>% group_by(group, position) %>%
        summarise(pos.mean = winsor.mean(cpm.x, trim=trim_pct, na.rm=TRUE),
                  neg.mean = winsor.mean(cpm.y, trim=trim_pct, na.rm=TRUE))
    smeta_group = stranded_meta(df=sdf_group) +
        facet_grid(~group)
    ggsave(smeta_group_out, plot=smeta_group, height=8, width=7*ngroups, units="cm")
    rm(smeta_group, sdf_group)
    
    sdf_sample = stranded %>% left_join(repl_df, by="sample") %>%
        group_by(group, sample, replicate, position) %>%
        summarise(pos.mean = winsor.mean(cpm.x, trim=trim_pct, na.rm=TRUE),
                  neg.mean = winsor.mean(cpm.y, trim=trim_pct, na.rm=TRUE))
    smeta_sample = stranded_meta(df=sdf_sample) +
        facet_grid(replicate~group) +
        theme(strip.text.y=element_text(angle=0))
    ggsave(smeta_sample_out, plot=smeta_sample, height=2+4.5*(max(repl_df$replicate)),
           width=7*ngroups, units="cm")
    rm(smeta_sample, sdf_sample, stranded)
    
    protection = import(intable_protection) %>%
        filter(sample %in% samplelist) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
    
    pdf_group = protection %>% group_by(group, position) %>%
        summarise(mean = winsor.mean(cpm, trim=trim_pct, na.rm=TRUE),
                  sem = winsor.sd(cpm, trim=trim_pct, na.rm=TRUE)/sqrt(n()))
    pmeta_group = protection_meta(df=pdf_group) +
        facet_grid(.~group)
    ggsave(pmeta_group_out, plot=pmeta_group, height=8, width=7*ngroups, units="cm")
    rm(pmeta_group)

    meta_heatmap_group = meta_heatmap(pdf_group) +
        facet_grid(group~., switch="y")
    ggsave(meta_heatmap_group_out, meta_heatmap_group,
           height=2.5+1.5*ngroups, width=14, units="cm")
    rm(meta_heatmap_group)
    
    pdf_sample = protection %>% left_join(repl_df, by="sample") %>% 
        group_by(group, sample, replicate, position) %>%
        summarise(mean = winsor.mean(cpm, trim=trim_pct, na.rm=TRUE),
                  sem = winsor.sd(cpm, trim=trim_pct, na.rm=TRUE)/sqrt(n()))
    pmeta_sample = protection_meta(df=pdf_sample) +
        facet_grid(replicate~group)
    ggsave(pmeta_sample_out, plot=pmeta_sample, height=2+4.5*(max(repl_df$replicate)),
           width=7*ngroups, units="cm")
    rm(pmeta_sample, protection)

    meta_heatmap_sample = meta_heatmap(pdf_sample) +
        facet_grid(sample~., switch="y")
    ggsave(meta_heatmap_sample_out, plot=meta_heatmap_sample,
           height=2.5+1.25*nsamples, width=14, units="cm")
    rm(meta_heatmap_sample)
    
    pmeta_goverlay = ggplot(data=pdf_group, aes(x=position, color=group, fill=group)) +
        geom_vline(xintercept = 0, size=1, color="grey65") 
    if(type=="scaled"){
        pmeta_goverlay = pmeta_goverlay +
            geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
    }
    pmeta_goverlay = pmeta_goverlay +
        geom_ribbon(aes(ymax=mean+1.96*sem, ymin=mean-1.96*sem),
                        alpha=0.3, size=0) +
        geom_line(aes(y=mean)) +
        scale_fill_ptol() + scale_color_ptol() +
        scale_y_continuous(limits=c(0, NA),
                           name="normalized coverage") +
        ggtitle(paste("mean", factor, "protection"),
                subtitle=paste(nindices, ylabel)) +
        theme_default +
        theme(legend.title = element_blank(),
              legend.text = element_text(size=12, face="bold"))
    pmeta_goverlay = x_label(pmeta_goverlay)
    ggsave(pmeta_goverlay_out, plot=pmeta_goverlay, width=14, height=8, units="cm")
    rm(pmeta_goverlay, pdf_group)
    
    pmeta_soverlay = ggplot(data=pdf_sample,
                            aes(x=position, group=sample, color=group, fill=group)) +
        geom_vline(xintercept = 0, size=1, color="grey65") 
    if(type=="scaled"){
        pmeta_soverlay = pmeta_soverlay +
            geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
    }
    pmeta_soverlay = pmeta_soverlay +
        geom_ribbon(aes(ymax=mean+1.96*sem, ymin=mean-1.96*sem),
                        alpha=0.2, size=0) +
        geom_line(aes(y=mean)) +
        scale_fill_ptol() + scale_color_ptol() +
        scale_y_continuous(limits=c(0, NA),
                           name="normalized coverage") +
        ggtitle(paste("mean", factor, "protection"),
                subtitle=paste(nindices, ylabel)) +
        theme_default +
        theme(legend.title = element_blank(),
              legend.text = element_text(size=12, face="bold"))
    pmeta_soverlay = x_label(pmeta_soverlay)
    ggsave(pmeta_soverlay_out, plot=pmeta_soverlay, width=14, height=8, units="cm")
    rm(pdf_sample)
    
    pmeta_soverlay_bygroup = pmeta_soverlay +
        facet_grid(group~., switch="y") +
        theme(legend.position="none",
              strip.placement="outside",
              strip.text.y=element_text(angle=180, hjust=1))
    ggsave(pmeta_soverlay_bygroup_out, plot=pmeta_soverlay_bygroup,
           width=14, height=2+4.5*(ngroups), units="cm")
}

main(intable_plus = args$intable_plus,
     intable_minus = args$intable_minus,
     intable_protection = args$intable_protection,
     samplelist = args$samplelist,
     type = args$type,
     upstream = args$upstream,
     dnstream= args$downstream,
     trim_pct = args$trim_pct,
     refptlabel = paste(args$refptlabel, collapse=" "),
     factor = paste(args$factor, collapse=" "),
     scaled_length = args$scaled_length,
     endlabel = paste(args$endlabel, collapse=" "),
     ylabel = paste(args$ylabel, collapse=" "),
     smeta_group_out = args$smeta_group_out,
     smeta_sample_out = args$smeta_sample_out,
     pmeta_group_out = args$pmeta_group_out,
     pmeta_sample_out = args$pmeta_sample_out,
     pmeta_goverlay_out = args$pmeta_goverlay_out,
     pmeta_soverlay_out = args$pmeta_soverlay_out,
     pmeta_soverlay_bygroup_out = args$pmeta_soverlay_bygroup_out,
     meta_heatmap_group_out = args$meta_heatmap_group_out,
     meta_heatmap_sample_out = args$meta_heatmap_sample_out)
