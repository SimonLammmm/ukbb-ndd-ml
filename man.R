require(qqman)
require(plyr)
require(dplyr)
require(data.table)

in_dict <- "Data/Axiom_UKB_WCSG.na34.annot.csv/Axiom_UKB_WCSG.na34.annot.csv"
in_ad <- "Results/GWAS/ad.tsv"
in_mnd <- "Results/GWAS/motor-neuron-disease.tsv"
in_mg <- "Results/GWAS/myasthenia-gravis.tsv"
in_pd <- "Results/GWAS/pd.tsv"

out_data_ad <- "Results/GWAS/ad_SNPs.tsv"
out_data_mnd <- "Results/GWAS/motor-neuron-disease_SNPs.tsv"
out_data_mg <- "Results/GWAS/myasthenia-gravis_SNPs.tsv"
out_data_pd <- "Results/GWAS/pd_SNPs.tsv"

out_plot_ad <- "Results/GWAS/ad_plot.png"
out_plot_mnd <- "Results/GWAS/motor-neuron-disease_plot.png"
out_plot_mg <- "Results/GWAS/myasthenia-gravis_plot.png"
out_plot_pd <- "Results/GWAS/pd_plot.png"

out_qq_ad <- "Results/GWAS/ad_qq.png"
out_qq_mnd <- "Results/GWAS/motor-neuron-disease_qq.png"
out_qq_mg <- "Results/GWAS/myasthenia-gravis_qq.png"
out_qq_pd <- "Results/GWAS/pd_qq.png"

dict <- fread(in_dict)
data_ad <- fread(in_ad)
data_mnd <- fread(in_mnd)
data_mg <- fread(in_mg)
data_pd <- fread(in_pd)

clean_for_manhattan <- function(data) {
    data$`dbSNP RS ID` <- NULL
    data$Chromosome <- NULL
    data$`Physical Position` <- NULL
    
    data <- data %>% left_join(dict, by = c("probeset_id" = "Probe Set ID"))
    data <- data %>% filter(!is.na(Chromosome)) %>%
        filter(Chromosome != "---") %>%
        filter(Chromosome != "") %>%
        filter(!is.na(test.statistic))
    data <- data %>% filter(!is.na(`dbSNP RS ID`))
    
    return(data)
    
}

plotdata_ad <- clean_for_manhattan(data_ad)
plotdata_mnd <- clean_for_manhattan(data_mnd)
plotdata_mg <- clean_for_manhattan(data_mg)
plotdata_pd <- clean_for_manhattan(data_pd)

fwrite(plotdata_ad, out_data_ad, sep = "\t")
fwrite(plotdata_mnd, out_data_mnd, sep = "\t")
fwrite(plotdata_mg, out_data_mg, sep = "\t")
fwrite(plotdata_pd, out_data_pd, sep = "\t")

plot_manhattan <- function(plotdata) {
    plotdata <- tibble(SNP = plotdata$`dbSNP RS ID`,
                       CHR = plotdata$Chromosome,
                       BP = as.numeric(plotdata$`Physical Position`),
                       P = as.numeric(plotdata$asymptotic.pvalue))
    
    unique(plotdata)
    smallest_nonzero_P <- plotdata %>% filter(P != 0) %>% select(P) %>% min()
    ylim <- round_any(-log10(smallest_nonzero_P), 10, ceiling)
    plotdata <- plotdata %>% mutate(P = case_when(P != 0 ~ P,
                                                  P == 0 ~ smallest_nonzero_P),
                                    CHR = case_when(CHR == "X" ~ 98,
                                                    #CHR == "Y" ~ 99,
                                                    T ~ as.numeric(CHR)))
    
    manhattan(plotdata, ylim = c(0,ylim), suggestiveline = F, genomewideline = F, chrlabs = c(1:22, "X"))
    
}

png(out_plot_ad, width = 15000, height = 3000, res = 600)
print({ plot_manhattan(plotdata_ad) })
dev.off()

png(out_plot_mnd, width = 15000, height = 3000, res = 600)
print({ plot_manhattan(plotdata_mnd) })
dev.off()

png(out_plot_mg, width = 15000, height = 3000, res = 600)
print({ plot_manhattan(plotdata_mg) })
dev.off()

png(out_plot_pd, width = 15000, height = 3000, res = 600)
print({ plot_manhattan(plotdata_pd) })
dev.off()

png(out_qq_ad, width = 3000, height = 3000, res = 600)
qq(plotdata_ad$asymptotic.pvalue)
dev.off()

png(out_qq_mnd, width = 3000, height = 3000, res = 600)
qq(plotdata_mnd$asymptotic.pvalue)
dev.off()

png(out_qq_mg, width = 3000, height = 3000, res = 600)
qq(plotdata_mg$asymptotic.pvalue)
dev.off()

png(out_qq_pd, width = 3000, height = 3000, res = 600)
qq(plotdata_pd$asymptotic.pvalue)
dev.off()

