require(CATTexact)
require(dplyr)
require(data.table)
require(foreach)
require(doParallel)

cl <- detectCores() - 2
registerDoParallel(cl)

studies <- c("ad", "pd", "motor-neuron-disease", "multiple-sclerosis", "myasthenia-gravis")

for (study in studies) {
    
    infile <- paste0("Results/axiom/", study, ".txt")
    infile_control <- "Results/axiom/control.txt"
    outdir <- "Results/GWAS/"
    dir.create(outdir, F)
    outfile <- paste0(outdir, study, ".tsv")
    
    cases <- fread(infile, fill = T)
    cases <- cases %>%
        filter(n_AA > 0 & n_AB > 0 & n_BB > 0) %>%
        mutate(n_AA = as.numeric(n_AA), n_AB = as.numeric(n_AB), n_BB = as.numeric(n_BB))
    
    controls <- fread(infile_control)
    controls <- controls %>%
        filter(n_AA > 0 & n_AB > 0 & n_BB > 0)
    
    
    cases_vs_controls <- inner_join(cases, controls, by = "probeset_id", suffix = c("_cases", "_controls"))
    cases_vs_controls <- cases_vs_controls %>%
        select(grep("probeset_id", names(cases_vs_controls)),
               grep("affy_snp_id", names(cases_vs_controls)),
               grep("n_AA", names(cases_vs_controls)),
               grep("n_AB", names(cases_vs_controls)),
               grep("n_BB", names(cases_vs_controls)),
               grep("affy_snp_id_annot", names(cases_vs_controls)),
               grep("Chromosome", names(cases_vs_controls)),
               grep("Chromosome Start", names(cases_vs_controls)),
               grep("Strand", names(cases_vs_controls)),
               grep("dbSNP RS ID", names(cases_vs_controls)))
    cat(names(cases_vs_controls), "\n")
    cat("Determining p-values for", study, "\n")
    pvalues <- foreach(i = 1:nrow(cases_vs_controls), .combine = "rbind") %dopar% {
        d <- c(0,1,2)
        n <- c(cases_vs_controls$n_AA_cases[i] + cases_vs_controls$n_AA_controls[i],
               cases_vs_controls$n_AB_cases[i] + cases_vs_controls$n_AB_controls[i],
               cases_vs_controls$n_BB_cases[i] + cases_vs_controls$n_BB_controls[i])
        r <- c(cases_vs_controls$n_AA_cases[i],
               cases_vs_controls$n_AB_cases[i],
               cases_vs_controls$n_BB_cases[i])
        catt_exact(d, n, r)
    }

    cat("Finished determining p-values for", study, "\n")    
    cases_vs_controls <- cbind(cases_vs_controls, pvalues)
    
    cat("Writing to file:", study, "\n")
    fwrite(cases_vs_controls, outfile, sep = "\t")
    
}

stopImplicitCluster()

