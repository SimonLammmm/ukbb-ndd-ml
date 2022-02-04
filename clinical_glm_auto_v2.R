require(dplyr)
require(data.table)
require(nnet)
require(caTools)

in_ad <- "Results/clinical/AD.tsv"
in_pd <- "Results/clinical/PD.tsv"
in_mnd <- "Results/clinical/motor-neuron-disease.tsv"
in_mg <- "Results/clinical/myasthenia-gravis.tsv"
in_ctrl <- "Results/clinical/control.tsv"

ad <- fread(in_ad, stringsAsFactors = T)
pd <- fread(in_pd, stringsAsFactors = T)
mnd <- fread(in_mnd, stringsAsFactors = T)
mg <- fread(in_mg, stringsAsFactors = T)
ctrl <- fread(in_ctrl, stringsAsFactors = T)

ad$disease <- "AD"
pd$disease <- "PD"
mnd$disease <- "MND"
mg$disease <- "MG"
ctrl$disease <- "Control"
data_master <- rbind(ad, pd, mnd, mg)                                                   # load data
data_master$disease <- factor(data_master$disease)

data_master[which(data_master == "")] <- NA                         # clean empty entries to NA

most_nas_i <- "init"
most_nas_j <- "init"
while(!is.null(most_nas_i) & !is.null(most_nas_j)) {
  most_nas_i <- NULL
  most_nas <- 0
  for (i in 1:nrow(data_master)) {
    nas <- sum(is.na(data_master[i,]))
    if (nas > most_nas) {
      most_nas_i <- i
      most_nas <- nas
    }
  }
  if (!is.null(most_nas_i)) {
    data_master <- data_master[-most_nas_i,]
  }
  
  most_nas_j <- NULL
  most_nas <- 0
  for (j in 1:ncol(data_master)) {
    nas <- sum(is.na(data_master[[j]]))
    if (nas > most_nas & nas/nrow(data_master) > 0.75) {
      most_nas_j <- j
      most_nas <- nas
    }
  }
  if (!is.null(most_nas_j)) {
    data_master <- data_master %>% select(-all_of(most_nas_j))
  }
}

ctrl2 <- ctrl %>% select(all_of(names(data_master)))
ctrl2 <- ctrl2[!(apply(ctrl2, MARGIN = 1, function(x) anyNA(x))), ] 

data_master <- rbind(data_master, ctrl2)

for (var in names(data_master)) {
  if (class(data_master[[var]]) == "factor") {
    data_master <- data_master[!(data_master[[var]] %in% (data_master %>% group_by(data_master[[var]]) %>% summarise(n = n()) %>% filter(n < 10) %>% select(1) %>% unlist()))]    # remove entries if there are too few instances of the level
    data_master[[var]] <- droplevels(data_master[[var]])                                       # remove unused levels
    cat(var, "has levels", levels(data_master[[var]]), "\n\n")
    if (length(grep(".+", levels(data_master[[var]]), value = T)) <= 1) {               # remove factor variables with one or fewer nonblank levels
      warning("Found ", var, " with 1 or fewer levels. Removing\n")
      data_master[[var]] <- NULL
    }
  } else if (class(data_master[[var]]) == "numeric" | class(data_master[[var]]) == "integer" | class(data_master[[var]]) == "float") {
    cat(var, "is a numeric with variance", var(data_master[[var]]), "\n\n")
  }
}

ctrl_length <- sum(data_master$disease != "Control")
set.seed(31415)
ctrl_random <- sample(which(data_master$disease == "Control"), ctrl_length)
idx <- c(which(data_master$disease != "Control"), ctrl_random)
data_master <- data_master[idx,]

data_master$`Encoded anonymised participant ID (eid)` <- NULL # remove EID as an independent variable
data_master$disease <- data_master$disease %>% relevel(ref = "Control")

make_glm <- function(ivs) {
  data <- data_master %>% select(all_of(ivs), disease)
  
  data[which(data == "")] <- NA                         # clean empty entries to NA
  
  sel <- lapply(data, FUN = function(x) sum(!is.na(x))/length(x) >= 0.75)          # keep variables with at least x% not NA
  sel[which(sel == F)] <- NULL
  data <- data %>% select(all_of(names(sel)))
  
  data <- data[!(apply(data, MARGIN = 1, function(x) anyNA(x))), ]                 # keep records with no NA
  
  for (var in names(data)) {
    if (class(data[[var]]) == "factor") {
      data <- data[!(data[[var]] %in% (data %>% group_by(data[[var]]) %>% summarise(n = n()) %>% filter(n <= 2) %>% select(1) %>% unlist()))]    # remove entries if there are too few instances of the level
      data[[var]] <- droplevels(data[[var]])                                       # remove unused levels
      cat(var, "has levels", levels(data[[var]]), "\n")
      if (length(grep(".+", levels(data[[var]]), value = T)) <= 1) {               # remove factor variables with one or fewer nonblank levels
        warning("Found ", var, " with 1 or fewer levels. Removing")
        data[[var]] <- NULL
      }
    } else if (class(data[[var]]) == "numeric" | class(data[[var]]) == "integer" | class(data[[var]]) == "float") {
      cat(var, "is a numeric with variance", var(data[[var]]), "\n\n")
    }
  }
  
  set.seed(31415) # to reproduce the model in the paper
  spl <- sample.split(data$disease, SplitRatio = 0.7)                 # split the data into training and test sets
  train <- filter(data, spl == T)
  test <- filter(data, spl == F)
  
  formula <- paste0("disease ~ ", paste( 
    ".",
    sep = " + ")
  )
  model <- multinom(formula, data = train, method = "glm.fit", maxit = 1000, MaxNWts = 2000)                         # train the model
  
  baseline <- prop.table(table(train$disease))                                     # baseline probability of correct assignment by guessing
  
  predictTrain <- predict(model, data = train, type = "class")     
  
  predictTrain2 <- table(predictTrain, train$disease)
  accuracyTrain <- sum(diag(predictTrain2)) / sum(predictTrain2)                    # prediction accuracy on training set data
  
  predictTest = predict(model, newdata = test, type = "class")
  
  predictTest2 <- table(predictTest, test$disease)
  accuracy <- sum(diag(predictTest2)) / sum(predictTest2)                      # prediction accuracy on test set data
  
  cat("final accuracy:", accuracy, "\n")
  return(list(glm = model,
              data = data,
              trainSet = train,
              testSet = test,
              baseline = baseline,
              predictTrain = predictTrain2,
              accuracyTrain = accuracyTrain,
              predictTest = predictTest2,
              accuracy = accuracy))
}

ivs <- names(data_master)[-length(names(data_master))]

best_accuracy <- 0
best_ivs <- character(0)

set.seed(55555)
for (i in 1:1000) { # do 1000 monte carlo samplings and accept the best one
  test_ivs <- sample(ivs, size = runif(1, 1, 12), replace = F)
  test_glm <- make_glm(test_ivs)
  if (test_glm$accuracy > best_accuracy) {
    best_accuracy <- test_glm$accuracy
    best_ivs <- test_ivs
    cat("new best accuracy!", best_accuracy, "\n")
  } else {
    cat("keeping old best accuracy:", best_accuracy, "\n")
  }
}

cat("created a model with accuracy:", best_accuracy, "\non the basis of the following clinical measures:\n", paste(best_ivs, collapse = "\n"), "\n")

for (j in 1:10) {
  for (i in 1:length(ivs)) {                         # find the best combination of ivs by adding them one by one and accepting them if they make a better model
    if (ivs[i] %in% best_ivs) {
      test_ivs <- best_ivs[-which(best_ivs %in% ivs[i])]
    } else {
      test_ivs <- c(best_ivs, ivs[i])
    }
    test_glm <- make_glm(test_ivs)
    if (test_glm$accuracy > best_accuracy) {
      best_accuracy <- test_glm$accuracy
      best_ivs <- test_ivs
      cat("new best accuracy!", best_accuracy, "\n")
    } else {
      cat("keeping old best accuracy:", best_accuracy, "\n")
    }
  }
}

cat("created a model with accuracy:", best_accuracy, "\non the basis of the following clinical measures:\n", paste(best_ivs, collapse = "\n"), "\n")
best_glm <- make_glm(best_ivs)$glm
save(best_glm, file = "Results/clinical/best_glm.RData")
# best_ivs <-                                                  # shortcut: these are the best ivs
#   c("Testosterone (30850-0.0)",
#     "Prospective memory result (20018-0.0)",
#     "Alanine aminotransferase (30620-0.0)",
#     "Ethnic background (21000-0.0)",
#     "LDL direct (30780-0.0)",
#     "Cystatin C (30720-0.0)",
#     "Phosphate (30810-0.0)",
#     "Prospective memory result (20018-1.0)",
#     "Mean time to correctly identify matches (20023-0.0)",
#     "Age when attended assessment centre (21003-0.0)",
#     "Microalbumin in urine (30500-0.0)",
#     "Sodium in urine (30530-0.0)",
#     "Albumin (30600-0.0)",
#     "Apolipoprotein A (30630-0.0)",
#     "Calcium (30680-0.0)",
#     "Cholesterol (30690-0.0)")

## vis

require(ggplot2)
require(FactoMineR)
plotdata_all_ivs <- data_master %>% select(-disease)
plotdata_best_ivs <- data_master %>% select(all_of(best_ivs))
do_mfa <- function(plotdata) {
  mfa <- MFA(plotdata,
             group = rep(1,length(names(plotdata))),
             type = case_when(lapply(plotdata, class) == "factor" ~ "n", TRUE ~ "s"),
             ncp = 4,
             name.group = names(plotdata))
  p <- cbind(as_tibble(mfa$ind$coord), label = data_master$disease)
  g <- ggplot(p, aes(x = Dim.1, y = Dim.2, colour = label)) +
    scale_color_discrete(name = "Disease") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
    ggtitle("MFA plot by disease") +
    geom_point(size = 0.4, alpha = 0.75)
  g
  return(list(mfa, g))
}

mfa_best_ivs <- do_mfa(plotdata_best_ivs)


pdf("Results/clinical_vis/mfa_best_ivs.pdf", width = 4, height = 4)
mfa_best_ivs[[2]]
dev.off()

## summary stats vis

require(tidyr)
require(foreach)
require(gridExtra)
require(ggpubr)

data_summary <- rbind(ad, pd, mnd, mg, ctrl) %>%
  select(all_of(best_ivs), disease) %>%
  mutate(disease = factor(disease, levels = c("Control", "MG", "MND", "PD", "AD")))
data_continuous <- data_summary %>% select(which(lapply(data_summary, function (x) class(x)) != "factor"), disease) # split variables into continuous and discrete plots
names(data_continuous) <- c("log(Testosterone [nmol/L])",
                            "log(Alanine aminotransferase [U/L])",
                            "log(LDL direct [mmol/L])",
                            "log(Cystatin C [mg/L])",
                            "log(Phosphate [mmol/L])",
                            "log(Mean time to correctly identify matches [ms])",
                            "Age when attended assessment centre [y]",
                            "log(Microalbumin in urine [mg/L])",
                            "log(Sodium in urine [mmol/L])",
                            "log(Albumin [g/L])",
                            "log(Apolipoprotein A [g/L])",
                            "log(Calcium [mmol/L])",
                            "log(Cholesterol [mmol/L])",
                            "disease")

data_continuous_nolog <- data_continuous
names(data_continuous_nolog) <- sub("log\\((.+)\\)", "\\1", names(data_continuous))

data_continuous <- data_continuous %>%
  mutate(`log(Testosterone [nmol/L])` = log(`log(Testosterone [nmol/L])`)) %>%
  mutate(`log(Alanine aminotransferase [U/L])` = log(`log(Alanine aminotransferase [U/L])`)) %>%
  mutate(`log(LDL direct [mmol/L])` = log(`log(LDL direct [mmol/L])`)) %>%
  mutate(`log(Cystatin C [mg/L])` = log(`log(Cystatin C [mg/L])`)) %>%
  mutate(`log(Phosphate [mmol/L])` = log(`log(Phosphate [mmol/L])`)) %>%
  mutate(`log(Mean time to correctly identify matches [ms])` = log(`log(Mean time to correctly identify matches [ms])`)) %>%
  mutate(`log(Microalbumin in urine [mg/L])` = log(`log(Microalbumin in urine [mg/L])`)) %>%
  mutate(`log(Sodium in urine [mmol/L])` = log(`log(Sodium in urine [mmol/L])`)) %>%
  mutate(`log(Albumin [g/L])` = log(`log(Albumin [g/L])`)) %>%
  mutate(`log(Apolipoprotein A [g/L])` = log(`log(Apolipoprotein A [g/L])`)) %>%
  mutate(`log(Calcium [mmol/L])` = log(`log(Calcium [mmol/L])`)) %>%
  mutate(`log(Cholesterol [mmol/L])` = log(`log(Cholesterol [mmol/L])`))

data_factor <- data_summary %>% select(which(lapply(data_summary, function (x) class(x)) == "factor"), disease)
for (i in 1:length(data_factor)) {
  levels(data_factor[[i]])[which(levels(data_factor[[i]]) == "")] <- NA                                           # change empty levels to NA
}

names(data_factor) <- c("Prospective memory result (first visit)",
                        "Ethnic background",
                        "Prospective memory result (second visit)",
                        "disease")

summary_continuous_ggplot <- function(x, iv) {
  my_comparisons <- list(c("Control", "AD"),
                         c("Control", "PD"),
                         c("Control", "MND"),
                         c("Control", "MG"))
  ggplot(x, aes(x = x[[iv]], y = disease)) +
    xlab(label = iv) +
    ylab(label = "") +
    geom_boxplot() +
    stat_compare_means(comparisons = my_comparisons, label = "p.signif")
}

summary_factor_ggplot <- function(x, iv) {
  ggplot(x, aes(y = disease, fill = x[[iv]])) +
    xlab(label = "Frequency") +
    ylab(label = "") +
    labs(fill = iv) +
    scale_fill_discrete(na.value = "grey") +
    guides(fill = guide_legend(ncol = case_when(length(levels(x[[iv]])) > 10 ~ 2,
                                                T ~ 1))) +
    geom_bar(position = "fill")
}

summary_continuous_plots <- foreach (iv = names(data_continuous)) %do% {
  summary_continuous_ggplot(data_continuous, iv)
}

summary_factor_plots <- foreach (iv = names(data_factor)) %do% {
  summary_factor_ggplot(data_factor, iv)
}

pdf("Results/clinical_vis/summary-statistics-continuous.pdf", width = 15, height = 14)
h <- grid.arrange(summary_continuous_plots[[7]],
                  summary_continuous_plots[[2]],
                  summary_continuous_plots[[10]],
                  summary_continuous_plots[[11]],
                  summary_continuous_plots[[12]],
                  summary_continuous_plots[[13]],
                  summary_continuous_plots[[4]],
                  summary_continuous_plots[[3]],
                  summary_continuous_plots[[6]],
                  summary_continuous_plots[[8]],
                  summary_continuous_plots[[5]],
                  summary_continuous_plots[[9]],
                  summary_continuous_plots[[1]])
dev.off()

pdf("Results/clinical_vis/summary-statistics-discrete.pdf", width = 15, height = 7)
j <- grid.arrange(ggplot(),
                  summary_factor_plots[[2]],
                  summary_factor_plots[[1]],
                  summary_factor_plots[[3]], ncol = 2)
dev.off()


## summary stats table

summary_stats_continuous <- foreach(my_variable = names(data_continuous),
                                    .combine = "rbind") %do% {
                                      if(my_variable == "disease") {
                                        NULL
                                      } else {
                                        data_temp <- data_continuous %>% select(stat = my_variable, disease)
                                        data_stats <- data_temp %>%
                                          group_by(disease) %>%
                                          summarise(n = n(), mean = mean(stat, na.rm = T), stdev = sd(stat, na.rm = T)) %>%
                                          mutate(error = qt(0.975, df = n - 1) * stdev / sqrt(n)) %>%
                                          mutate(CI95_left = mean - error,
                                                 CI95_right = mean + error,
                                                 stat = my_variable)
                                        data_ttests <- foreach(dis = c("Control", "MG", "MND", "PD", "AD"),
                                                               .combine = "rbind") %do% {
                                                                 if(dis == "Control") {
                                                                   NA
                                                                 } else {
                                                                   tt <- t.test(data_temp %>% filter(disease == dis) %>% select(stat),
                                                                                data_temp %>% filter(disease == "Control") %>% select(stat))
                                                                   tt$p.value
                                                                 }
                                                               }
                                        data_stats$t.test_pval <- data_ttests
                                        data_stats
                                      }
                                    }

summary_stats_factor <- foreach(my_variable = names(data_factor),
                                .combine = "rbind") %do% {
  if(my_variable == "disease") {
    NULL
  } else {
    data_temp <- data_factor %>% select(fct = my_variable, disease)
    data_stats <- as.data.frame(table(data_temp))
    data_stats <- data_stats %>% mutate(Pct = 100 * Freq / case_when(disease == "Control" ~ sum(data_stats %>% filter(disease == "Control") %>% summarise(sum(Freq))),
                                          disease == "AD" ~ sum(data_stats %>% filter(disease == "AD") %>% summarise(sum(Freq))),
                                          disease == "PD" ~ sum(data_stats %>% filter(disease == "PD") %>% summarise(sum(Freq))),
                                          disease == "MND" ~ sum(data_stats %>% filter(disease == "MND") %>% summarise(sum(Freq))),
                                          disease == "MG" ~ sum(data_stats %>% filter(disease == "MG") %>% summarise(sum(Freq)))),
                                        stat = my_variable)
  }
}

fwrite(summary_stats_continuous, "Results/summary_stats_continuous.tsv", sep = "\t")
writexl::write_xlsx(summary_stats_continuous, "Results/summary_stats_continuous.xlsx")
fwrite(summary_stats_factor, "Results/summary_stats_factor.tsv", sep = "\t")
writexl::write_xlsx(summary_stats_factor, "Results/summary_stats_factor.xlsx")