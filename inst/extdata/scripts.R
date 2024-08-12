library(usethis)
library(devtools)
library(roxygen2)
library(BiocCheck)
load_all();document();roxygenize()
library(scQTLtools)
check()
check(vignettes = FALSE)
devtools::build()
BiocCheck::BiocCheck('scQTLtools' = TRUE)
BiocCheck::BiocCheckGitClone()
readCitationFile("inst/CITATION")
data(testEQTL)
Gene <- rownames(testEQTL@rawData$normExpMat)
SNP <- rownames(testEQTL@rawData$snpMat)
linearResult <- linearModel(
   eQTLObject = testEQTL,
   geneIDs = Gene,
   snpIDs = SNP,
   biClassify = FALSE,
   pAdjustMethod = "bonferroni",
   pAdjustThreshold = 0.05,
   logfcThreshold = 0.025
 )

data(testEQTL)
eqtl <- callQTL(eQTLObject = testEQTL, useModel = "linear")
visualizeQTL(eQTLObject = eqtl,
             SNPid = "1:632647",
             Geneid = "RPS27",
             groupName = NULL,
             plottype = "violin",
             removeoutlier = TRUE)

 rm(list = ls())

data(testEQTL)
Gene <- rownames(testEQTL@filterData$expMat)
SNP <- rownames(testEQTL@filterData$snpMat)
eQTLObject <- testEQTL
geneIDs = Gene
snpIDs = SNP
biClassify = FALSE
pAdjustMethod = "bonferroni"
pAdjustThreshold = 0.05
logfcThreshold = 0.02

result <- data.frame(
  SNPid = character(),
  group = character(),
  Geneid = character(),
  pvalue = double(),
  adjusted_pvalue = double(),
  b = double(),
  abs_b = double(),
  Remark = character(),
  stringsAsFactors = FALSE
)
expressionMatrix <- round(eQTLObject@filterData$expMat * 1000)
snpMatrix <- eQTLObject@filterData$snpMat
unique_group <- unique(eQTLObject@groupBy$group)

result_all <- data.frame()

for (k in unique_group) {
  k <- "CMP"
  split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == k]
  expressionMatrix_split <- expressionMatrix[, split_cells]
  snpMatrix_split <- snpMatrix[, split_cells]

  message(k, ":")
  message("0%   10   20   30   40   50   60   70   80   90   100%")
  message("[----|----|----|----|----|----|----|----|----|----|")
  pb <- progress_bar$new(
    total = length(snpIDs),
    format = "[:bar]",
    clear = FALSE,
    width = 51
  )
  for (i in seq_len(length(snpIDs))) {
    snpid <- snpIDs[1]
    snp_mat <- snpMatrix_split[snpid, ]
    snp_mat <- as.data.frame(snp_mat)
    snp_mat$cells <- rownames(snp_mat)


    for (j in seq_len(length(geneIDs))) {
      gene_id <- geneIDs[47]
      gene_mat <- expressionMatrix_split[gene_id, ]
      gene_mat <- as.data.frame(gene_mat)
      gene_mat$cells <- rownames(gene_mat)

      combined_df <- merge(snp_mat, gene_mat, by = "cells")
      combined_df <- subset(combined_df, snp_mat != 0)

      lmodel <- stats::glm(combined_df$gene_mat ~ combined_df$snp_mat,
                           family = poisson()
      );summary(lmodel)$coefficients

      if(length(summary(lmodel)$coefficients[, "Pr(>|z|)"]) >= 2) {
        lmout_pvalue <- summary(lmodel)$coefficients[2, "Pr(>|z|)"]
        lmout_b <- summary(lmodel)$coefficients[2, "Estimate"]
        new_row <- data.frame(
          SNPid = snpid,
          group = k,
          Geneid = geneIDs[47],
          pvalue = lmout_pvalue,
          b = lmout_b
        )
        result <- rbind(result, new_row)
      }

      new_row <- data.frame(
        SNPid = snpid,
        group = k,
        Geneid = geneIDs[j],
        pvalue = lmout_pvalue,
        b = lmout_b
      )
      result <- rbind(result, new_row)
    }
    pb$tick()
  }
  message("finished!")

  if (!pAdjustMethod %in% c("bonferroni",
                            "holm",
                            "hochberg",
                            "hommel",
                            "BH")) {
    stop("Invalid p-adjusted method.
         Please choose from 'bonferroni', 'holm', 'hochberg', 'hommel',
           or'fdr or BH'.")
  }

  # adjust p-value
  result[, "adjusted_pvalue"] <- p.adjust(result[, "pvalue"], method = "BH")
  result <- result[order(result[, "adjusted_pvalue"]), ]
  rownames(result) <- NULL
  result <- result[result$adjusted_pvalue <= pAdjustThreshold, ]

  # abs_b
  result <- result %>%
    mutate(abs_b = abs(result[, "b"]))

  result <- result[result$abs_b >= logfcThreshold, ]

  result_all <- rbind(result_all, result)
}
return(result_all)


#' visualizeQTL: Visualize the gene-snp pairs by group.
#' @param SNPid ID of SNP.
#' @param Geneid ID of Gene.
#' @param plottype Types of plot,one of "QTLplot","violin","boxplot" or
#' "histplot".
#' @param eQTLObject An S4 object of class eQTLObject.
#' @param groupName Users can choose one or more than one single cell groups.
#' @param removeoutlier Whether identify and remove the outliers.
#' Default by FALSE.
#'
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_annotation
#' @importFrom stats median
#' @importFrom graphics title
#'
#' @return list
#' @export
#' @examples
#' data(testEQTL)
#' ## We have to call the eQTLs firstly using `callQTL()`.
#' eqtl <- callQTL(eQTLObject = testEQTL, useModel = "linear")
#' visualizeQTL(eQTLObject = eqtl,
#' SNPid = "1:632647",
#' Geneid = "RPS27",
#' groupName = NULL,
#' plottype = "QTLplot",
#' removeoutlier = FALSE)
visualizeQTL <- function(eQTLObject,
                         SNPid,
                         Geneid,
                         groupName = NULL,
                         plottype = "QTLplot",
                         removeoutlier = FALSE) {
  eQTLresult <- eQTLObject@eQTLResult
  expressionMatrix <- eQTLObject@filterData$expMat
  snpMatrix <- eQTLObject@filterData$snpMat
  biClassify <- eQTLObject@biClassify

  if (is.null(groupName)) {
    unique_group <- unique(eQTLObject@groupBy$group)
  } else {
    unique_group <- groupName
  }

  df_all <- data.frame()

  for (i in unique_group) {
    split_cells <- rownames(eQTLObject@groupBy)[eQTLObject@groupBy$group == i]
    split_expressionMatrix <- expressionMatrix[, split_cells]
    snpMatrix_split <- snpMatrix[, split_cells]

    result <- eQTLresult[(eQTLresult$group == i) &
                           (eQTLresult$SNPid == SNPid) &
                           (eQTLresult$Geneid == Geneid), ]

    if (biClassify == TRUE) {
      result_split1 <- result$group == i
      result_split <- result[result_split1, ]

      adjust.pvalue <- result_split$adjusted_pvalue

      snpMatrix_split[snpMatrix_split == 3] <- 2

      ref_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid, ] == 1]
      alt_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid, ] == 2]

      if (removeoutlier) {
        remove_outliers(split_expressionMatrix,
                        Geneid,
                        ref_cells,
                        alt_cells)
      }

      counts_Ref <- unlist(split_expressionMatrix[Geneid, ref_cells])
      counts_Alt <- unlist(split_expressionMatrix[Geneid, alt_cells])

      df <- data.frame(
        expression = c(counts_Ref, counts_Alt),
        snp = c(
          rep("REF", length(counts_Ref)),
          rep("ALT", length(counts_Alt))), group = i)
      df$snp <- factor(df$snp, levels = c("REF", "ALT"))
      title <- paste(i)

      df_all <- rbind(df_all, df)
    } else if (biClassify == FALSE) {
      result_split1 <- result$group == i
      result_split <- result[result_split1, ]

      adjust.pvalue <- result_split$adjusted_pvalue

      AA_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid, ] == 1]
      Aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid, ] == 3]
      aa_cells <- colnames(snpMatrix_split)[snpMatrix_split[SNPid, ] == 2]

      if (removeoutlier) {
        remove_outliers(split_expressionMatrix,
                        Geneid,
                        AA_cells,
                        Aa_cells,
                        aa_cells)
      }

      counts_AA <- unlist(split_expressionMatrix[Geneid, AA_cells])
      counts_Aa <- unlist(split_expressionMatrix[Geneid, Aa_cells])
      counts_aa <- unlist(split_expressionMatrix[Geneid, aa_cells])

      df <- data.frame(
        expression = c(counts_AA, counts_Aa, counts_aa),
        snp = c(rep("ref/ref", length(counts_AA)),
                rep("ref/alt", length(counts_Aa)),
                rep("alt/alt", length(counts_aa))),group = i)
      df$snp <- factor(df$snp, levels = c("ref/ref", "ref/alt", "alt/alt"))
      snp <- df$snp
      title <- paste(i)

      df_all <- rbind(df_all, df)
    } else {
      stop("biClassify can only be selected as 'TRUE' or 'FALSE'")
    }
  }

  title_all <- paste("Plot", "of", Geneid, "and", SNPid)

  options(warn = -1)
  plot_list <- list()
  if (plottype == "QTLplot") {
    for (j in seq_len(length(unique_group))) {
      group_split <- unique_group[j]
      title <- group_split
      df_split <- df_all[df_all$group == group_split, ]
      plot <- draw_QTLplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  } else if (plottype == "violin") {
    for (j in seq_len(length(unique_group))) {
      group_split <- unique_group[j]
      title <- group_split
      df_split <- df_all[df_all$group == group_split, ]
      plot <- draw_violinplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  } else if (plottype == "boxplot") {
    for (j in seq_len(length(unique_group))) {
      group_split <- unique_group[j]
      title <- group_split
      df_split <- df_all[df_all$group == group_split, ]
      plot <- draw_boxplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  } else if (plottype == "histplot") {
    for (j in seq_len(length(unique_group))) {
      group_split <- unique_group[j]
      title <- group_split
      df_split <- df_all[df_all$group == group_split, ]
      plot <- draw_histplot(df_split, unique_group[j])
      plot_list[[j]] <- plot
    }
  } else {
    stop("Invalid plottype,
         Please choose from 'QTLplot', 'violin' , 'boxplot' or 'histplot'.")
  }

  combined_plot <- wrap_plots(plot_list)
  title_annotation <- plot_annotation(
    title = title_all,
    theme = theme(
      plot.title = element_text(
        size = 16, hjust = 0.5, vjust = 1
      )
    )
  )
  combined_plot <- combined_plot + title_annotation
  message(combined_plot)
  options(warn = 0)
}



QTLplot <- function(df, unique_group){
  ggplot(data = df,aes(x = snp,
                       y = expression ,
                       fill = factor(snp)))+
    scale_fill_manual(values = c("#D7AA36", "#D85356", "#94BBAD")) +
    geom_violin(alpha = 0.7, position = position_dodge(width = .75),
                size = 0.8, color="black") +
    geom_boxplot(notch = FALSE, outlier.size = -1,
                 color="black", lwd=0.5, alpha = 0.7) +
    geom_point(shape = 21, size=1.8, stroke=NA,
               position = position_jitterdodge(jitter.width = 0.8),
               alpha = 1.2) +
    theme_bw() +
    labs(title = unique_group, y = "Expression", x = "") +
    theme(axis.text.x = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))
}


drawviolinplot <- function(df, unique_group){
  ggplot(data = df, aes(x = snp,
                        y = expression ,
                        fill = factor(snp)))+
    scale_fill_manual(values = c("#D7AA36", "#D85356", "#94BBAD")) +
    geom_violin(alpha = 0.7, position = position_dodge(width = .75),
                size = 0.8, color="black") +
    theme_bw() +
    labs(title = unique_group, y = "Expression", x = '') +
    theme(axis.text.x = element_text(size = 15, color = "black"),
          axis.text.y = element_text(size = 12, color = "black"),
          axis.ticks = element_line(size=0.2, color="black"),
          axis.ticks.length = unit(0.2, "cm"),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "none",
          panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 12))
}



drawhistplot <- function(df, unique_group){
  ggplot(df, aes(expression, fill = factor(snp)))+
    geom_histogram()+
    facet_grid(snp ~ ., margins=FALSE, scales="free_y") +
    scale_fill_brewer(palette = "Pastel1")+
    labs(title = unique_group,
         x = "Expression",
         y = "Count")+
    theme_minimal()+
    ggtitle(title)+
    theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_text(vjust = 0.5, face = "bold"),
          plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = "none")+
    guides(fill = guide_legend(title = NULL), color = guide_legend(title = NULL))
}
