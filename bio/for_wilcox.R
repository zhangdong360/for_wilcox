## differential analysis ----
for_wilcox <- function(data_1, data_2, group,by ="group",
                       wilcox_exact = NULL,use_adjust_p = T, 
                       p_adj_method = "fdr", logFC_cut_off = 0, dir = getwd()) {
  library(progress)
  data_1 <- as.data.frame(data_1)
  data_2 <- as.data.frame(data_2)

  data_1$id <- rownames(data_1)
  data_2$id <- rownames(data_2)
  
  group$group <- factor(group[,by])
  
  merge_data <- merge(data_1,data_2,by = 'id')
  
  rownames(merge_data) <- merge_data$id
  
  merge_data <- subset(merge_data,select = -c(id))
  
  merge_data <- as.data.frame(t(merge_data))
  
  merge_data$id <- rownames(merge_data)
  
  group$id <- rownames(group)
  
  col_num <- length(colnames(group))+1
  
  merge_data_all <- merge(group,merge_data,by = 'id')
  
  merge_data_1 <- subset(merge_data_all,
                         merge_data_all$group == levels(merge_data_all$group)[1])
  
  merge_data_2 <- subset(merge_data_all,
                         merge_data_all$group == levels(merge_data_all$group)[2])
  
  length_result <- length(merge_data_all)-length(colnames(group))
  
  result_all <- data.frame(id = 1:length_result,
                           gene = rep(0,length_result),
                           mean_data_1 = rep(0,length_result),
                           mean_data_2 = rep(0,length_result),
                           logFC = rep(0,length_result),
                           pvalue = rep(0,length_result),
                           adjp = rep(0,length_result)
                           )

  pb <- progress_bar$new(
    format = "  wilcox.test [:bar] :current/:total :percent eta: :eta",
    total = length(merge_data_all)-length(colnames(group)), 
    clear = FALSE, 
    width= 60)
  # wilcox.test
  for(i in col_num:length(merge_data_all)){
    a <- wilcox.test(merge_data_all[,i] ~ group,
                     data = merge_data_all,
                     exact = wilcox_exact)
    j <- i-length(colnames(group))
    result_all$gene[j] <- colnames(merge_data_all)[i]
    result_all$pvalue[j] <- a[["p.value"]]
    result_all$mean_data_1[j] <- mean(merge_data_1[,i])
    result_all$mean_data_2[j] <- mean(merge_data_2[,i])
    result_all$logFC[j] <- result_all$mean_data_1[j]-result_all$mean_data_2[j]
    pb$tick()
    Sys.sleep(col_num/length(merge_data_all))
  }
  result_all <- result_all[order(result_all$pvalue),]
  if (use_adjust_p == T){
    # p value adjust
    adjp <- p.adjust(result_all$pvalue,method = p_adj_method)
    result_all$adjp <- adjp
    result <- list(all = result_all,
                   filter_pvalue = subset(result_all,
                                          result_all$adjp < 0.05),
                   filter_pvalue_cutoff = subset(result_all,
                                                 result_all$adjp < 0.05&abs(result_all$logFC) > logFC_cut_off),
                   up = subset(result_all,
                               result_all$adjp < 0.05&result_all$logFC > 0),
                   down = subset(result_all,
                                 result_all$adjp < 0.05&result_all$logFC < 0),
                   data = list(data_1 = data_1,
                               data_2 = data_2),
                   group = group,
                   wilcox.test_exact = wilcox_exact,
                   use_adjust_p = T,
                   logFC_cutoff = logFC_cut_off,
                   p_adj_method = p_adj_method
    )
    cat(paste0("up:    ",
               nrow(subset(result_all,
                           result_all$logFC>0&result_all$adjp < 0.05)),
               "\n",
               "down:  ",
               nrow(subset(result_all,
                           result_all$logFC<0& result_all$adjp < 0.05)),
               "\n",
               "not:   ",
               nrow(subset(result_all,
                           result_all$adjp>0.05))))
  }
  if (use_adjust_p == F){
    result <- list(all = result_all,
                   filter_pvalue = subset(result_all,
                                          result_all$pvalue < 0.05),
                   filter_pvalue_cutoff = subset(result_all,
                                                 result_all$pvalue < 0.05&abs(result_all$logFC) > logFC_cut_off),
                   up = subset(result_all,
                               result_all$pvalue < 0.05&result_all$logFC > 0),
                   down = subset(result_all,
                                 result_all$pvalue < 0.05&result_all$logFC < 0),
                   data = list(data_1 = data_1,
                               data_2 = data_2),
                   group = group,
                   wilcox.test_exact = wilcox_exact,
                   logFC_cutoff = logFC_cut_off,
                   use_adjust_p = F
    )
    cat(paste0("up:    ",
               nrow(subset(result_all,
                           result_all$logFC>0&result_all$pvalue < 0.05)),
               "\n",
               "down:  ",
               nrow(subset(result_all,
                           result_all$logFC<0& result_all$pvalue < 0.05)),
               "\n",
               "not:   ",
               nrow(subset(result_all,
                           result_all$pvalue>0.05))))
  }
  write.csv(result$all,file = paste0(dir,'result_all.csv'))
  write.csv(result$fill_pvalue,file = paste0(dir,'result.csv'))
  save(result,file = paste0(dir,'result.rdata'))
  return(result)
}

