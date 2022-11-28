## differential analysis ----
for_wilcox <- function(data_1, data_2, group,by ="group",dir = getwd()) {
  data_1 <- as.data.frame(t(data_1))
  data_2 <- as.data.frame(t(data_2))

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
  
  merge_data_1 <- subset(merge_data_all,merge_data_all$group == levels(merge_data_all$group)[1])
  
  merge_data_2 <- subset(merge_data_all,merge_data_all$group == levels(merge_data_all$group)[2])
  
  
  result_all <- data.frame(id = 1:(length(merge_data_all)-length(colnames(group))))
  
  for(i in col_num:length(merge_data_all)){
    a <- wilcox.test(merge_data_all[,i] ~ group,data = merge_data_all)
    result_all$gene[i-length(colnames(group))] <- colnames(merge_data_all)[i]
    result_all$pvalue[i-length(colnames(group))] <- a[["p.value"]]
    result_all$mean_data_1[i-length(colnames(group))] <- mean(merge_data_1[,i])
    result_all$mean_data_2[i-length(colnames(group))] <- mean(merge_data_2[,i])
    result_all$logFC[i-length(colnames(group))] <- result_all$mean_data_1[i-length(colnames(group))]-result_all$mean_data_2[i-length(colnames(group))]
    print(i/length(merge_data_all))
    }
  
  nrow(subset(result_all,abs(result_all$logFC) > sd(result_all$logFC)*3))
  result <- list(all = result_all,
                 fill_pvalue = subset(result_all,
                                      result_all$pvalue < 0.05),
                 up = subset(result_all,
                             result_all$pvalue < 0.05&result_all$logFC > 0),
                 down = subset(result_all,
                             result_all$pvalue < 0.05&result_all$logFC < 0))
  
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
  write.csv(result$all,file = paste0(dir,'result_all.csv'))
  write.csv(result$fill_pvalue,file = paste0(dir,'result.csv'))
  return(result)
}

