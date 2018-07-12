##### Hotspot程式為算出hotspot的範圍並輸出圖片, 輸入檔案為處理過的mutation資訊, 會利用biomaRt讀取ensembl資料庫中的資料, 輸出檔案只有突變 #####

library(grid)
library(ggplot2)
library(ggrepel)
library(shape)
library("e1071")
library("plotrix")
#library("biomaRt")
library(plyr)

getwd()

source("hotspot_fun4.r")

##### 參數設定 #####
### FC1 刪減 hotspot 的第一條件 (移除較少突變數量的 hotspot) ###
total_prop = 0.05
### FC2 刪減 hotspot 的第二條件 (比較 hotspot 間的突變數量) ###
compare_hotspot = 0.2
### FC3 刪減 hotspot 的第三條件 (each hotspot mutation count / hotspot patients)  ###
patient_density <- 0.1
### FC4 刪減 hotspot 的第四條件 (each hotspot / hotspot length) ###
hotspot_density = 2

### Score 參數 (決定考慮的範圍) ###
alph1 = 2.5
### 山岳刪減分群法 參數(決定刪減分數的範圍) ###
beta = 0.005
### 延伸距離的最小值 ###
peak_region_min = 3
### 延伸距離的最大值 ###
peak_region_max = 10
### 剩餘多少 gap 資訊 ###
remain_gap = 0.2
### 山岳刪減分群法 參數(決定找多少個群中心) ###
num_center = 10

###Input###
input <- ("input_blca.txt")

###Output###
output_pdf <- ("blca.pdf")
file_output <- ("blca.txt")
filter_output <-("blca_filter.txt")

M_hotspot(input,output_pdf,file_output,filter_output,total_prop,compare_hotspot , patient_density ,hotspot_density,alph1 ,beta,peak_region_min,peak_region_max ,remain_gap ,num_center)

M_hotspot <- function(input,output_pdf,file_output,filter_output,total_prop,compare_hotspot , patient_density ,hotspot_density,alph1 ,beta,peak_region_min,peak_region_max ,remain_gap ,num_center){
  b <- read.delim(input, header = F)
  c <- read.delim("ensembl-81-final.txt", header = F)
  qax <- as.matrix(unique(b$V1))
  #Output
  pdf(output_pdf)
  par(mfrow = c(2, 1))
  #write.table(input_name,file_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")
  Parameters_name <- c("peak_region_min\tpeak_region_max\tremain_gap\talpha\tnum_center\tbeta\tFC1\tFC2\tFC3\tFC4")
  write.table(Parameters_name,file_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")
  parameter <- cbind(peak_region_min, peak_region_max, remain_gap, alph1, num_center, beta, total_prop, compare_hotspot,patient_density,hotspot_density)
  write.table(parameter, file_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")
  out_col_name <- c("Cancer\tSymbol\tENSP\tAmino_Acid_length\tTotal_mutation_count\tMaximum_muta_count_per\tpeak_region\tHotspot_mutation_count\tHotspot_mutation_region\tlength\t#position")
  write.table(out_col_name , file_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")  
  filter_name <- c("total_prop\tcompare_hotspot\thotspot_density\tpatient_density\talph1\tbeta\tpeak_region_min\tpeak_region_max\tremain_gap")
  write.table(filter_name,filter_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")
  Filter <- cbind(total_prop, compare_hotspot,hotspot_density,patient_density,alph1,beta,peak_region_min, peak_region_max, remain_gap)
  write.table(Filter,filter_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")
  qwe <- apply(qax, 1, Main_fun)
  
  dev.off()
  
}


  print(which(qax == x))
  ALL <- c[which(c$V3 ==x),]
  Length <- as.numeric(ALL[5])
  genes <-x
  symbol <- as.matrix(ALL[4])
  bb <- b[b$V1 ==x,]
  bb <- bb[(which((as.numeric(bb[, 2]) <= Length) == TRUE)),]
  patient <- unique(bb$V3)
  patient <- data.frame(patient,stringsAsFactors = F)
  #write.table(symbol,file_output1, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = FALSE, na = "")
  if (nrow(bb) != 0) {
    colnames(bb)[2] <- "position"
    NAME <-x
    final_P <- NULL
    final <- NULL
    hotspot_position <- NULL
    hotspot_region <- NULL
    tmp1 = NULL
    rCount <- each_position(Length, bb)
    colnames(rCount) <- c("position", "mutation")
    ##### mutation為有發生突變的位點 #####		
    mutation <- which(rCount[, 2] != 0)
    Total <- sum(rCount[, 2])
    name <- paste(NAME, symbol, sep = "_")
    Name <- paste(NAME, symbol, Total, sep = "_")
    new_rCount <- rCount[!rCount$mutation %in% 0,]
    ##### peak_region_min為延伸距離自訂的最小值, peak_region_max為延伸距離自訂的最大值, remain_gap為取前x%的 gap 資訊 #####
    peak_region <- Extend_region(peak_region_min, peak_region_max, remain_gap, mutation)
    
    ##### alph1為轉換成密度資訊時的參數 #####
    #final <- new_Density_information(Length, alph1, rCount)
    final <- new_Density_information(Length, alph1, new_rCount)
    ##### num_center為山岳刪減分群法想找到多少個群中心, beta為可調整之參數 #####
    final_P <- mountain_clusting(mutation, final, num_center, Length, beta)
    
    if(length(unique(bb[,2]))==1)
    {
      final_count<-t(as.matrix(c(nrow(bb),unique(bb$position))))
      final_count <- data.frame(final_count)
      fc1_value <- as.numeric(final_count[,1])/Total
      fc2_value <- as.numeric(final_count[,1])/max(as.numeric(final_count[,1]))
      fc3_value <- as.numeric(final_count[,1])/nrow(patient)
      fc4_value <- as.numeric(final_count[,1])
      fc1 <- as.numeric(final_count[,1])/Total > total_prop
      fc2 <- as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) > compare_hotspot
      fc3 <- as.numeric(final_count[,1])/nrow(patient) > patient_density
      fc4 <- as.numeric(final_count[,1]) > hotspot_density
      
      total <- Total
      gene <-  symbol[1]
      filter_result <- cbind(final_count,fc1_value,fc2_value,fc3_value,fc1,fc2,fc4,fc3,total,gene)
      #colnames(filter_result) <- c("mutation count","region","fc5_value","fc1","fc2","fc4","fc3")
      
      write.table(filter_result,filter_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = TRUE, na = "")
      final_region<-final_count[as.numeric(final_count[,1])/Total> total_prop & as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) > compare_hotspot & as.numeric(final_count[,1]) > hotspot_density & as.numeric(final_count[,1])/nrow(patient) > patient_density ,]
      if (length(final_region) == 2) {
        final_region <- (as.matrix(final_region))
      }
      
    }else{
      final_count<-extend(final_P,peak_region,mutation,rCount)
      fc1_value <- as.numeric(final_count[,1])/Total
      fc2_value <- as.numeric(final_count[,1])/max(as.numeric(final_count[,1]))
      fc3_value <- as.numeric(final_count[,1])/nrow(patient)
      fc1 <- as.numeric(final_count[,1])/Total > total_prop
      fc2 <- as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) > compare_hotspot
      fc3 <- as.numeric(final_count[,1])/nrow(patient) > patient_density
      fc4 <- as.numeric(final_count[,4]) > hotspot_density
      #fc5 <- as.numeric(final_count[,1])/all_p_num >= FC5
      filter_result <- cbind(final_count,fc1_value,fc2_value,fc3_value,fc1,fc2,fc4,fc3)
      colnames(filter_result) <- c("mutation count","region","hotspot rate","fc4_value","fc1_value","fc2_value","fc3_value","fc1","fc2","fc4","fc3")
      filter_result <- data.frame(filter_result)
      filter_result$Total <- "" 
      filter_result[,'Total'] <- Total
      filter_result$gene <- ""
      filter_result[,'gene'] <- symbol[1]
      #filter_result$cancer <- ""
      #filter_result[,'cancer'] <- bb$V5
      write.table(filter_result,filter_output, sep = "\t", quote = FALSE, row.name = FALSE, append = TRUE, col.name = TRUE, na = "")
      final_region_count<-final_count[ as.numeric(final_count[,1])/Total > total_prop & as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) > compare_hotspot & as.numeric(final_count[,4]) > hotspot_density & as.numeric(final_count[,1])/nrow(patient) > patient_density,1:4]
      final_region<-final_count[as.numeric(final_count[,1])/Total > total_prop & as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) > compare_hotspot & as.numeric(final_count[,4]) > hotspot_density & as.numeric(final_count[,1])/nrow(patient)> patient_density,1:2]
      if(length(final_region_count)==4){
        final_region_count <- t(as.matrix(final_region_count))
      }
      if (length(final_region) == 2) {
        final_region <- t(as.matrix(final_region))
      }
    }
    
    #if (length(final_region) == 2) {
    #  final_region <- (as.matrix(final_region))
    #}
    
    if(nrow(final_region) == 0){
      #new_final_region_count<-final_count[ as.numeric(final_count[,1]) & as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) & as.numeric(final_count[,3]) & as.numeric(final_count[,4]) ,1:4]
      new_final_region_count<-final_count[ as.numeric(final_count[,1]) & as.numeric(final_count[,1])/max(as.numeric(final_count[,1])) ,1:2]
      
      #new_final_region<-final_count[as.numeric(final_count[,1]) & as.numeric(final_count[,1])/max(as.numeric(final_count[,1]))  &  as.numeric(final_count[,3])  & as.numeric(final_count[,4]) ,1:2]
      new_final_region<-final_count[as.numeric(final_count[,1]) & as.numeric(final_count[,1])/max(as.numeric(final_count[,1])), 1:2]
      
      plot(rCount, main = Name, type = "h", xlab = "position", ylab = "mutation_count", col = "azure4", yaxt = "n")
      axis(2, at = seq(1, max(rCount[, 2]), ceiling((max(rCount[, 2]) - 0) / 3)))
      plot(final, main = name, type = "h", xlab = "position", ylab = "score")
    }   
    
    if (nrow(final_region) != 0) {
      plot(rCount, main = Name, type = "h", xlab = "position", ylab = "mutation_count", col = "azure4", yaxt = "n")
      axis(2, at = seq(1, max(rCount[, 2]), ceiling((max(rCount[, 2]) - 0) / 3)))
      plot(final, main = name, type = "h", xlab = "position", ylab = "score")
      #pvalue_averge_plot <- permutation(bb, Length, rCount, final, final_region, per_time, alph1)
      #pvalue_averge <- paste(pvalue_averge_plot, collapse = ",")
      
      tmp1 <- matrix(ncol = 2, nrow = (length(final_region) / 2))
      g = 1
      if (length(unique(bb[, 2])) == 1) {
        tmp1[1, 1] <- as.numeric(final_region[1, 2])
        tmp1[1, 2] <- as.numeric(final_region[1, 2])
        hotspot_position[1] <- length(which((rCount[tmp1[1, 1]:tmp1[1, 2], 2]) != 0))
        hotspot_region[1] <- ((tmp1[1, 2] - tmp1[1, 1]) + 1)
      } else {
        for (g in 1:nrow(final_region)) {
          tmp <- final_region[, 2]
          tmp1[g, 1] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[1]
          tmp1[g, 2] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[2]
          hotspot_position[g] <- length(which((rCount[tmp1[g, 1]:tmp1[g, 2], 2]) != 0))
          hotspot_region[g] <- ((tmp1[g, 2] - tmp1[g, 1]) + 1)
        }
      }
      
      tmp2 = NULL
      tmp2 <- matrix(ncol = 2, nrow = (length(final_region) / 2))
      tmp <- final_region[, 2]
      if (length(unique(bb[, 2])) == 1) {
        tmp2[1, 1] <- as.numeric(final_region[1, 2])
        tmp2[1, 2] <- as.numeric(final_region[1, 2])
        new_final_region <- final_region
      } else {
        g = 1
        for (g in 1:(length(final_region) / 2)) {
          tmp <- final_region[, 2]
          tmp2[g, 1] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[1]
          tmp2[g, 2] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[2]
          new_final_region <- final_region
        }
      }
      if (length(unique(bb[, 2])) != 1) {
        tmp3 <- matrix(ncol = 2, nrow = (length(final_region) / 2))
        new_final_region <- matrix(ncol = 2, nrow = (length(final_region) / 2))
        new_final_region[,1] <- final_region[,1]
        new_final_region[,2] <- final_region[,2]
        tmp <- final_region[, 2]
        g = 1
        for (g in 1:nrow(final_region)){
          tmp3[g, 1] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[1]
          tmp3[g, 2] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[2]
          if (tmp3[g,1] == tmp3[g,2]){
            #final_region[g,2]<- tmp3[g,1]
            new_final_region[g,2]<- tmp3[g,1]
          }
        }
      }
      
      #hotspot_information <- paste(hotspot_region, hotspot_position, pvalue_averge_plot, sep = ',')
      hotspot_information <- paste(hotspot_region, hotspot_position, sep = ',')
      detail <- as.matrix(hotspot_information)
      hotspot_information1 <- cbind(hotspot_region,hotspot_position)
      detail1 <- as.matrix(hotspot_information1)
      Total <- rep(sum(rCount[, 2]),nrow(final_region))
      maximun_count <- rep(max(rCount[, 2]),nrow(final_region))
      Amino_Acid_length <- rep(ALL$V5,nrow(final_region))
      cancer <- rep(as.character(bb$V5[1]),nrow(final_region))
      name1 <- rep(NAME,nrow(final_region))
      peak <- rep(peak_region,nrow(final_region))
      symbol1 <- rep(symbol,nrow(final_region))
      out_put <- cbind(cancer,symbol1,name1, Amino_Acid_length, Total, maximun_count,peak,new_final_region,detail1)
      write.table(out_put, file_output, col.name = F, row.name = F, append = TRUE, quote = F, sep = "\t")
      final_region[, 1] <- paste(paste(final_region[, 1], hotspot_information, sep = "("), "", sep = ")")
      plot_result(final, final_region, final_P, Length, rCount, bb, Name)
    }
  }
  #dev.off()
}


