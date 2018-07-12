######MAIN FUNCTION#######
##### M-hotspot #####
Main_fun <- function(x) {
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



###### 為計算每個位點上之突變數目 ######
### Length 為蛋白質的長度、bb為突變資料 ###

each_position <- function(Length, bb) {
    tt <- Length
    TL <- matrix(1:Length, ncol = 1)
    TL <- as.matrix(TL, ncol = 2)
    zxc <- function(y) {
        num <- length(which(bb[, 2] %in% y))
        return(num)
    }
    Count <- apply(TL, 1, zxc)
    Count <- as.matrix(Count, ncol = 1)
    rCount <- as.data.frame(cbind(TL, Count))
    return(rCount)
}

##############################

###### Score 轉換(density information)######
### Length為蛋白質的總長、rCount為每個位點的突變個數 ###

#Density_information <- function(Length, alph1, rCount) {
#    final <- NULL
#    i = 1
#    for (i in 1:Length) {
#        final_score = 0
#        k = 1
#        for (k in k:Length) {
#            score <- exp((-abs(i - k)) / (alph1 ^ 2)) * rCount[k, 2]
#            final_score <- score + final_score
#        }
#        final[i] <- final_score
#    }
#    return(final)
#}

new_Density_information <- function(Length, alph1, new_rCount) {
    final <- NULL
    i = 1
    for (i in 1:Length) {
        final_score = 0
        k = 1
        for (k in k:nrow(new_rCount)) {
            score <- exp((-abs(i - new_rCount[k, 1])) / (alph1 ^ 2)) * new_rCount[k, 2]
            final_score <- score + final_score
        }
        final[i] <- final_score
    }
    return(final)
}
#new_Density_information <- function(Length,alph1,rCount) {
#    test <- as.matrix(rCount$position)
#    final <- NULL
#    #new_final <- matrix(ncol = 1, nrow = Length)
#    density <- function(k){
#        count_num <- rCount[k,2]
#        score <- exp((-abs(i - k)) / (alph1 ^ 2)) * count_num
#        #final_score <- sum(score)
#    }
#    for (i in 1:Length) {
#        final_score = 0
#        #k = 1
#        score <- apply(test, 1, density)
#        final_score <- sum(score)
#        final[i] <- final_score
#    }
#    return(final)
#}
##############################

##### 延伸距離 #####
### Length為蛋白質的總長、mutation為有發生突變的蛋白質位點 ###

Extend_region <- function(peak_region_min, peak_region_max, remain_gap, mutation) {
    if (length(mutation) == 1) {
        peak_region <- 10
    } else {
        count <- mutation[-length(mutation)]
        count1 <- mutation[-c(1)]
        gap1 <- count1 - count
        gap <- sort(gap1)
        new_gap <- gap[1:(ceiling(length(gap) * remain_gap))]
        if (length(new_gap) == 1) {
            peak_region <- ceiling(mean(new_gap))
        } else {
            peak_region <- ceiling(mean(new_gap) + sd(new_gap))
        }
    }
    
    ### 判斷是否有落入我們自定的範圍,超出的話就給予最大或最小值 ###
    
    if (peak_region < peak_region_min) {
        peak_region = peak_region_min
    }
    if (peak_region > peak_region_max) {
        peak_region = peak_region_max
    }
    return(peak_region)
}

##############################

##### 山岳刪減分群法 #####
### Length為蛋白質的總長、mutation為有發生突變的蛋白質位點、final為經過轉換過後各位點的分數 ###

mountain_clusting <- function(mutation, final, num_center, Length, beta) {
    beta_value = beta / (alph1 ^ 2)
    final_P = NULL
    final_P <- as.list(final_P)
    P1 <- which(final == max(final))
    w = 1
    PP = NULL
    P_score = NULL
    for (w in 1:num_center) {
        if (w == 1) {
            ### Score 最高的點當做第一群中心 ###		
            PP[w] <- P1
        } else {
            j = 1
            for (j in 1:Length) {
                if (w == 2) {
                    ### 刪減第一群中心的影響找出第二群中心 ###	
                    P_score[j] <- (final[j] - ((max(final)) * (exp(-(beta_value * (P1 - j) ^ 2)))))
                } else {
                    ### 以此類推, 找出第剩餘的群中心 ###	
                    P_score[j] <- (P_score[j] - ((max(P_score)) * (exp(-(beta_value * (P - j) ^ 2)))))
                }
            }
            P <- which(P_score == max(P_score))
            PP[w] <- P
        }
        final_P[[1]] <- PP
    }
    
    ### 找出的群中心一定要是在有發生突變的位點上 ###	
    final_P[[1]] <- unique(final_P[[1]])
    final_P[[1]] <- final_P[[1]][which(final_P[[1]] %in% mutation)]
    return(final_P[[1]])
}

##############################

##### 從中心點延伸並找出hotspt #####
### Length為蛋白質的總長、mutation為有發生突變的蛋白質位點、final為經過轉換過後各位點的分數 ###
### grade為與其他的突變位置距離群中心多遠 ###

extend <- function(final_P, peak_region, mutation, rCount) {
    position_star <- NULL
    position_end <- NULL
    e = 1
    for (e in 1:(length(final_P))) {
        grade <- mutation - final_P[e]
        center <- as.numeric(which(grade == 0))
        
        ### 如果群中心點在發生突變為點的最兩側(第一位點或最後一個位點), 會另外處理(只向其中一邊做延伸) ###
        
        if (center == 1 | center == (length(grade))) {
            if (center == 1) {
                position_star[e] <- mutation[center]
                o = 1
                for (o in 1:(length(mutation) - center)) {
                    ### 利用 mutation 資訊(有發生突變的蛋白質位點) 做延伸 ###
                    if (o == 1) {
                        if (abs(grade[center + o]) <= peak_region) {
                            position_end[e] <- mutation[center + o]
                        } else {
                            position_end[e] <- mutation[center]
                            break
                        }
                    } else {
                        if (abs(abs(grade[center + o]) - abs(grade[center + o - 1])) <= peak_region) {
                            position_end[e] <- mutation[center + o]
                        } else {
                            position_end[e] <- mutation[center + o - 1]
                            break
                        }
                    }
                }
            } else {
                position_end[e] <- mutation[center]
                u = 1
                for (u in 1:(center - 1)) {
                    if (u == 1) {
                        if ((abs(grade[center - u])) <= peak_region) {
                            position_star[e] <- mutation[center - u]
                        } else {
                            position_star[e] <- mutation[center]
                            break
                        }
                    } else {
                        if (abs(abs(grade[center - u]) - abs(grade[center - u + 1])) <= peak_region) {
                            position_star[e] <- mutation[center - u]
                        } else {
                            position_star[e] <- mutation[center - u + 1]
                            break
                        }
                    }
                }
            }
        } else {
            u = 1
            for (u in 1:(center - 1)) {
                if (u == 1) {
                    if ((abs(grade[center - u])) <= peak_region) {
                        position_star[e] <- mutation[center - u]
                    } else {
                        position_star[e] <- mutation[center]
                        break
                    }
                } else {
                    if (abs(abs(grade[center - u]) - abs(grade[center - u + 1])) <= peak_region) {
                        position_star[e] <- mutation[center - u]
                    } else {
                        position_star[e] <- mutation[center - u + 1]
                        break
                    }
                }
            }
            o = 1
            for (o in 1:(length(mutation) - center)) {
                if (o == 1) {
                    if (abs(grade[center + o]) <= peak_region) {
                        position_end[e] <- mutation[center + o]
                    } else {
                        position_end[e] <- mutation[center]
                        break
                    }
                } else {
                    if (abs(abs(grade[center + o]) - abs(grade[center + o - 1])) <= peak_region) {
                        position_end[e] <- mutation[center + o]
                    } else {
                        position_end[e] <- mutation[center + o - 1]
                        break
                    }
                }
            }
        }
    }
    
    final_center <- unlist(final_P)
    final_region <- paste(position_star, final_center, position_end, sep = "-")
    test_region <- paste(position_star, position_end, sep = "-")
    uniq_region <- unique(test_region)
    uniq_star <- unique(position_star)
    uniq_end <- unique(position_end)
    s = 1
    m_count = NULL
    for (s in 1:length(uniq_star)) {
        m_count[s] <- sum(rCount[uniq_star[s]:uniq_end[s], 2])
    }
    final_count <- cbind(as.matrix(m_count), as.matrix(uniq_region))
    final_count <- cbind(final_count, as.numeric(final_count[, 1]) / sum(as.numeric(final_count[, 1])), as.numeric(final_count[, 1]) / (uniq_end - uniq_star + 1))
    return(final_count)
}

##############################

##### 利用permutation方法算出p-value值 #####
### Length為蛋白質的總長、rCount和Count為每個位點的突變個數、final為經過轉換過後各位點的分數 ###

#permutation <- function(bb, Length, rCount, final, final_region, per_time, alph1) {
#    per_final = as.list(1:per_time)
#    per_final_region <- NULL
#    l = 1
#    
#    ### 固定position 的個數並隨機放入相同數量的突變資訊進去 ###
#    
#    for (l in 1:per_time) {
#        mm <- cbind(rCount[which(rCount[, 2] != 0),])
#        permutation <- cbind(sort(sample(1:Length, nrow(mm))), mm$mutation)
#        mm_count <- sum(rCount[, 2])
#        if (nrow(permutation) != 1) {
#            per_random_mutation <- sort(sample(permutation[, 1], mm_count, replace = TRUE))
#            i = 1
#            for (i in 1:nrow(permutation)) {
#                permutation[i, 2] <- length(which(per_random_mutation == permutation[i, 1]))
#            }
#        }
#        
#        per_Count <- rCount
#        per_Count[, 2] <- 0
#        per_Count[which(per_Count[, 1] %in% permutation[, 1]), 2] <- permutation[, 2]
#        
#        ### 對Random 資料做Score 轉換 ###
#        
#        per_final_P = NULL
#        if (sum((rCount[, 2] == 0) == FALSE) == 1) {
#            i = 1
#           for (i in 1:nrow(per_Count)) {
#               per_final_score = 0
#                k = 1
#                for (k in k:nrow(per_Count)) {
#                    per_score <- exp((-abs(i - k)) / (alph1 ^ 2)) * per_Count[k, 2]
#                    per_final_score <- per_score + per_final_score
#                }
#                per_final[[l]][i] <- per_final_score
#            }
#            per_final_P[[1]] <- which(per_final[[l]] == max(per_final[[l]]))
#        } else {
#            i = 1
#            for (i in 1:nrow(per_Count)) {
#                per_final_score = 0
#                k = 1
#                for (k in k:nrow(per_Count)) {
#                    per_score <- exp((-abs(i - k)) / (alph1 ^ 2)) * per_Count[k, 2]
#                    per_final_score <- per_score + per_final_score
#                }
#                per_final[[l]][i] <- per_final_score
#            }
#            per_P1 <- which(per_final[[l]] == max(per_final[[l]]))
#       }
#    }
#    tmp4 = matrix(ncol = per_time, nrow = nrow(final_region))
    
    ### 比較Random 資料與原始資料之 hotspot 內的平均分數高低 ###
    
#    options(digits = 3)
#    w = 1
#    for (w in 1:per_time) {
#        tmp1 = NULL
#        tmp1 <- matrix(ncol = 2, nrow = (length(final_region) / 2))
#        if (length(unique(bb[, 2])) == 1) {
#            tmp1[1, 1] <- as.numeric(final_region[1, 2])
#            tmp1[1, 2] <- as.numeric(final_region[1, 2])
#        } else {
#            g = 1
#            for (g in 1:(length(final_region) / 2)) {
#                tmp <- final_region[, 2]
#                tmp1[g, 1] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[1]
#                tmp1[g, 2] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[2]
#            }
#        }
#        i = 1
#        for (i in 1:(length(final_region) / 2)) {
#            REGION <- final[tmp1[i, 1]:tmp1[i, 2]]
#            per_REGION <- per_final[[w]][tmp1[i, 1]:tmp1[i, 2]]
#            if (mean(per_REGION) > mean(REGION)) {
#                tmp4[i, w] <- 1
#            }
#        }
#        tmp4[is.na(tmp4)] <- 0
#        pvalue_averge <- (rowSums(tmp4)) / per_time
#        final_tmp4 <- cbind(final_region, pvalue_averge, rowSums(tmp4))
#    }
#    
#    
#    return(pvalue_averge)
#}

#########################
##### 畫結果圖 #####

plot_result <- function(final, final_region, final_P, Length, rCount, bb, Name) {
    Max <- which(final == max(final))
    d2 <- cbind(as.data.frame(final_P), which(final_P == final_P))
    colnames(d2) <- c("position", "mutation")
    mid <- NULL
    #ellipse <- NULL
    if (length(unique(bb[,2])) == 1) {
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(2, nrow(final_region))))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        #ellipse <- NULL
        mid <- (Max - 0.5)
        #test_ellipse <- as.data.frame(getellipse(rx = 0.5, ry = (max(rCount[, 2]) / 6), mid = c((Max - 0.5), (max(rCount[, 2]) / 4))))
        #ellipse <- rbind(ellipse, test_ellipse)
        #colnames(ellipse) <- c("x", "y")
        mid <- as.data.frame(mid, ncol = 1)
        tmp <- final_region[, 2]
        #tmp1 <- as.numeric(tmp)
        #tmp1[2] <- tmp1[1]-1
        #tmp1[1] <- tmp1[1]+5
        bb$hotspot=''
        bb[bb$position %in% tmp,6] <- tmp
        
        mid <- as.data.frame(mid, ncol = 1)
        colnames(final_region) <- c("V1","V2")
        GGG <- ggplot() +
            geom_histogram(data = bb, aes(x = position, fill =hotspot), bins = Length, breaks = 1:length(rCount[, 2]))+
            #geom_histogram(data = bb, aes(x = position), bins = Length, breaks = 1:length(rCount[, 2]), fill = "black") +
            scale_y_continuous(limits = c(0, max(rCount[, 2])), breaks = seq(min(rCount[which(rCount[, 2] != 0), 2]), max(rCount[which(rCount[, 2] != 0), 2]), ceiling((max(rCount[which(rCount[, 2] != 0), 2]) - min(rCount[which(rCount[, 2] != 0), 2])) / 3))) +
            geom_text_repel(data = mid, aes(x = mid, y = (max(rCount$mutation) / 4) + (max(rCount[, 2]) / 4), label = (as.data.frame(final_region))$V1), fontface = 'bold', color = 'black', size = 5) +
            #annotate("rect", xmin = re1[1,2], xmax = re1[1,1], ymin =  re1[1,3], ymax =  re1[1,4] ,fill="transparent",color="red",size = 0.5)+
            #geom_point(data = ellipse, aes(x = x, y = y), colour = "red", size = 0.5) +
            theme(axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 9, color = "black", face = "bold"),panel.background = element_rect(fill = "white", colour = "grey50")) +
            ggtitle(Name)
        print(GGG, vp = vplayout(1, 1))
        
        GG <- ggplot() +
            geom_histogram(data = bb, aes(x = position,fill =hotspot), bins = Length, breaks = (Max - 10):(Max + 10)) +
            #geom_histogram(data = bb, aes(x = position), bins = Length, breaks = (Max - 10):(Max + 10), fill = "black") +
            scale_y_continuous(limits = c(0, max(rCount[, 2])), breaks = seq(min(rCount[which(rCount[, 2] != 0), 2]), max(rCount[which(rCount[, 2] != 0), 2]), ceiling((max(rCount[which(rCount[, 2] != 0), 2]) - min(rCount[which(rCount[, 2] != 0), 2])) / 3))) +
            geom_text(data = as.data.frame(mid), aes(x = as.data.frame(mid), y = (max(rCount$mutation) / 4) + (max(rCount[, 2]) / 4), label = (as.data.frame(final_region))$V1[1]), fontface = 'bold', color = 'black', size = 5) +
            #annotate("rect", xmin = re1[1,2], xmax = re1[1,1], ymin =  re1[1,3], ymax =  re1[1,4] ,fill="transparent",color="red",size = 0.5)+
            #geom_point(data = ellipse, aes(x = x, y = y), colour = "red", size = 0.5) +
            theme(axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 9, color = "black", face = "bold"),panel.background = element_rect(fill = "white", colour = "grey50")) +
            ggtitle(final_region[1, 2])+guides(fill=FALSE)
        print(GG, vp = vplayout(2, 1))
        
    } else {
        
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(2, nrow(final_region))))
        vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
        
        tmp1 <- matrix(nrow = nrow(final_region),ncol = 2)
        g = 1
        for (g in 1:nrow(final_region)) {
            tmp <- final_region[, 2]
            tmp1[g,] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))
            tmp1[g,1] <- tmp1[g,1]
            mid[g] <- tmp1[g,1] + ((tmp1[g,2] - tmp1[g,1]) / 2)
            
            #re <- as.data.frame(tmp1)
            #re_y <-  as.data.frame(t(rep(c(0,10),1)))
            #re1 <- cbind(re,re_y)
            #colnames(re1) <- c("x", "y")
            #test_ellipse <- as.data.frame(getellipse(rx = ((tmp1[2] - tmp1[1]) / 2), ry = (max(rCount[, 2]) / 6), mid = c(tmp1[1] + ((tmp1[2] - tmp1[1]) / 2), (max(rCount[, 2]) / 4))))
            #ellipse <- rbind(ellipse, test_ellipse)
        }
        #colnames(ellipse) <- c("x", "y")
        bb$hotspot=''
        i=1
        for (i in 1:nrow(final_region)){
            if (tmp1[i,1] == tmp1[i,2]){
                #final_region[g,2]<- tmp3[g,1]
                bb[bb$position %in% c(tmp1[i,1]:tmp1[i,2]),6] <- tmp1[i,1]
            } else { 
                bb[bb$position %in% c(tmp1[i,1]:tmp1[i,2]),6] <- tmp[i]
            }
        }
        bb1 <- bb
        bb1[which(bb1[,6]==""),6] <-"nohotspot"
        mid <- as.data.frame(mid, ncol = 1)
        #re <- as.data.frame(tmp1)
        #re_y <- as.data.frame(c(0,10), ncol = 1)
        #re1 <- cbind(re,re_y)
        #colnames(re1) <- c("x", "y")
        GGG <- #ggplot() +
            #geom_histogram(data = bb1, aes(x = position, fill =hotspot), bins = Length, breaks = 1:length(rCount[, 2]))+
            ggplot(bb1,aes(x=position)) +
            geom_histogram(data=subset(bb1,hotspot !="nohotspot"),aes(x = position,fill = hotspot), bins = Length, breaks = 1:length(rCount[, 2])) +
            geom_histogram(data=subset(bb1,hotspot == "nohotspot"),aes(x = position),fill = "gray", bins = Length, breaks = 1:length(rCount[, 2]))+
            #geom_histogram(data = bb, aes(x = position), bins = Length, breaks = 1:length(rCount[, 2]), fill = "black") +
            scale_y_continuous(limits = c(0, max(rCount[, 2])), breaks = seq(min(rCount[which(rCount[, 2] != 0), 2]), max(rCount[which(rCount[, 2] != 0), 2]), ceiling((max(rCount[which(rCount[, 2] != 0), 2]) - min(rCount[which(rCount[, 2] != 0), 2])) / 3))) +
            #geom_ribbon(aes(x = mid,ymin=tmp1[1], ymax=tmp1[1]))
            #geom_text_repel(data = mid, aes(x = mid, ymin=tmp1[1], ymax=tmp1[1], label = (as.data.frame(final_region))$V1), fontface = 'bold', color = 'blue', size = 5) +
            geom_text_repel(data = mid, aes(x = mid, y = (max(rCount$mutation) / 4) + (max(rCount[, 2]) / 4), label = (as.data.frame(final_region))$V1), fontface = 'bold', color = 'black', size = 5) +
            #geom_ribbon(aes(xmin = tmp1[1],xmax = tmp1[2],ymin = mid[1]-5,ymax = mid[1]+5), fill = "red")
            #geom_crossbar(data = re1,aes(xmin = re1[1,1],xmax = re1[1,2] ,ymin = re1[1,2],ymax = re1[2,2]), width = 0.2)
            #annotate("rect", xmin = re1[,1], xmax = re1[,2], ymin =  re1[,3], ymax =  re1[,4] ,fill="transparent",color="red",size = 0.5)+
            #geom_ribbon(data = re1,aes(x = x ,ymin = re1[1,2],ymax = re1[2,2]),fill = "red",alpha = .5)
            #geom_rect(xmin = tmp1[1],xmax = tmp1[2],ymin = mid[1]-5,ymax = mid[1]+5)
            #geom_point(data = ellipse, aes(x = x, y = y), colour = "red", size = 0.5) +
            theme(axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 9, color = "black", face = "bold"),panel.background = element_rect(fill = "white", colour = "grey50")) +
            ggtitle(Name)
        
        
        print(GGG, vp = vplayout(1, (1:nrow(final_region))))
        
        tmp2 = NULL
        tmp2 <- matrix(ncol = 2, nrow = (length(final_region) / 2))
        r = 1
        for (r in 1:nrow(final_region)) {
            mid <- NULL
            #ellipse <- NULL
            tmp <- final_region[, 2]
            tmp1 <- as.numeric(unlist(strsplit(tmp[r], split = "-", fix = T)))
            mid[r] <- tmp1[1] + ((tmp1[2] - tmp1[1]) / 2)
            #test_ellipse <- as.data.frame(getellipse(rx = (((tmp1[2] - tmp1[1]) / 2) + 0.5), ry = (max(rCount[, 2]) / 6), mid = c(tmp1[1] + (((tmp1[2] - tmp1[1]) / 2) - 0.5), (max(rCount[, 2]) / 4))))
            #ellipse <- rbind(ellipse, test_ellipse)
            tmp <- final_region[, 2]
            if (length(unique(bb[, 2])) == 1) {
                tmp2[1, 1] <- as.numeric(final_region[1, 2])
                tmp2[1, 2] <- as.numeric(final_region[1, 2])
            } else {
                g = 1
                for (g in 1:(length(final_region) / 2)) {
                    tmp <- final_region[, 2]
                    tmp2[g, 1] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[1]
                    tmp2[g, 2] <- as.numeric(unlist(strsplit(tmp[g], split = "-", fix = T)))[2]
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
            GG <- #ggplot() +
                #geom_histogram(data = bb1, aes(x = position,fill = hotspot), bins = Length, breaks =(tmp2[r, 1] - 10):(tmp2[r, 2] + 10))+
                ggplot(bb1,aes(x=position)) +
                geom_histogram(data=subset(bb1,hotspot !="nohotspot"),aes(x = position,fill = hotspot), bins = Length, breaks =(tmp2[r, 1] - 10):(tmp2[r, 2] + 10)) +
                geom_histogram(data=subset(bb1,hotspot == "nohotspot"),aes(x = position),fill = "gray", bins = Length, breaks =(tmp2[r, 1] - 10):(tmp2[r, 2] + 10))+
                #geom_histogram(data = bb, aes(x = position), bins = Length, breaks = (tmp2[r, 1] - 10):(tmp2[r, 2] + 10), fill = "black") +
                scale_y_continuous(limits = c(0, max(rCount[, 2])), breaks = seq(min(rCount[which(rCount[, 2] != 0), 2]), max(rCount[which(rCount[, 2] != 0), 2]), ceiling((max(rCount[which(rCount[, 2] != 0), 2]) - min(rCount[which(rCount[, 2] != 0), 2])) / 3))) +
                geom_text(data = as.data.frame(mid), aes(x = as.data.frame(mid), y = (max(rCount$mutation) / 4) + (max(rCount[, 2]) / 4), label = (as.data.frame(final_region))$V1[r]), fontface = 'bold', color = 'black', size = 5) +
                #geom_point(data = ellipse, aes(x = V1, y = V2), colour = "red", size = 0.5) +
                #annotate("rect", xmin = re1[r,1], xmax = re1[r,2], ymin =  re1[r,3], ymax =  re1[r,4] ,fill="transparent",color="red",size = 0.5)+
                theme(axis.text = element_text(size = 9, color = "black"), axis.title = element_text(size = 9, color = "black", face = "bold"),panel.background = element_rect(fill = "white", colour = "grey50")) +
                ggtitle(new_final_region[r, 2])+guides(fill=FALSE)
            print(GG, vp = vplayout(2, r))
        }
    }
}

##############################
