selection <- function(df,row_file,col_file,number) {
  row <- read.table(row_file)
  col <- read.table(col_file)
  
  id <- 1:nrow(df)
  df <- cbind(id=id, df) #add id column to dataframe
  
  row_vector <- row[,1] #convert row indices dataframe col into vector
  df_rows <- df[df$id %in% row_vector,] #select rows with not empty values
  
  col_vector <- col[,1] #dataframe column into vector
  full_df <- df_rows[,names(df) %in% col_vector] #select columns with not empty values
  #print(full_df)
  #keep subset with numeric values
  sub_df <- subset(full_df[-which(colnames(full_df)=='asthma')])
  sub_df <- subset(sub_df[-which(colnames(sub_df)=='currentsmoking')])
  sub_df <- subset(sub_df[-which(colnames(sub_df)=='ics.use')])

  sub_df <- normal(sub_df) #normalize data 
  print (dim(sub_df))
  
  #create Ë and Z matrices for each bicluster
  #tmp <- t(sub_df)
  #setwd("C:/Users/Kikh/Desktop/bl_10")
  
  #for (i in 1:4) {
   # f_matr <- c()
    #s_matr <- c()
    #print (i)
    #feat <- read.table(file=paste("1_features_",i,".txt"))
    #sam <- read.table(file=paste("1_samples_",i,".txt"))
    #for (f in feat$V1) {
     # print(f)
    #  f_matr <- rbind(f_matr,tmp[f,])
   # }
  #  for (s in sam$V1) {
      #s <- toString(s)
     # s_matr <- rbind(s_matr,tmp[,s])
    #}
    #rownames(f_matr) <- feat$V1
   # rownames(s_matr) <- sam$V1
  #  write.table (f_matr, file=paste("feat_",i,".txt"),row.names = FALSE,col.names = FALSE)
   # write.table (s_matr, file=paste("sam_",i,".txt"),row.names = FALSE,col.names = FALSE)
  #}

  #tmp <- t(sub_df)
  #write.table(tmp[,"13", drop=FALSE],"col.txt",row.names = FALSE,col.names = FALSE)
  
  #write label columns in out files
  #new <- data.frame(rownames(full_df),full_df$asthma)
  #write.table(new, "asthma.txt", append = FALSE, sep = " ", dec = ".",
  #            row.names = FALSE, col.names = FALSE)
  
  #new <- data.frame(rownames(full_df),full_df$currentsmoking)
  #write.table(new, "smoking.txt", append = FALSE, sep = " ", dec = ".",
  #            row.names = FALSE, col.names = FALSE)
  
  #first run - counter matrices
  #initialize counter matrix to zero
  #zero_matr <- matrix(0L, nrow = ncol(sub_df), ncol = nrow(sub_df)) #transpose
  #row_names <- head(col, -3)
  #row_names <- unlist(row_names, use.names=FALSE)
  #rownames(zero_matr) <- row_names
  #col_names <- unlist(row, use.names=FALSE)
  #colnames(zero_matr) <- col_names
  
  #next runs
  #zero_matr<- read.table("counter.txt")
  #col_names <- colnames(zero_matr)
  #col_names <- sub('.', '', col_names)
  #colnames(zero_matr) <- col_names

  run_fabia(t(sub_df),number,zero_matr,full_df) #transpose because we want conditions (measurements) as rows
  #and samples(patients) as cols
}

normal<-function(m){ #z-score normalization
  print (class(m))
  mean = colMeans(m)
  sd = sapply(m, sd, na.rm = TRUE)
  for (j in 1:ncol(m)){
    for (i in 1:nrow(m)) {
      m[i,j] <- (m[i,j]-mean[j])/sd[j]
    }
  }
  return(m)
}

create_plot <- function(x,y1,y2,y3,y4,tl) { #plot the measurements of biclusters
  
  df <- data.frame("y1"=y1,"y2"=y2,"y3"=y3,"y4"=y4)
  df <- normal(df)
  
  y1 <- df$y1
  y2 <- df$y2
  y3 <- df$y3
  y4 <- df$y4
  
  plot(x, y1, main=tl, xlab="Number of biclusters",
       ylab="Quality measures", type="o", col="blue",
       pch="o", lty=1, ylim = c(-2,2))
  points(x, y2, col="red", pch="*")
  lines(x, y2, col="red",lty=2)
  points(x, y3, col="dark red",pch="+")
  lines(x, y3, col="dark red", lty=3)
  points(x, y4, col="dark green",pch="-")
  lines(x, y4, col="dark green", lty=4)
  op <- par(cex = 0.7)
  legend(1, 2, legend=c("Information content", "Variance", "MSR", "Virtual error"),
         col=c("blue", "red", "dark red", "dark green"), lty=1:4)
}


run_fabia <- function(data,number,zero_matr,full_df) {
  print("edw")
  print(dim(data))
  res <- fabia(data, number, 0.1)
  #summary(res)
  #extractPlot(res,ti="FABIA")
  #plot(res)
  print (data)
  rb <- extractBic(res)
  write.table(rb$bic[1,]$bixv,file="featvalues.txt",row.names = FALSE,col.names = FALSE)
  write.table(rb$bic[1,]$biypv,file="samvalues.txt",row.names = FALSE,col.names = FALSE)

  #print (rb$bic[1,]$biypn)
  #print (rb$bic[2,]$bixn)
  #plotBicluster(rb,5)
  details <- rb$bic #get the details of all biclusters
  genes <- details[,'bixn'] #get the names of genes in each bicluster
  samples <- details[,'biypn'] #get the names of samples in each bicluster
  
  index <- 1
  new <- data.frame(full_df$asthma)
  row.names(new) <- rownames(full_df)
  news <- data.frame(full_df$currentsmoking)
  row.names(news) <- rownames(full_df)
  #print(typeof(new))
  for (g in genes) { #loop over each bicluster
    lbla <- c()
    lbls <- c()
    if (!identical(g, character(0))) {
      sam <- samples[[index]] #the samples of the bicluster
      for (s in sam) {
        lbla <- rbind(lbla,toString(new[s,]))
        lbls <- rbind(lbls,toString(news[s,]))
      }
      write.table(lbla,file=paste("asthma_",index,".txt"))
      write.table(lbls,file=paste("smoke_",index,".txt"))
      write.table(rb$bic[index,]$bixv,file=paste("featvalues_",index,".txt"))
      write.table(rb$bic[index,]$biypv,file=paste("samvalues_",index,".txt"),row.names = FALSE,col.names = FALSE)
    }
    index <- index + 1
  }  
      
  measures(genes,res,samples,data,number,zero_matr)
}


measures <- function(genes,res,samples,data,number,zero_matr) {
  index <- 1
  matr <- c()
  sum_info <- 0
  sum_var <- 0
  sum_msr <- 0
  sum_ve <- 0
  for (g in genes) { #loop over each bicluster
    if (!identical(g, character(0))) { #for every non-empty biclusters (values over default threshold)
      print (paste0("BICLUSTER ",index))
      
      sam <- samples[[index]] #the samples of the bicluster
      write (sam, file=paste("1_samples_",index,".txt")) #write the indices for the first run
      write (g, file=paste("1_features_",index,".txt"))
      print (g)
    
      #for every new bicluster
#      for (o in 1:4) { #for every old bicluster
 #       count = 0
  #      features_old <- read.table(file=paste("1_features_",o,".txt"))
   #     samples_old <- read.table(file=paste("1_samples_",o,".txt"))
    #    for (f in g) {
     #     if (f %in% features_old$V1) {
            #print (f)
#            for (s in sam) {
 #            if (s %in% samples_old$V1) {
                #print(s)
  #              count = count + 1
   #           }
    #        }
     #     }
      #  }
  #      print (paste("New bicluster ",index," with old bicluster ",o,"overlaps:"))
   #     print (count/(length(features_old$V1)*length(samples_old$V1)))
    #  }
      
      info_cont <- res@avini
      print (paste0("Information content: ",info_cont[index]))
      sum_info <- sum_info + info_cont[index]
      
      #create matrix represantation for the bicluster
      for (i in g) {
        l <- c()
        for (j in sam) {
          l <- c(l,data[i,j]) #the value of (gene,sample) of the initial matrix
          #zero_matr[i,j] <- zero_matr[i,j]+1 
        }
        matr <- rbind(matr,l) #combine all the rows to a whole matrix
      }
      #print("matr")
      #print(matr)
      
      mean = mean(matr) #compute mean of whole matrix
      #compute variance of whole matrix/bicluster
      sum1 = 0
      for (i in 1:nrow(matr)) {
        sum2 = 0
        for (j in 1:ncol(matr)) {
          sum2 = sum2 + (matr[i,j]-mean)^2
        }
        sum1 = sum1 + sum2
      }
      variance = sum1
      print (paste0('Variance: ',variance))
      sum_var <- sum_var + variance
      
      #compute Mean Squared Residue (MSR)
      div = nrow(matr)*ncol(matr)
      row_mean = rowMeans(matr, na.rm = FALSE, dims = 1)
      col_mean = colMeans(matr, na.rm = FALSE, dims = 1)
      sum1 = 0
      for (i in 1:nrow(matr)) {
        sum2 = 0
        for (j in 1:ncol(matr)) {
          sum2 = sum2 + (matr[i,j]-row_mean[i]-col_mean[j]+mean)^2
        }
        sum1 = sum1 + sum2
      }
      msr = sum1/div
      print (paste0('MSR: ',msr))
      sum_msr <- sum_msr + msr
      
      #compute virtual error
      
      #virtual pattern
      p_list <- c()
      for (i in 1:nrow(matr)) {
        sum = 0
        for (j in 1:ncol(matr)) {
          sum = sum + matr[i,j]
        }
        p = sum/ncol(matr)
        p_list <- c(p_list,p)
      }
      #print('virtual pattern')
      plist <- as.list(p_list)
      #print (plist)
      
      #standardized bicluster matrix
      
      stand_matr <- c()
      mean_matr <- mean(matr)
      sd_matr <- sd(matr)
      #print("sd_matr")
      #print(sd_matr)
      if (!is.na(sd_matr) && !(sd_matr==0)) { #sd is nan when the matrix has only one element
        for (i in 1:nrow(matr)) {
          l <- c()
          for (j in 1:ncol(matr)) {
            l <- c(l,(matr[i,j]-mean_matr)/sd_matr) #standardization type
          }
          stand_matr <- rbind(stand_matr,l) #combine all the rows to a whole matrix
        }
      }
      else { #when the matrix has only one element then the standardized matrix
        stand_matr <- 0 #has only the value 0 (since standardized matrices have zero mean)
      }
      #print("standardization matrix")
      #print(stand_matr)
      
      #standardized p values
      mean_p <- mean(p_list)
      sd_p <- sd(p_list)
      stand_p <- c()
      #print("sd_p")
      #print (sd_p)
      if (!is.na(sd_p) && !(sd_p==0)) { #if not the matrix has only one element
        for (i in 1:length(plist)) {
          s <- as.numeric(unlist(plist[i]))
          stand_p <- c(stand_p,(s-mean_p)/sd_p)
        }
      }
      else {
        stand_p <- 0 #standardized matrices with one element have the value 0
      }
      #print("stand p")
      #print (stand_p)
      
      #virtual error (VE)
      if (stand_p[1]==0 && length(stand_p)==1) { #matrix with one element (ve = (0-0)/1 = 0)
        ve = 0
      }
      else {
        sum1 = 0
        for (i in 1:nrow(stand_matr)) {
          sum2 = 0
          for (j in 1:ncol(stand_matr)) {
            sum2 = sum2 + abs(stand_matr[i,j]-as.numeric(unlist(stand_p[i])))
          }
          sum1 = sum1 + sum2
        }
        ve = sum1 / (nrow(stand_matr)*ncol(stand_matr))
      }
      print (paste0("VE: ",ve))
      cat ("\n")
      sum_ve <- sum_ve + ve
      
      matr <- c() #initialize for next bicluster
    } #end of non-empty biclusters
    
    index <- index+1
    #break
  } #end of all biclusters
  #print (zero_matr)
  #write.table(zero_matr, file = "counter.txt", append = FALSE, quote = TRUE, sep = " ",
  #            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
  #            col.names = TRUE, qmethod = c("escape", "double"),
  #            fileEncoding = "")
  print("SUMS")
  print (paste0(sum_info," ",sum_var," ",sum_msr," ",sum_ve))
}

setwd("C:/Users/Kikh/Desktop")
library(fabia)
#library(ggplot2)

#read data tables and indices text files
df <- read.table("Database_biopten_v2_6.csv", header=TRUE, sep=";")
#ge <- read.table("GE.txt",header=TRUE,sep="\t")
#ge <- head(ge,200)
#lung
#selection(df,"row_lung.txt","col_lung.txt",8)
#blood
#selection(df,"row_blood.txt","col_blood.txt",10)
#sputum
#selection(df,"row_sputum.txt","col_sputum.txt",10)
#biopsy
selection(df,"row_biopsy.txt","col_biopsy.txt",6)
#lung & sputum
#selection(df,"row_lung_sputum.txt","col_lung_sputum.txt",10)
#blood & lung
#selection(df,"row_lung_blood.txt","col_lung_blood.txt",8)
#blood & sputum
#selection(df,"row_sputum_blood.txt","col_sputum_blood.txt",8)
#lung & blood & sputum
#selection(df,"row_lung_sputum_blood.txt","col_lung_sputum_blood.txt",10)

#create plots for quality measures of each category
#LUNG
x <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)
y1 <- c(739,1309,1476,2062,2253,2958,3216,3250,3557,3963,3835,4872,6800,7852,
        9928,8838,10515,10379,10461,8733,8406,9644,8547,7868,8349)
y2 <- c(285,384,612,661,850,899,967,1038,1191,1441,1357,1965,2913,3213,4314,3773,4940,
        4632,4777,3877,4058,4298,4207,3992,4542)
y3 <- c(0.37,0.58,0.93,1.03,1.34,1.6,1.86,1.78,1.98,2.04,2.1,2.25,2.54,2.64,
        2.8,2.96,3.43,3.4,3.58,2.97,3.24,3.36,3.69,3.23,3.66)
y4 <- c(1.33,2.89,3.06,3.71,4.24,4.7,5.77,5.17,5.76,5.69,5.99,6.4,7.64,8.8,9.71,
        9.2,10.63,10.95,10.91,9.56,9.37,10.4,10.66,9.77,10.5)
create_plot(x,y1,y2,y3,y4,"Lung function")

#SPUTUM
x <- c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)
y1 <- c(236,184,236,184,184,236,236,212,293,321,481,559,620,532,487,584)
y2 <- c(378,355,403,355,355,403,403,390,515,544,661,694,788,733,697,753)
y3 <- c(0.73,0.73,0.91,0.73,0.73,0.91,0.91,0.96,1.29,1.31,1.99,2.26,2.63,2.04,2.07,2.53)
y4 <- c(0.72,0.62,0.75,0.62,0.62,0.75,0.75,0.79,1.13,1.09,1.65,1.9,2.4,1.76,1.77,2.12)
create_plot(x,y1,y2,y3,y4,"Sputum")

#BLOOD
x <- c(2,4,6,8,10)
y1 <- c(0,192,1912,2781,2726)
y2 <- c(0,92.4,318,478,320)
y3 <- c(0,0.11,0.51,0.78,0.3)
y4 <- c(0,0.14,0.8,0.95,0.37)
create_plot(x,y1,y2,y3,y4,"Blood")

#LUNG - SPUTUM
x <- c(4,6,8,10,12,14,16)
y1 <- c(883,1388,1692,1951,1898,2164,2115)
y2 <- c(363,469,929,1141,1151,1218,1266)
y3 <- c(0.7,1.03,2.19,2.66,2.79,3.22,3.48)
y4 <- c(2.51,3.17,4.58,5.31,5.6,6.57,6.54)
create_plot(x,y1,y2,y3,y4,"Lung function - Sputum")

#BLOOD - SPUTUM
x <- c(4,6,8,10,12,14,16)
y1 <- c(1534,2251,3331,3864,4428,5406,6023)
y2 <- c(436,736,1016,1173,1344,1596,1705)
y3 <- c(0.6,1.24,1.72,1.94,2.27,2.58,3.07)
y4 <- c(3.27,4.28,6.31,6.75,6.57,7.8,8.46)
create_plot(x,y1,y2,y3,y4,"Blood - Sputum")

#BIOPSY
x <- c(2,4,6,8,10,12,14)
y1 <- c(1272,2842,3718,4321,5029,5205,3865)
y2 <- c(294,612,910,1060,1120,1343,3865)
y3 <- c(1.14,2.11,3.11,3.24,6.37,6.91,3.6)
y4 <- c(0.49,1.48,2.69,2.93,3.38,3.3,2.96)
create_plot(x,y1,y2,y3,y4,"Biopsy")

#barplots
#LUNG
#2 biclusters
b1 <- c(72,97,72,97,72,100,97,72,97,52)
b2 <- c(0,0,0,0,0,100,0,0,0,0)
test <- cbind(b1,b2)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#4 biclusters
b1 <- c(82,84,94,84,82,84,98,98,86,82)
b2 <- c(73,73,0,75,73,75,6,91,0,73)
b3 <- c(54,4,79,89,54,89,0,68,93,54)
b4 <- c(90,90,0,0,94,0,0,0,0,90)
test <- cbind(b1,b2,b3,b4)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#6 biclusters
b1 <- c(80,97,71,56,100,53,52,52,97,53)
b2 <- c(34,56,52,0,100,0,49,85,0,87)
b3 <- c(74,0,74,94,97,100,74,0,0,94)
b4 <- c(74,71,72,98,74,0,0,71,93,98)
b5 <- c(0,0,0,0,0,0,0,0,0,88)
b6 <- c(0,0,0,0,0,0,0,0,100,0)
test <- cbind(b1,b2,b3,b4,b5,b6)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#8 biclusters
b1 <- c(98,98,84,84,86,96,79,86,98,100)
b2 <- c(95,0,100,41,55,65,67,41,55,0)
b3 <- c(21,21,86,63,75,61,61,63,78,5)
b4 <- c(0,90,96,90,97,72,98,93,98,90)
b5 <- c(0,92,96,92,100,72,97,94,98,92)
b6 <- c(98,0,89,0,100,0,0,0,0,0)
b7 <- c(97,0,91,0,98,0,0,0,0,0)
b8 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')


#barplots 
#SPUTUM
#2 biclusters
b1 <- c(100,100,100,0,0,100,0,100,0,100)
b2 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#4 biclusters
b1 <- c(0,0,100,100,100,100,100,0,100,100)
b2 <- c(0,0,0,0,0,0,0,0,0,0)
b3 <- c(0,0,0,0,0,0,0,0,0,0)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#6 biclusters
b1 <- c(100,62,100,100,100,100,100,100,62,100)
b2 <- c(0,0,0,0,0,0,0,0,0,0)
b3 <- c(0,0,0,0,0,0,0,0,0,0)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
b5 <- c(0,0,0,0,0,0,0,0,0,0)
b6 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#8 biclusters
b1 <- c(100,100,100,100,100,100,100,100,100,100)
b2 <- c(0,0,0,0,0,0,0,0,0,0)
b3 <- c(0,0,0,0,0,0,0,0,0,0)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
b5 <- c(0,0,0,0,0,0,0,0,0,0)
b6 <- c(0,0,0,0,0,0,0,0,0,0)
b7 <- c(0,0,0,0,0,0,0,0,0,0)
b8 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')


#barplots 
#BLOOD
#4 biclusters
b1 <- c(0,0,100,100,100,0,0,0,0,0)
b2 <- c(0,0,0,0,0,0,0,0,0,0)
b3 <- c(0,0,0,0,0,0,0,0,0,0)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#6 biclusters
b1 <- c(97,0,100,97,98,100,97,97,97,97)
b2 <- c(0,0,100,0,43,43,0,0,0,0)
b3 <- c(0,0,0,0,0,0,0,0,0,0)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
b5 <- c(0,0,0,0,0,0,0,0,0,0)
b6 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#8 biclusters
b1 <- c(88,100,100,98,98,84,92,100,100,97)
b2 <- c(90,100,100,98,98,100,100,100,100,0)
b3 <- c(90,100,100,0,0,96,98,98,96,98)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
b5 <- c(0,0,0,0,0,0,0,0,0,0)
b6 <- c(0,0,0,0,0,0,0,0,0,0)
b7 <- c(0,0,0,0,0,0,0,0,0,0)
b8 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#10 biclusters
b1 <- c(92,92,98,100,92,96,98,92,70,100)
b2 <- c(94,98,100,100,98,100,100,98,98,77)
b3 <- c(98,98,100,100,98,100,100,100,98,75)
b4 <- c(0,0,0,0,0,0,0,0,0,0)
b5 <- c(0,0,0,0,0,0,0,0,0,0)
b6 <- c(0,0,0,0,0,0,0,0,0,0)
b7 <- c(0,0,0,0,0,0,0,0,0,0)
b8 <- c(0,0,0,0,0,0,0,0,0,0)
b9 <- c(0,0,0,0,0,0,0,0,0,0)
b10 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#BIOPSY
#6 biclusters
b1 <- c(100,100,100,100,100,100,100,100,100,100)
b2 <- c(100,100,100,100,93,100,100,100,100,100)
b3 <- c(100,100,100,66,100,100,69,48,66,100)
b4 <- c(100,100,96,100,100,96,100,45,0,45)
b5 <- c(92,91,100,100,0,92,43,100,100,100)
b6 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6)
barplot(test,main="Biopsy",beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#LUNG FUNCTION - SPUTUM
#10 biclusters
b1 <- c(64,56,56,56,58,73,51,78,52,54)
b2 <- c(40,73,25,29,39,86,73,29,97,34)
b3 <- c(100,100,70,100,100,54,100,100,100,100)
b4 <- c(60,72,60,75,75,75,72,87,90,11)
b5 <- c(100,96,96,96,83,96,73,96,100,26)
b6 <- c(98,72,65,72,73,54,72,100,100,16)
b7 <- c(94,0,48,0,0,44,72,97,97,0)
b8 <- c(0,0,0,0,0,0,0,0,0,0)
b9 <- c(0,0,0,0,0,0,0,0,0,0)
b10 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#LUNG FUNCTION - BLOOD
#8 biclusters
b1 <- c(53,53,52,53,78,93,69,54,50,66)
b2 <- c(94,92,100,86,92,92,98,92,88,90)
b3 <- c(100,0,100,0,0,0,100,0,100,100)
b4 <- c(0,84,0,56,56,90,0,90,68,52)
b5 <- c(66,66,66,64,94,37,92,71,66,89)
b6 <- c(7,64,5,62,41,86,5,84,79,62)
b7 <- c(68,73,10,64,73,96,7,0,7,68)
b8 <- c(0,0,0,0,0,0,0,0,0,0)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

#find the number of each class in every bicluster
for (n in 4:5) {
  print (paste("BICLUSTER ",n))
  #features <- read.table(file=paste("featvalues_",n,".txt"))
  samples <- read.table(file=paste("1_samples_",n,".txt"))
  #features <- features$x
  label_asthma <- df$asthma
  label_smoking <- df$currentsmoking
  label_med <- df$ics.use
  count_current <- 0
  count_clinical <- 0
  count_complete <- 0
  count_healthy <- 0
  count_smoking <- 0
  count_nonsmoking <- 0
  count_ics <- 0
  for (s in samples$V1) {
    lbl1 <- label_asthma[s]
    lbl2 <- label_smoking[s]
    lbl3 <- label_med[s]
    if (lbl1 == 'current asthma') {
      count_current = count_current + 1
    }
    else if (lbl1 == 'clinical remission') {
      count_clinical = count_clinical + 1
    }
    else if (lbl1 == 'complete remission') {
      count_complete = count_complete + 1
    }
    else if (lbl1 == 'healthy') {
      count_healthy = count_healthy + 1
    }
    
    if (lbl2 == "smoker") {
      count_smoking = count_smoking + 1
    }
    else if (lbl2 == "nonsmoker") {
      count_nonsmoking = count_nonsmoking + 1
    }

    if (!is.na(lbl3)) {
      if (lbl3 == 'ics use') {
        count_ics = count_ics + 1
      }
    }
  }
  print (paste('*****ASTHMA*****',count_current+count_clinical+count_complete))
  print (paste("Current asthma: ",count_current))
  print (paste("Clinical remission: ",count_clinical))
  print (paste("Complete remission: ",count_complete))
  print (paste("*****HEALTHY*****",count_healthy))
  cat("\n")
  print (paste("Smoking: ",count_smoking))
  print (paste("Not smoking: ",count_nonsmoking))
  cat ("\n")
  print (paste("Ics use: ",count_ics))
  cat ("\n")
}



#res <- fabia(lung, 6, 0.1, 400)
#summary(res)
#extractPlot(res,ti="FABIA")
#plot(res)
#rb <- extractBic(res)
#plotBicluster(rb,5)


