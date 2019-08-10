setwd("C:/Users/Kikh/Desktop")
library(fabia)

ge <- read.table("GE.txt",header=TRUE,sep="\t")
#ge <- head(ge,200)
ge <- subset(ge[-which(colnames(ge)=='geneid')])
ge <- data.matrix(ge)
lbl <- read.table("gene.csv",header=TRUE,sep="\t")

res <- fabia(ge, 10, 0.1, cyc=400)
rb <- extractBic(res)

lbl <- read.table("gene.csv",header=TRUE,sep="\t")
new <- lbl$id
newa <- lbl$asthma
news <- lbl$currentsmoking
index <- 1
for (g in genes) { #loop over each bicluster
  lbla <- c()
  lbls <- c()
  if (!identical(g, character(0))) {
    print (index)
    sam <- samples[[index]] #the samples of the bicluster
    for (s in sam) {
      tmp <- which(new==s)
      lbla <- rbind(lbla,toString(newa[tmp]))
      lbls <- rbind(lbls,toString(news[tmp]))
    }
    write.table(lbla,file=paste("asthma_",index,".txt"))
    write.table(lbls,file=paste("smoke_",index,".txt"))
    write.table(rb$bic[index,]$bixv,file=paste("featvalues_",index,".txt"),row.names = FALSE,col.names = FALSE)
    write.table(rb$bic[index,]$biypv,file=paste("samvalues_",index,".txt"),row.names = FALSE,col.names = FALSE)
  }
  index <- index + 1
}  
details <- rb$bic #get the details of all biclusters
genes <- details[,'bixn'] #get the names of genes in each bicluster
samples <- details[,'biypn'] #get the names of samples in each bicluster

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
    write (sam, file=paste("heatsamples_",index,".txt")) #write the indices for the first run
    write (g, file=paste("heatfeatures_",index,".txt"))
    #print (g)
    
    #for every new bicluster
        for (o in 1:10) { #for every old bicluster
          count = 0
          features_old <- read.table(file=paste("1_features_",o,".txt"))
          samples_old <- read.table(file=paste("1_samples_",o,".txt"))
          for (f in g) {
            if (f %in% features_old$V1) {
                #print (f)
                for (s in sam) {
                  if (s %in% samples_old$V1) {
                    #print(s)
                    count = count + 1
                  }
                }
            }
          }
          print (paste("New bicluster ",index," with old bicluster ",o,"overlaps:"))
          print (count/(length(features_old$V1)*length(samples_old$V1)))
        }
    
    info_cont <- res@avini
    print (paste0("Information content: ",info_cont[index]))
    sum_info <- sum_info + info_cont[index]
    
    #create matrix representation for the bicluster
    for (i in g) {
      l <- c()
      for (j in sam) {
        l <- c(l,ge[i,j]) #the value of (gene,sample) of the initial matrix
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

#10 biclusters
b1 <- c(0,16,13,0,0,0,0,0,0,0)
b2 <- c(43,97,97,43,27,65,43,0,20,24)
b3 <- c(39,0,0,35,0,46,36,23,13,26)
b4 <- c(63,64,0,0,67,83,65,67,68,0)
b5 <- c(0,96,80,10,0,15,0,0,0,64)
b6 <- c(15,0,63,16,0,0,14,0,0,0)
b7 <- c(0,83,85,29,0,13,0,0,0,0)
b8 <- c(52,97,95,51,40,67,52,43,27,35)
b9 <- c(37,90,33,33,0,48,31,0,12,73)
b10 <- c(25,87,86,28,23,39,26,12,14,14)
test <- cbind(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10)
barplot(test,beside=T,xlab = 'Number of biclusters',ylab='Percentage of overlap')

setwd("C:/Users/Kikh/Desktop/new/featsam")
lbl <- read.table("gene.csv",header=TRUE,sep="\t")
for (n in 4:4) {
  print (paste("BICLUSTER ",n))
  #features <- read.table(file=paste("featvalues_",n,".txt"))
  samples <- read.table(file=paste("heatsamples_",n,".txt"))
  samples <- unlist(samples, recursive = TRUE, use.names = TRUE)
  #features <- features$x
  label_asthma <- lbl$asthma
  label_smoking <- lbl$currentsmoking
  label_med <- lbl$ics.use
  count_current <- 0
  count_clinical <- 0
  count_complete <- 0
  count_healthy <- 0
  count_smoking <- 0
  count_nonsmoking <- 0
  count_ics <- 0
  for (s in samples) {
    lbl1 <- which(label_asthma==s)
    lbl2 <- which(label_smoking==s)
    lbl3 <- which(label_med==s)
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

