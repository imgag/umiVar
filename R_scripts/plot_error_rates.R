suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(seqinr)) 


parser <- ArgumentParser()

# setting parameters
parser$add_argument("-bc1", "--bc1_file", type="character", help="Barcode corrected 1", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-bc2", "--bc2_file", type="character", help="Barcode corrected 2", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-bc3", "--bc3_file", type="character", help="Barcode corrected 3", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-bc4", "--bc4_file", type="character", help="Barcode corrected 4 or more", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-out", "--out_folder", type="character", help="out_folder", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()

bc1 <- args$bc1_file

bc2 <- args$bc2_file

bc3 <- args$bc3_file

bc4 <- args$bc4_file


OUTF <- args$out_folder



#### SOME FUNCTIONS
ER_NT <- function(x){
  NT <- c("A","C","T","G")
  TAB <- NULL
  DP <- sum(as.numeric(x$DP))
  for (i in NT){
    X <- x[x$REF == i,]
    Z <- NT[!(NT == i)]
    for (y in Z){
      COUNT <- tryCatch(as.numeric(sum(X[y])), error = function(e) NA)
      TAB2 <- data.frame(i,y,sum(X$DP),COUNT,DP,COUNT/DP,COUNT/sum(X$DP))
      colnames(TAB2) <- c("REF","ALT","DP_REF","ALT_COUNT","DP_TOTAL", "ER_TOTAL","ER_REF")
      TAB <- rbind(TAB,TAB2)
    }
  }
  return(TAB)
}

NT_CHANGE <- function(x){
  x$Change <- apply(x,1,function(x) ifelse(x["REF"] %in% c("A","C"), paste(x["REF"],x["ALT"],sep='>'),paste(comp(x['REF'],forceToLower = F),comp(x['ALT'],forceToLower = F),sep='>')) )
  return(x)
}

LOADING <- function(DATA){
  NUCLEOTIDES <- c('A','C','T','G', 'INS', 'INSo', 'DEL', 'DELo')
  
  DATA <- as.data.frame(fread(DATA))

  DATA$ALT_COUNT <- rowSums(DATA[,NUCLEOTIDES])-DATA$REFf-DATA$REFr
  DATA$DP <- rowSums(DATA[,NUCLEOTIDES])
  DATA$AB <- DATA$ALT_COUNT/DATA$DP_HQ

  #DATA <- DATA[DATA$AB < 0.05 & DATA$DP > 20,]
  DATA <- DATA[!(is.na(DATA$ALT_COUNT)),]

  return(DATA)
}

ERROR_RATE <- function(DATA){
  #DATA$AB <- DATA$ALT_COUNT/DATA$DP
  DATA <- DATA[DATA$DP > 20,]
  DATA <- DATA[!(is.na(DATA$ALT_COUNT)),]
  MEDIAN_DP <- median(DATA$DP)
  SUM_DP <- sum(as.numeric(DATA$DP))
  SUM_ALT <- sum(as.numeric(DATA$ALT_COUNT))    
  ER <- SUM_ALT/SUM_DP
  OUT <- data.frame(ER, nrow(DATA), SUM_ALT, SUM_DP)
  colnames(OUT) <- c('ER', 'Sites', 'ALT_COUNT', 'DP')
  return(OUT)
}

COLLAPSED_ER <- function(DATA){
  DATA[is.na(DATA)] <- 0
  COLLAPSED <- aggregate(cbind(ALT_COUNT,DP_REF) ~ Change, data = DATA, sum)
  COLLAPSED$ER <- COLLAPSED$ALT_COUNT/COLLAPSED$DP_REF
  return(COLLAPSED)
}

give.n <- function(x){
  return(c(y = 0, label = length(x)))
}


## Loading data
# Barcode 1x1
BC1 <- LOADING(bc1)
BC1$Correction_level <- '1x'

# Barcode 1x2
BC2 <- LOADING(bc2)
BC2$Correction_level <- '2x'

# Barcode 2x1
BC3 <- LOADING(bc3)
BC3$Correction_level <- '3x'

# Barcode 2x2
BC4 <- LOADING(bc4)
BC4$Correction_level <- '4x'

## Coverage plots
ALL <- rbind(BC1, BC2, BC3, BC4)

COV_PLOT <- ggplot(ALL, aes(x = Correction_level, y = DP))+
  geom_boxplot(outlier.size = -1, color = "dodgerblue3", fill = "lightskyblue2") +
  scale_fill_grey() +
  stat_summary(aes(x = Correction_level),position=position_dodge(0.90),
               fun.data = give.n, geom = "text",
               vjust = +1, size = 3.5) +
  theme_bw() +
  coord_cartesian(ylim=c(0,max((aggregate(DP ~ Correction_level, data = ALL, function(x) quantile(x,0.75) + 1.5* IQR(x)))$DP)))+
  labs(x='Correction level', y='Depth of coverage', col='Correction level')+
  ggtitle("Depth of coverage per correction level")

## Finding potential germline calls
ALL2 <- aggregate(cbind(DP,ALT_COUNT) ~ CHROM + POS, 
                data = ALL, sum)

ALL2$AB <- ALL2$ALT_COUNT/ALL2$DP
EXCLUDE <- ALL2[(ALL2$AB >= 0.05 & ALL2$DP > 100) | ALL2$AB >= 0.1,]

rm(ALL)
rm(ALL2)

## Error Rates
# Barcode 1
BC1 <- BC1[!paste(BC1$CHROM, BC1$POS, sep = ';') %in% paste(EXCLUDE$CHROM, EXCLUDE$POS, sep = ';') & BC1$DP > 20 & BC1$AB < 0.1,]
ER1 <- ERROR_RATE(BC1)
ER1$Correction_level <- '1x'

ER1_NT <- ER_NT(BC1)
ER1_NT <- NT_CHANGE(ER1_NT)
ER1_NT <- COLLAPSED_ER(ER1_NT)
ER1_NT$Correction_level <- '1x' 

# Barcode 2
BC2 <- BC2[!paste(BC2$CHROM, BC2$POS, sep = ';') %in% paste(EXCLUDE$CHROM, EXCLUDE$POS, sep = ';') & BC2$DP > 20 & BC2$AB < 0.1,]
ER2 <- ERROR_RATE(BC2)
ER2$Correction_level <- '2x'

ER2_NT <- ER_NT(BC2)
ER2_NT <- NT_CHANGE(ER2_NT)
ER2_NT <- COLLAPSED_ER(ER2_NT)
ER2_NT$Correction_level <- '2x' 

# Barcode 3
BC3 <- BC3[!paste(BC3$CHROM, BC3$POS, sep = ';') %in% paste(EXCLUDE$CHROM, EXCLUDE$POS, sep = ';') & BC3$DP > 20 & BC3$AB < 0.1,]
ER3 <- ERROR_RATE(BC3)
ER3$Correction_level <- '3x'

ER3_NT <- ER_NT(BC3)
ER3_NT <- NT_CHANGE(ER3_NT)
ER3_NT <- COLLAPSED_ER(ER3_NT)
ER3_NT$Correction_level <- '3x' 

# Barcode 4
BC4 <- BC4[!paste(BC4$CHROM, BC4$POS, sep = ';') %in% paste(EXCLUDE$CHROM, EXCLUDE$POS, sep = ';') & BC4$DP > 20 & BC4$AB < 0.1,]
ER4 <- ERROR_RATE(BC4)
ER4$Correction_level <- '4x'

ER4_NT <- ER_NT(BC4)
ER4_NT <- NT_CHANGE(ER4_NT)
ER4_NT <- COLLAPSED_ER(ER4_NT)
ER4_NT$Correction_level <- '4x' 


### Getting error rates
ER <- rbind(ER1,ER2,ER3,ER4)
ER$FILTER <- (ER$DP > 1e5)*1

FILTER <- ER[,c('Correction_level', 'FILTER')]

ER_PASS <- ER[ER$FILTER == 1,]
ymax <- round(max(ER_PASS$ER)*1.25,digits = 5)+2E-5
ymax <- ifelse(is.finite(ymax), ymax, 1e-2)
PLOT_ER <- ggplot(ER_PASS,aes(x=Correction_level,y=ER))+
  geom_point(color="blue",size=2)+
  scale_x_discrete(limits=ER$Correction_level)+
  theme_bw()+
  #ylim(min(ER$Error_rate)/5,max(ER$Error_rate)*1.25)+
  scale_y_continuous(breaks = seq(0,ymax, 2.5E-5))+
  ggtitle("Error rate per analysis type")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Error correction type")+
  ylab("Error rate")



### Error rates per nucleotide
### Grouping all together
ER_NT_GROUP <- rbind(ER1_NT,ER2_NT,
                     ER3_NT,ER4_NT)


ER_NT_GROUP$FILTER <- (ER_NT_GROUP$DP_REF > 1e5)*1

ER_NT_GROUP_PASS <- ER_NT_GROUP[ER_NT_GROUP$FILTER == 1,]

ymax <- round(max(ER_NT_GROUP_PASS$ER),digits = 5)+2E-5
ymax <- ifelse(is.finite(ymax), ymax, 1e-2)
PLOT_ER_NT <- ggplot(ER_NT_GROUP_PASS,aes(x=Correction_level,y=ER, color = Change, group = Change))+
  geom_point(size = 2) + geom_line()+
  scale_y_continuous(breaks = seq(0,ymax, 5E-5))+
  scale_x_discrete(limits=ER$Correction_level) +
  theme_bw()+
  ggtitle("Error rate per analysis type")+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  xlab("Error correction level")+
  ylab("Error rate")


### PLOTING INFORMATION ERROR RATES
setwd(OUTF)
ggsave(file="Coverage.pdf",plot=COV_PLOT, useDingbats=FALSE)
ggsave(file="Error_rates_per_nucleotide.pdf",plot=PLOT_ER_NT, useDingbats=FALSE)
ggsave(file='Error_rates.pdf',plot=PLOT_ER, useDingbats=FALSE)
write.table(ER,"Error_rates.txt",row.names = F,col.names = T,sep='\t', quote = F)
write.table(ER_NT_GROUP,"Error_rates_per_nucleotide.txt",row.names = F,col.names = T,sep='\t', quote = F)


# ### COUNT DISTRIBUTION PLOTS
### COUNT DISTRIBUTION PLOTS
BC1_1000 <- rbinom(n = 100000, 1000, ER[ER$Correction_level == '1x',]$ER)

BC2_1000 <- rbinom(n = 100000, 1000, ER[ER$Correction_level == '2x',]$ER)

BC3_1000 <- rbinom(n = 100000, 1000, ER[ER$Correction_level == '3x',]$ER)

BC4_1000 <- rbinom(n = 100000, 1000, ER[ER$Correction_level == '4x',]$ER)

DATA <- data.frame(c(BC1_1000, BC2_1000, BC3_1000, BC4_1000), 
                   c(rep('1x',length(BC1_1000)), rep('2x',length(BC2_1000)), rep('3x',length(BC3_1000)),rep('4x',length(BC4_1000))))

colnames(DATA) <- c('ALT_COUNT_1000', 'Correction_type')

ALT_COUNT_PLOT <- ggplot(DATA, aes(x=as.numeric(ALT_COUNT_1000),  group=Correction_type)) +
  geom_histogram(binwidth = 1, size = 0.25,
                 aes(y = ..density..), color="black", fill="lightblue") +
  facet_grid(~Correction_type) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw()+
  coord_cartesian(xlim=c(-1, 5), ylim=c(0, 1)) +
  #scale_x_discrete(limits=seq(0, 8)) +
  labs(title = "Alt counts per num. of duplicates. Coverage: 1000", x = "Alt count", y = "Fraction", fill = "Alt count")



ALT_COUNT_PLOT_ZOOM <- ggplot(DATA, aes(x=as.numeric(ALT_COUNT_1000),  group=Correction_type)) +
  geom_histogram(binwidth = 1, size = 0.25,
                 aes(y = ..density..), color="black", fill="lightblue") +
  facet_grid(~Correction_type) +
  scale_y_continuous(labels = percent_format()) +
  theme_bw()+
  coord_cartesian(xlim=c(-1, 5), ylim=c(0,0.25)) +
  labs(title = "Alt counts per num. of duplicates. Coverage: 1000", x = "Alt count", y = "Fraction", fill = "Alt count")
  #scale_x_discrete(limits=seq(0, 8))

## SAVING LAST PLOTS
ggsave(file='Alternative_Count.pdf',plot=ALT_COUNT_PLOT,useDingbats=FALSE)
ggsave(file='Alternative_Count_zoom.pdf',plot=ALT_COUNT_PLOT_ZOOM,useDingbats=FALSE)



# # ### COUNT DISTRIBUTION PLOTS
# ### COUNT DISTRIBUTION PLOTS
# BC1$Correction_type <- "1x"
# BC1$ALT_COUNT_1000 <- (BC1$ALT_COUNT/BC1$DP)
# 
# BC2$Correction_type <- "2x"
# BC2$ALT_COUNT_1000 <- (BC2$ALT_COUNT/BC2$DP)
# 
# BC3$Correction_type <- "3x"
# BC3$ALT_COUNT_1000 <- (BC3$ALT_COUNT/BC3$DP)
# 
# BC4$ALT_COUNT_1000 <- (BC4$ALT_COUNT/BC4$DP)
# BC4$Correction_type <- "4x"
# 
# 
# DATA <- rbind(BC1,BC2,BC3,BC4)
# 
# ALT_COUNT_PLOT <- ggplot(DATA, aes(x=as.numeric(ALT_COUNT_1000),  group=Correction_type)) +
#   geom_histogram(binwidth = 0.005, size = 0.25,
#                  aes(y = ..density..), color="black", fill="lightblue") +
#   facet_grid(~Correction_type) +
#   scale_y_continuous(labels = percent_format()) +
#   theme_bw()+
#   coord_cartesian(xlim=c(-0.001, 0.05)) +
#   scale_x_discrete(limits=seq(0, 8)) +
#   labs(title = "Alt counts per num. of duplicates. Coverage: 1000", x = "Alt count", y = "Fraction", fill = "Alt count")
# 
# 
# 
# ALT_COUNT_PLOT_ZOOM <- ggplot(DATA, aes(x=as.numeric(ALT_COUNT_1000),  group=Correction_type)) +
#   geom_histogram(binwidth = 1, size = 0.25,
#                  aes(y = ..density..), color="black", fill="lightblue") +
#   facet_grid(~Correction_type) +
#   scale_y_continuous(labels = percent_format()) +
#   theme_bw()+
#   coord_cartesian(xlim=c(-1, 8), ylim=c(0,0.15)) +
#   scale_x_discrete(limits=seq(0, 8))
#   
# ## SAVING LAST PLOTS
# ggsave(file='Alternative_Count.pdf',plot=ALT_COUNT_PLOT,useDingbats=FALSE)
# ggsave(file='Alternative_Count_zoom.pdf',plot=ALT_COUNT_PLOT_ZOOM,useDingbats=FALSE)
# 
