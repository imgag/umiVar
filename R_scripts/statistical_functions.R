### Deduplicated tsv files are passed as arguments as specified later


## set options
options(warn=-1) #all warnings are ignored
options(scipen=999) #prevents scientific notation
set.seed(07031992)


## suppress start messages of packages
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(survcomp)) # Bioconductor 


#####################################################
# Functions
#####################################################

#--------------------
# Functions to obtain p-value
# Based on parameters alpha and beta provided
# For beta-binomial distribution
#--------------------
### Compute p-value
P_VAL <- function(ALT_COUNT, DP_HQ, a, b){
  if (length(ALT_COUNT) > 0){
    p_val <- pzoibetabinom.ab(ALT_COUNT-1, DP_HQ, a, b, lower.tail=F)
    p_val[p_val < 0] <- 0
    return (p_val)
  } 
}

### Compute p-value of variant
P_VAL_Var <- function(SNP, a, b){
  if (nrow(SNP) > 0){
    p_val <- pzoibetabinom.ab(SNP$VariantCountBias-1, SNP$DP_HQ, a, b, lower.tail=F)
    p_val[p_val < 0] <- 0
    return (p_val)
  }
}

#--------------------
# Functions to estimate the distance between SNVs
#--------------------
# Sub-function
dist_SNP2 <- function(CHROM,POS,CALLS){
  distances <- abs(CALLS$POS - POS)
  distances <- distances[distances > 0] # Ignore the variant we are analysing
  CLOSEST_SNP <- min(distances, na.rm = T) # Get the distance to the closest potential mutation
  return(CLOSEST_SNP)
}

# Main function
dist_SNP <- function(CALLS, SOMATIC_LIKE){
  FINAL <- ''
  NAMES <- names(table(CALLS$CHROM))
  # Computed for each chromosome
  for (chrom in NAMES){
    CALLS2 <- CALLS[CALLS$CHROM == chrom,]
    # Get the coordinates of the potential somatic calls
    SOMATIC_LIKE2 <- SOMATIC_LIKE[SOMATIC_LIKE$CHROM == chrom,] 
    # For each candidate site, get the distance to the closest potential SNV
    CALLS2$DIST <- apply(CALLS2, 1, function(x)  dist_SNP2(x["CHROM"], as.numeric(x["POS"]), SOMATIC_LIKE2))
    # Append the results of the chromosome to a new df
    FINAL <- rbind(FINAL, CALLS2)  
  }
  
  FINAL <- FINAL[!FINAL$CHROM == '',]
  
  return(FINAL)
}

#--------------------
# Function to perform the counting of the reference bases
#--------------------
COUNT <- function(DP1,DP2,DP3,DP4){
  # Sites
  # columns of CRHROM,POS and REF from de-duplication files are stored in unique variables
  unique1 <- unique(DP1[,c('CHROM', 'POS', 'REF')])
  unique2 <- unique(DP2[,c('CHROM', 'POS', 'REF')])
  unique3 <- unique(DP3[,c('CHROM', 'POS', 'REF')])
  unique4 <- unique(DP4[,c('CHROM', 'POS', 'REF')])
  
  
  # Variant covered
  ALL <- (rbind(unique1, unique2, unique3, unique4)) #join (unique) vectors to one matrix ALL
  ALL <- ALL[nchar(ALL$REF) == 1,] # all rows in reference that contain only one character
  ALL <- unique(ALL) # remove duplicated elements/rows
  
  # T and A in nT, C and G in nG
  # Consider in the same count the base and its complementary
  nT <- nrow(ALL[(ALL$REF %in% c('T','A')),])
  nG <- nrow(ALL[(ALL$REF %in% c('G','C')),])
  
  return(c(nT,nG))
}


### TODO: add documentation
COLLAPSE <- function(DATA){
  # Sites with multiple 
  MULTIPLE <- table(paste(DATA$CHROM, DATA$POS, sep = '-')) 
  MULTIPLE <- MULTIPLE[MULTIPLE > 1]
  SITES <- names(MULTIPLE)
  
  # Sites with single calls
  SINGLE <- DATA[!paste(DATA$CHROM, DATA$POS, sep = '-') %in% SITES,]
  
  # Data.frame to save the collapsed variants
  Multiple <- NULL
  
  for (i in SITES){
    # SITE
    CASE <- DATA[paste(DATA$CHROM, DATA$POS, sep = '-') == i,]
    CASE$SUM <- CASE$ALT_COUNT + CASE$Lq_alt_count
    CASE$SUM2 <- rowSums(CASE[,c('P1','P2', 'P3', 'P4')] < 1)
    
    # Determining the most common variants 
    MAX <- CASE[CASE$SUM == max(CASE$SUM),]
    
    # If there are different sites
    if (nrow(MAX) > 1){
      MAX <- MAX[MAX$SUM2 == max(MAX$SUM2),]
    }
    if (nrow(MAX) > 1){
      MAX <- MAX[MAX$Alt_qual == max(MAX$Alt_qual),]
    } 
    if (nrow(MAX) > 1){
      MAX <- MAX[1,]
    }
    
    # Getting values
    NO_MAX <- CASE[!paste(CASE$CHROM, CASE$POS, CASE$REF, CASE$ALT, sep = '-') %in% paste(MAX$CHROM, MAX$POS, MAX$REF, MAX$ALT, sep = '-'),]
    VariantCountBias <- sum(NO_MAX$ALT_COUNT)
    Ns <- sum(NO_MAX$N_count)
    
    DP_HQ <- sum(NO_MAX$DP_HQ) 
    DP <- sum(NO_MAX$DP) 
    Ref_fwd <- sum(NO_MAX$Ref_fwd)
    Ref_rev <- sum(NO_MAX$Ref_rev)
    
    # Recalculating values
    MAX$VariantCountBias <- MAX$VariantCountBias + VariantCountBias
    MAX$DP <- MAX$DP + DP - MAX$N_count - Ns
    MAX$DP_HQ <- MAX$DP_HQ + DP_HQ - MAX$N_count - Ns
    MAX$Ref_fwd <- MAX$Ref_fwd + Ref_fwd
    MAX$Ref_rev <- MAX$Ref_rev + Ref_rev
    
    # Final data.frame
    COMBINE <- MAX[colnames(MAX) %in% colnames(SINGLE)]
    
    Multiple <- rbind(Multiple, COMBINE)
  }
  
  # Merging single and multiple sites
  FINAL <- rbind(SINGLE,Multiple)
  
  return(FINAL)
}


#--------------------
# Main function for beta-binomial computation
# This function is run independently for each DP*
#--------------------
CALLING <- function(DATA, params, ALL_SITES){
  ALT_COLS <- c('A','C','T','G','DEL','INS')
  
  DATA$MM <- rowSums(DATA[,ALT_COLS]) 
  DATA$MM <- DATA$MM - DATA$REFf - DATA$REFr 
  
  CALLS <- DATA[DATA$MM > 0,]
  CALLS_ZERO <- DATA[DATA$MM <= 0,]
  
  # Split sites based on the reference base
  # Different nucleotide changes have different beta-binomial parameters 
  CALLS_A <- CALLS[CALLS$REF == 'A',]
  CALLS_C <- CALLS[CALLS$REF == 'C',]
  CALLS_T <- CALLS[CALLS$REF == 'T',]
  CALLS_G <- CALLS[CALLS$REF == 'G',]
  
  # FOR A in REFERENCE
  if (nrow(CALLS_A) > 0){
    CALLS_A$PA <- NA
    CALLS_A$PC <- P_VAL(CALLS_A$C, CALLS_A$DP_HQ, params[params$Param == 'Alpha','TG'], params[params$Param == 'Beta','TG'])
    CALLS_A$PT <- P_VAL(CALLS_A$T, CALLS_A$DP_HQ, params[params$Param == 'Alpha','TA'], params[params$Param == 'Beta','TA'])
    CALLS_A$PG <- P_VAL(CALLS_A$G, CALLS_A$DP_HQ, params[params$Param == 'Alpha','TC'], params[params$Param == 'Beta','TC'])
    CALLS_A$PINS <- P_VAL(CALLS_A$INS, CALLS_A$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
    CALLS_A$PDEL <- P_VAL(CALLS_A$DEL, CALLS_A$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
    CALLS_A$PINSo <- P_VAL(CALLS_A$INSo, CALLS_A$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
    CALLS_A$PDELo <- P_VAL(CALLS_A$DELo, CALLS_A$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
  } else {
    ## assigns null to each variable if CALLS_A is empty
    CALLS_A$PA <- numeric(0)
    CALLS_A$PC <- numeric(0)
    CALLS_A$PT <- numeric(0)
    CALLS_A$PG <- numeric(0)
    CALLS_A$PINS <- numeric(0)
    CALLS_A$PDEL <- numeric(0)
    CALLS_A$PINSo <- numeric(0)
    CALLS_A$PDELo <- numeric(0) 
  }
  
  # FOR C in REFERENCE
  if (nrow(CALLS_C) > 0){
    CALLS_C$PA <- P_VAL(CALLS_C$A, CALLS_C$DP_HQ, params[params$Param == 'Alpha','GT'], params[params$Param == 'Beta','GT'])
    CALLS_C$PC <- NA
    CALLS_C$PT <- P_VAL(CALLS_C$T, CALLS_C$DP_HQ, params[params$Param == 'Alpha','GA'], params[params$Param == 'Beta','GA'])
    CALLS_C$PG <- P_VAL(CALLS_C$G, CALLS_C$DP_HQ, params[params$Param == 'Alpha','GC'], params[params$Param == 'Beta','GC'])
    CALLS_C$PINS <- P_VAL(CALLS_C$INS, CALLS_C$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
    CALLS_C$PDEL <- P_VAL(CALLS_C$DEL, CALLS_C$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
    CALLS_C$PINSo <- P_VAL(CALLS_C$INSo, CALLS_C$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
    CALLS_C$PDELo <- P_VAL(CALLS_C$DELo, CALLS_C$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
  } else {
    CALLS_C$PA <- numeric(0)
    CALLS_C$PC <- numeric(0)
    CALLS_C$PT <- numeric(0)
    CALLS_C$PG <- numeric(0)
    CALLS_C$PINS <- numeric(0)
    CALLS_C$PDEL <- numeric(0)
    CALLS_C$PINSo <- numeric(0)
    CALLS_C$PDELo <- numeric(0) 
  }
  # FOR T in REFERENCE
  if(nrow(CALLS_T) > 0){
    CALLS_T$PA <- P_VAL(CALLS_T$A, CALLS_T$DP_HQ, params[params$Param == 'Alpha','TA'], params[params$Param == 'Beta','TA'])
    CALLS_T$PC <- P_VAL(CALLS_T$C, CALLS_T$DP_HQ, params[params$Param == 'Alpha','TC'], params[params$Param == 'Beta','TC'])
    CALLS_T$PT <- NA
    CALLS_T$PG <- P_VAL(CALLS_T$G, CALLS_T$DP_HQ, params[params$Param == 'Alpha','TG'], params[params$Param == 'Beta','TG'])
    CALLS_T$PINS <- P_VAL(CALLS_T$INS, CALLS_T$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
    CALLS_T$PDEL <- P_VAL(CALLS_T$DEL, CALLS_T$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
    CALLS_T$PINSo <- P_VAL(CALLS_T$INSo, CALLS_T$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
    CALLS_T$PDELo <- P_VAL(CALLS_T$DELo, CALLS_T$DP_HQ, params[params$Param == 'Alpha','Indel_A'], params[params$Param == 'Beta','Indel_A'])
  } else {
    CALLS_T$PA <- numeric(0)
    CALLS_T$PC <- numeric(0)
    CALLS_T$PT <- numeric(0)
    CALLS_T$PG <- numeric(0)
    CALLS_T$PINS <- numeric(0)
    CALLS_T$PDEL <- numeric(0)
    CALLS_T$PINSo <- numeric(0)
    CALLS_T$PDELo <- numeric(0) 
  }
  
  # FOR G in REFERENCE
  if (nrow(CALLS_G) > 0){
    CALLS_G$PA <- P_VAL(CALLS_G$A, CALLS_G$DP_HQ, params[params$Param == 'Alpha','GA'], params[params$Param == 'Beta','GA'])
    CALLS_G$PC <- P_VAL(CALLS_G$C, CALLS_G$DP_HQ, params[params$Param == 'Alpha','GC'], params[params$Param == 'Beta','GC'])
    CALLS_G$PT <- P_VAL(CALLS_G$T, CALLS_G$DP_HQ, params[params$Param == 'Alpha','GT'], params[params$Param == 'Beta','GT'])
    CALLS_G$PG <- NA
    CALLS_G$PINS <- P_VAL(CALLS_G$INS, CALLS_G$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
    CALLS_G$PDEL <- P_VAL(CALLS_G$DEL, CALLS_G$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
    CALLS_G$PINSo <- P_VAL(CALLS_G$INSo, CALLS_G$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
    CALLS_G$PDELo <- P_VAL(CALLS_G$DELo, CALLS_G$DP_HQ, params[params$Param == 'Alpha','Indel_G'], params[params$Param == 'Beta','Indel_G'])
  } else {
    CALLS_G$PA <- numeric(0)
    CALLS_G$PC <- numeric(0)
    CALLS_G$PT <- numeric(0)
    CALLS_G$PG <- numeric(0)
    CALLS_G$PINS <- numeric(0)
    CALLS_G$PDEL <- numeric(0)
    CALLS_G$PINSo <- numeric(0)
    CALLS_G$PDELo <- numeric(0) 
  }
  
  
  # Filling zero table
  if (nrow(CALLS_ZERO) > 0){
    CALLS_ZERO$PA <- 1
    CALLS_ZERO$PC <- 1
    CALLS_ZERO$PT <- 1
    CALLS_ZERO$PG <- 1
    CALLS_ZERO$PINS <- 1
    CALLS_ZERO$PDEL <- 1
    CALLS_ZERO$PINSo <- 1
    CALLS_ZERO$PDELo <- 1
  } else {
    CALLS_ZERO$PA <- numeric(0)
    CALLS_ZERO$PC <- numeric(0)
    CALLS_ZERO$PT <- numeric(0)
    CALLS_ZERO$PG <- numeric(0)
    CALLS_ZERO$PINS <- numeric(0)
    CALLS_ZERO$PDEL <- numeric(0)
    CALLS_ZERO$PINSo <- numeric(0)
    CALLS_ZERO$PDELo <- numeric(0)
  }
  
  VARIANTS <- rbind(CALLS_A, CALLS_C, CALLS_T, CALLS_G, CALLS_ZERO)
  
  return (VARIANTS)
}

#--------------------
# Functions to estimate beta-binomial parameters
#--------------------
# 1. Main alpha-beta parameters estimation function
# Alternative counts versus number of trials (high-quality coverage = DP_HQ)
MODEL_BETABIN_DP_DP <- function(DATA){
  result <- tryCatch({mle2(ALT_COUNT~dbetabinom.ab(size=DP_HQ,shape1,shape2), ## mle2: maximum likelihood estimation n=1
                           data=DATA,
                           method="Nelder-Mead", #optimization method
                           skip.hessian=TRUE, ##bypasses hessian calculation
                           start=list(shape1=1,shape2=round(mean(DATA$DP_HQ))),
                           control=list(maxit=1000))}, 
                     error = function(e) {(estBetaParams(mean(DATA$ALT_COUNT/DATA$DP_HQ, na.rm = T),var(DATA$ALT_COUNT/DATA$DP_HQ, na.rm = T)))})
  PARAM1 <- ifelse(is.null(coef(result)[[1]]),result[[1]], coef(result)[[1]])
  PARAM2 <- ifelse(is.null(coef(result)[[2]]),result[[2]], coef(result)[[2]])
  
  return (c(PARAM1, PARAM2))  
} 

# 2. Alternative alpha-beta parameters estimation function
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#---------------------------
# Function to perform the last step of the variant calling
# Merge all results and obtain values for the vcf contruction
#---------------------------
CALL <- function(LINE, INDEL, num_sites, Ps){
  alt <- c("A", "C", "T", "G", "INS", "DEL")
  P <- c("PA", "PC", "PT", "PG", "PINS", "PDEL", "PINSo", "PDELo")
  P2 <- c("PA", "PC", "PT", "PG", "PINS", "PDEL")
  P_adj <- c("PA_adj", "PC_adj", "PT_adj", "PG_adj", "PINS_adj", "PDEL_adj")
  FOR <- c("Af", "Cf", "Tf", "Gf", "INSf", "DELf")
  REV <- c("Ar", "Cr", "Tr", "Gr", "INSr", "DELr")
  ALL_id <- seq(1,length(P))
  
  num_sites <- max(as.numeric(num_sites), 0)
  
  # Reference info
  REF <- as.character(LINE['REF'])
  REFf <- as.numeric(LINE['REFf'])
  REFr <- as.numeric(LINE['REFr'])
  
  
  ALL_ALT_COUNT <- c('A','C','T','G', 'INS', 'DEL', 'INSo','DELo')
  P_ADJ <- as.numeric(LINE[P_adj]) # Extract adjusted p-values from the LINE
  p2 <- as.numeric(LINE[P2]) # Extract p-values from the LINE
  
  # What indexes have an adjusted p-value < 0.05
  # In next steps, we will use these indexes for several computations
  id <- which(P_ADJ <= 0.05)
  
  # Taking always one, even though none of them is significant (length(id) == 0)
  if (length(id) == 0){
    id <- which(p2 == min(p2, na.rm = T))[[1]]
  }
  
  # If more than one position significant
  if (length(id) > 0){
    call <- NULL
    for (i in id){
      ALT <- alt[i]
      
      # Get genotypes
      if (ALT %in% c('DEL', 'INS')){
        GENOTYPE <- INDEL[paste(LINE['CHROM'], as.numeric(LINE['POS']), sep = '-') == paste(INDEL$CHROM, INDEL$POS, sep = '-') & 
                            INDEL$Type == ALT,]$Genotype
      } else {
        GENOTYPE <- paste(REF, ALT, sep = '>')
      }
      
      # Different values for constructing the info for VCF
      Ps <- as.numeric(LINE[P]) # All probabilities
      FORWARD <- as.numeric(LINE[FOR[i]]) # Forward-read counts
      REVERSE <- as.numeric(LINE[REV[i]]) # Reverse-read counts
      ALT_COUNT <- as.numeric(FORWARD+REVERSE) # Sum: Fwd+Rev
      ALT_COUNT_all <- as.numeric(LINE[ALL_ALT_COUNT]) # All allele counts
      ALT_COUNT_padj <- formatC(as.numeric(P_ADJ[i]), format = "e", digits = 2)
      ALT_COUNT_p <- formatC(as.numeric(Ps[i]), format = "e", digits = 2)
      
      # Fisher strand values values
      FISHER <- (fisher.test(matrix(c(FORWARD, REVERSE, REFf, REFr), nrow = 2))[[1]])
      FISHER_adj <- formatC(p.adjust(FISHER, n = num_sites, method = 'fdr'), format = "e", digits = 2)
      
      # Label variant as significant or not, based on the corrected p-value (0.1)
      if (as.numeric(ALT_COUNT_padj) < 0.1){
        FILTER <- 'SIG'
      } else {
        FILTER <- 'NO_SIG'
      }
      
      # Calculate Allele balance
      AB <- round(ALT_COUNT/as.numeric(LINE['DP_HQ']), 6)
      
      # Obtain the alternative counts and probabilities
      # Other significant alternative alleles, such as scenarios found in biallelic sites, are not considered for this variable
      id_o <- ALL_id[!ALL_id %in% id]
      Ps_out <- Ps[ALL_id[!ALL_id %in% id]]
      OUT <- (prod(Ps_out, na.rm = T))
      OUT_adj <- formatC(p.adjust(c(OUT, Ps), method = 'fdr')[1], format = "e", digits = 2)
      ALT_COUNT_o <- as.numeric(sum(ALT_COUNT_all[id_o]) - as.numeric(LINE['REFf']) - as.numeric(LINE['REFr']))
      
      # Collapse all information in a string variable
      call_temp <- paste(GENOTYPE, FILTER, ALT_COUNT, AB, ALT_COUNT_p, ALT_COUNT_padj, paste(FORWARD, REVERSE, REFf, REFr, sep = '-'), FISHER_adj, ALT_COUNT_o, OUT_adj, sep = ':')
      
      # Append to all possible somatic calls
      call <- c(call, call_temp)
    }
    call <- paste(call, collapse = '|')
    
  } else {
    call <- 'No_call'
  }
  return(call)
}

#####################################################
#####################################################


###################
# Arguments
###################

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-t1", "--tumor_file1", type="character", help="Tumor - Read count input file for duplicates = 1", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-t2", "--tumor_file2", type="character", help="Tumor - Read count input file for duplicates = 2", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-t3", "--tumor_file3", type="character", help="Tumor - Read count input file for duplicates = 3", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-t4", "--tumor_file4", type="character", help="Tumor - Read count input file for duplicates = 4", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-n", "--normal_file", type="character", help="Normal - Read count input file", nargs=1, default = NULL)
parser$add_argument("-o", "--out_file", type="character", help="output_file", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-params", "--parameters", type="character", help="Table with parameters", nargs=1, required=TRUE)
parser$add_argument("-num", "--num_sites", default = NULL, type="character", help="Total number of sites analysed. Important for the p-value correction", nargs=1, required=FALSE)

# Reading parameters
args <- parser$parse_args()

# Assign names to files
TD1 <- args$tumor_file1
TD2 <- args$tumor_file2
TD3 <- args$tumor_file3
TD4 <- args$tumor_file4
parameters <- args$parameters
num_sites <- as.numeric(args$num_sites)

ND <- args$normal_file

OUTF <- args$out_file

# print (paste("Normal file:", ND,sep=" "))
# print (paste("Tumor file (DP1):", TD1, sep=" "))
# print (paste("Tumor file (DP2):", TD2, sep=" "))
# print (paste("Tumor file (DP3):", TD3, sep=" "))
# print (paste("Tumor file (DP4):", TD4, sep=" "))
# print (paste("Out file:", OUTF))

#------------------
# Running
#------------------

# Pass different tumor files to DP variables
# as.data.frame:check if an object is a data frame, if not coerce it
# Fread for regular delimited files
DP1 <- as.data.frame(fread(TD1)) 
DP2 <- as.data.frame(fread(TD2))
DP3 <- as.data.frame(fread(TD3))
DP4 <- as.data.frame(fread(TD4))
PARAMS <- as.data.frame(fread(parameters))

# Get total number of unique potential somatic coordinates to be analysed
ALT_COLS <- c('A','C','T','G','DEL','INS')
ALL_SITES <- rbind(DP1, DP2, DP3, DP4)
ALL_SITES$MM <- rowSums(ALL_SITES[,ALT_COLS]) 
ALL_SITES$MM <- ALL_SITES$MM - ALL_SITES$REFf - ALL_SITES$REFr 
ALL_SITES <- ALL_SITES[ALL_SITES$MM > 0,]
ALL_SITES <- unique(ALL_SITES[,c('CHROM', 'POS')])

## perform variant calling
# DP1
CALLING1 <- CALLING(DP1, PARAMS[PARAMS$BARCODE == 'DP1',], ALL_SITES)
rm (DP1)

# DP2
CALLING2 <- CALLING(DP2, PARAMS[PARAMS$BARCODE == 'DP2',], ALL_SITES)
rm(DP2)

# DP3
CALLING3 <- CALLING(DP3, PARAMS[PARAMS$BARCODE == 'DP3',], ALL_SITES)
rm(DP3)

# DP4
CALLING4 <- CALLING(DP4, PARAMS[PARAMS$BARCODE == 'DP4',], ALL_SITES)
rm(DP4)

# Merging results
# Merge results obtained previously
if (sum(c(nrow(CALLING1),nrow(CALLING2),nrow(CALLING3),nrow(CALLING4))) > 0){
  VARIANT0 <- aggregate(cbind(DP, DP_HQ, REFf, REFr, A, Af, Ar, Aq, Ao, C, Cf, Cr, Cq, Co, `T`, Tf, Tr, Tq, To,
                              G, Gf, Gr, Gq, Go, INS, INSf, INSr, INSo, DEL, DELf, DELr, DELo) ~ 
                          CHROM+POS+REF, 
                        rbind(CALLING1, CALLING2, CALLING3, CALLING4), 
                        FUN = sum, na.rm = F, na.action= NULL)
  
  LQ <- aggregate(cbind(Lq_alt_count,MM) ~ 
                    CHROM+POS+REF, 
                  rbind(CALLING1, CALLING2, CALLING3, CALLING4), 
                  FUN = median, na.rm = F, na.action= NULL)
  LQ$Ratio <- LQ$Lq_alt_count/(LQ$Lq_alt_count+LQ$MM)
  LQ <- LQ[,c('CHROM', 'POS', 'REF', 'Ratio')]
  
  colnames(LQ) <- c('CHROM', 'POS', 'REF', 'Lq_alt_count')
  
  VARIANT0 <- merge(VARIANT0, LQ, by = c('CHROM', 'POS', 'REF'), all.x = T)
} else {
  VARIANT0 <- rbind(CALLING1, CALLING2, CALLING3, CALLING4)
}


# Duplicate1
# Rename colnames to merged them aftewards
P_VAL1 <- CALLING1[,c('CHROM', 'POS', 'REF', 'PA','PC', 'PT', 'PG', 'PINS', 'PDEL', 'PINSo', 'PDELo')]
colnames(P_VAL1) <- c('CHROM', 'POS', 'REF', 'PA1','PC1', 'PT1', 'PG1', 'PINS1', 'PDEL1', 'PINSo1', 'PDELo1')


# Duplicate2
# Rename colnames to merged them aftewards
P_VAL2 <- CALLING2[,c('CHROM', 'POS', 'REF', 'PA','PC', 'PT', 'PG', 'PINS', 'PDEL', 'PINSo', 'PDELo')]
colnames(P_VAL2) <- c('CHROM', 'POS', 'REF', 'PA2','PC2', 'PT2', 'PG2', 'PINS2', 'PDEL2', 'PINSo2', 'PDELo2')


# Duplicate3
# Rename colnames to merged them aftewards
P_VAL3 <- CALLING3[,c('CHROM', 'POS', 'REF', 'PA','PC', 'PT', 'PG', 'PINS', 'PDEL', 'PINSo', 'PDELo')]
colnames(P_VAL3) <- c('CHROM', 'POS', 'REF', 'PA3','PC3', 'PT3', 'PG3', 'PINS3', 'PDEL3', 'PINSo3', 'PDELo3')


# Duplicate4
# Rename colnames to merged them aftewards
P_VAL4 <- CALLING4[,c('CHROM', 'POS', 'REF', 'PA','PC', 'PT', 'PG', 'PINS', 'PDEL', 'PINSo', 'PDELo')]
colnames(P_VAL4) <- c('CHROM', 'POS', 'REF', 'PA4','PC4', 'PT4', 'PG4', 'PINS4', 'PDEL4', 'PINSo4', 'PDELo4')


# Merging results
# Append columns 
VARIANT1 <- merge(VARIANT0, P_VAL1, by = c('CHROM', 'POS', 'REF'), all.x = T)
VARIANT2 <- merge(VARIANT1, P_VAL2, by = c('CHROM', 'POS', 'REF'), all.x = T)
VARIANT3 <- merge(VARIANT2, P_VAL3, by = c('CHROM', 'POS', 'REF'), all.x = T)
VARIANT4 <- merge(VARIANT3, P_VAL4, by = c('CHROM', 'POS', 'REF'), all.x = T)

rm(VARIANT1)
rm(VARIANT2)
rm(VARIANT3)

# Change NA by "1"
VARIANT4[is.na(VARIANT4)] <- 1

# Merging p-values by fisher test
if (nrow(VARIANT4) > 0){
  VARIANT4$PA <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'A', NA, 
                                                       combine.test(c(as.numeric(x['PA1']),as.numeric(x['PA2']), as.numeric(x['PA3']), as.numeric(x['PA4'])), 
                                                                    method = 'fisher')))
  VARIANT4$PC <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'C', NA, 
                                                       combine.test(c(as.numeric(x['PC1']),as.numeric(x['PC2']), as.numeric(x['PC3']), as.numeric(x['PC4'])), 
                                                                    method = 'fisher')))
  VARIANT4$PT <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'T', NA, 
                                                       combine.test(c(as.numeric(x['PT1']),as.numeric(x['PT2']), as.numeric(x['PT3']), as.numeric(x['PT4'])), 
                                                                    method = 'fisher'))) 
  VARIANT4$PG <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'G', NA, 
                                                       combine.test(c(as.numeric(x['PG1']),as.numeric(x['PG2']), as.numeric(x['PG3']), as.numeric(x['PG4'])), 
                                                                    method = 'fisher')))
  VARIANT4$PINS <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'IND', NA, 
                                                         combine.test(c(as.numeric(x['PINS1']),as.numeric(x['PINS2']), as.numeric(x['PINS3']), as.numeric(x['PINS4'])), 
                                                                      method = 'fisher')))
  VARIANT4$PDEL <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'IND', NA, 
                                                         combine.test(c(as.numeric(x['PDEL1']),as.numeric(x['PDEL2']), as.numeric(x['PDEL3']), as.numeric(x['PDEL4'])), 
                                                                      method = 'fisher')))
  VARIANT4$PINSo <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'IND', NA, 
                                                          combine.test(c(as.numeric(x['PINSo1']),as.numeric(x['PINSo2']), as.numeric(x['PINSo3']), as.numeric(x['PINSo4'])), 
                                                                       method = 'fisher')))
  VARIANT4$PDELo <- apply(VARIANT4, 1, function(x) ifelse(x['REF'] == 'IND', NA, 
                                                          combine.test(c(as.numeric(x['PDELo1']),as.numeric(x['PDELo2']), as.numeric(x['PDELo3']), as.numeric(x['PDELo4'])), 
                                                                       method = 'fisher')))
}


# Collapsed calls
VARIANTS <- VARIANT4
rm(VARIANT4)

# P-value correction (FDR method)
if (nrow(VARIANTS) > 0){
  if (is.null(num_sites) || num_sites <= 0){
    counts <- nrow(ALL_SITES)
    
    VARIANTS$PA_adj <- p.adjust(VARIANTS$PA, n = max(counts, length(VARIANTS$PA)),method = 'fdr')
    VARIANTS$PC_adj <- p.adjust(VARIANTS$PC, n = max(counts, length(VARIANTS$PC)),method = 'fdr')
    VARIANTS$PT_adj <- p.adjust(VARIANTS$PT, n = max(counts, length(VARIANTS$PT)),method = 'fdr')
    VARIANTS$PG_adj <- p.adjust(VARIANTS$PG, n = max(counts, length(VARIANTS$PG)),method = 'fdr')
    VARIANTS$PINS_adj <- p.adjust(VARIANTS$PINS, n = max(counts, length(VARIANTS$PINS)),method = 'fdr')
    VARIANTS$PDEL_adj <- p.adjust(VARIANTS$PDEL, n = max(counts, length(VARIANTS$PDEL)),method = 'fdr')
    VARIANTS$PINSo_adj <- p.adjust(VARIANTS$PINSo, n = max(counts, length(VARIANTS$PINSo)),method = 'fdr')
    VARIANTS$PDELo_adj <- p.adjust(VARIANTS$PDELo, n = max(counts, length(VARIANTS$PDELo)),method = 'fdr')
  } else {
    VARIANTS$PA_adj <- p.adjust(VARIANTS$PA, n = max(num_sites, length(VARIANTS$PA)),method = 'fdr')
    VARIANTS$PC_adj <- p.adjust(VARIANTS$PC, n = max(num_sites, length(VARIANTS$PC)),method = 'fdr')
    VARIANTS$PT_adj <- p.adjust(VARIANTS$PT, n = max(num_sites, length(VARIANTS$PT)),method = 'fdr')
    VARIANTS$PG_adj <- p.adjust(VARIANTS$PG, n = max(num_sites, length(VARIANTS$PG)),method = 'fdr')
    VARIANTS$PINS_adj <- p.adjust(VARIANTS$PINS, n = max(num_sites, length(VARIANTS$PINS)),method = 'fdr')
    VARIANTS$PDEL_adj <- p.adjust(VARIANTS$PDEL, n = max(num_sites, length(VARIANTS$PDEL)),method = 'fdr')
    VARIANTS$PINSo_adj <- p.adjust(VARIANTS$PINSo, n = max(num_sites, length(VARIANTS$PINSo)),method = 'fdr')
    VARIANTS$PDELo_adj <- p.adjust(VARIANTS$PDELo, n = max(num_sites, length(VARIANTS$PDELo)),method = 'fdr')
  }
} else {
  VARIANTS$PA_adj <- numeric(0)
  VARIANTS$PC_adj <- numeric(0)
  VARIANTS$PT_adj <- numeric(0)
  VARIANTS$PG_adj <- numeric(0)
  VARIANTS$PINS_adj <- numeric(0)
  VARIANTS$PDEL_adj <- numeric(0)
  VARIANTS$PINSo_adj <- numeric(0)
  VARIANTS$PDELo_adj <- numeric(0)
  
}


## Decide genotypes
#-----------------
## Indel analysis
#-----------------
ALL <- rbind(CALLING1, CALLING2, CALLING3, CALLING4)
ALL$COUNT <- 1

# Deletion
ALL_del <- ALL[ALL$DEL > 0,]
if (nrow(ALL_del) > 0){
  Del <- aggregate(cbind(DEL, COUNT) ~ CHROM + POS + REF + DELg, ALL_del, sum)
  
  # first decide for the one more represented in more duplicate groups
  Del_max_count <- aggregate(COUNT ~ CHROM + POS, Del, max)
  colnames(Del_max_count) <- c('CHROM','POS','Max_count')
  
  Del <- merge(Del, Del_max_count, by = c('CHROM', 'POS'))
  Del <- Del[Del$COUNT == Del$Max_count,]
  
  # second choose for more represented variant
  Del_max_del <- aggregate(DEL ~ CHROM + POS, Del, max)
  colnames(Del_max_del) <- c('CHROM','POS','Max_del')
  
  Del <- merge(Del, Del_max_del, by = c('CHROM', 'POS'))
  Del <- Del[Del$DEL == Del$Max_del,]  
  
  Del <- Del[,c('CHROM', 'POS','DELg')]
  colnames(Del) <- c('CHROM', 'POS','Genotype')
  Del$Type <- 'DEL'
  
} else {
  Del <- numeric(0)
}

# Insertion
ALL_ins <- ALL[ALL$INS > 0,]
if (nrow(ALL_ins) > 0){
  Ins <- aggregate(cbind(INS, COUNT) ~ CHROM + POS + REF + INSg, ALL_ins, sum)
  
  # first decide for the one more represented in more duplicate groups
  Ins_max_count <- aggregate(COUNT ~ CHROM + POS, Ins, max)
  colnames(Ins_max_count) <- c('CHROM','POS','Max_count')
  
  Ins <- merge(Ins, Ins_max_count, by = c('CHROM', 'POS'))
  Ins <- Ins[Ins$COUNT == Ins$Max_count,]
  
  # second choose for more represented variant
  Ins_max_Ins <- aggregate(INS ~ CHROM + POS, Ins, max)
  colnames(Ins_max_Ins) <- c('CHROM','POS','Max_Ins')
  
  Ins <- merge(Ins, Ins_max_Ins, by = c('CHROM', 'POS'))
  Ins <- Ins[Ins$INS == Ins$Max_Ins,]  
  
  Ins <- Ins[,c('CHROM', 'POS','INSg')]
  colnames(Ins) <- c('CHROM', 'POS','Genotype')
  Ins$Type <- 'INS'
  
} else {
  Ins <- numeric(0)
}

INDEL <- rbind(Ins, Del)

rm(CALLING1)
rm(CALLING2)
rm(CALLING3)
rm(CALLING4)

## Get the calls for each site
VARIANTS$DP_HQ <- rowSums(VARIANTS[,c('A','C','T','G', 'INS', 'DEL', 'INSo','DELo')])
VARIANTS$MM <- VARIANTS$DP_HQ - VARIANTS$REFf - VARIANTS$REFr

# Take p-values for correction, remove NAs and random select num_sites
if (nrow(VARIANTS) > 0){
  Ps <- c(VARIANTS$PA, VARIANTS$PG, VARIANTS$PC, VARIANTS$PT)
  Ps <- sample(Ps, size = max(num_sites, nrow(VARIANTS)))
  VARIANTS$CALL <- apply(VARIANTS, 1, function(x) CALL(x, INDEL, num_sites = max(num_sites, nrow(VARIANTS)), Ps)) 
} else {
  VARIANTS$FORMAT <- numeric(0)
  VARIANTS$CALL <- numeric(0)
}

# Final calls
FINAL <- VARIANTS[,c('CHROM','POS','REF','DP','DP_HQ', 'REFf', 'REFr','MM','CALL')]
FINAL$POS <- as.numeric(as.character(FINAL$POS))

# Likely somatic
if (nrow(FINAL) > 0){
  SOMATIC_LIKE <- FINAL[!grepl('NO_SIG', FINAL$CALL) & grepl('SIG', FINAL$CALL),]
  # Distance between variants
  FINAL <- dist_SNP(CALLS = FINAL, SOMATIC_LIKE = SOMATIC_LIKE)
} else {
  FINAL$DIST <- numeric(0)
}

# Ordering tsv file
FINAL <- FINAL[with(FINAL, order(CHROM, POS)),]

# Printing final table
FORMAT <- paste('##CALL COLUMN = GT', 'SIG', 'ALT_COUNT', 'AB', 'ALT_COUNT_p', 'ALT_COUNT_padj', 'ALTf-ATLr-REFf-REFr', 'FISHER_adj', 'ALT_COUNT_o', 'OUT_adj\n', sep = ':')
cat(FORMAT, file = OUTF)
write.table(x = FINAL, file = OUTF, quote = F, sep = '\t', col.names = T, row.names = F, append = T)
