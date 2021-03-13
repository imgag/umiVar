options(warn=-1)
options(scipen=999)

set.seed(07031992)

suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))



# Variant calling functions
# FUNCTION

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

MODEL_BETABIN_DP_DP <- function(DATA){
  result <- tryCatch({mle2(ALT_COUNT~dbetabinom.ab(size=DP_HQ,shape1,shape2),
                           data=DATA,
                           method="Nelder-Mead",
                           skip.hessian=TRUE,
                           start=list(shape1=1,shape2=round(mean(DATA$DP_HQ))),
                           control=list(maxit=1000))}, 
                     error = function(e) {(estBetaParams(mean(DATA$ALT_COUNT/DATA$DP_HQ, na.rm = T),var(DATA$ALT_COUNT/DATA$DP_HQ, na.rm = T)))})
  PARAM1 <- ifelse(is.null(coef(result)[[1]]),result[[1]], coef(result)[[1]])
  PARAM2 <- ifelse(is.null(coef(result)[[2]]),result[[2]], coef(result)[[2]])
  
  return (c(PARAM1, PARAM2))  
} 

SAMPLE <- function(TD){
  NUCLEOTIDES <- c('A','C','T','G', 'INS', 'INSo', 'DEL', 'DELo')
  
  # Training set
  CONTROL <- TD
  CONTROL$ALT_COUNT <- rowSums(CONTROL[,NUCLEOTIDES])-CONTROL$REFf-CONTROL$REFr
  CONTROL$DP_HQ <- rowSums(CONTROL[,NUCLEOTIDES])
  CONTROL$AB <- CONTROL$ALT_COUNT/CONTROL$DP_HQ
  CONTROL <- CONTROL[CONTROL$AB < 0.1 & CONTROL$DP_HQ >= 10 & CONTROL$ALT_COUNT >= 0,] # Focusing only to non-germline sites
  
  if (nrow(CONTROL) > 500000) {
    CONTROL <- CONTROL[sample(nrow(CONTROL), 500000),]
  } else {
    CONTROL <- CONTROL[sample(nrow(CONTROL), 500000, replace = T),]
  }
  
  CONTROL <- CONTROL[CONTROL$AB < 0.1 & CONTROL$DP_HQ >= 10,]
  
  return(CONTROL)
}

PARAMS <- function(CONTROL){
  ## ONLY INCLUDING VARIANTS BUT INDELS SEPARATED
  # Calculate indels
  Indel <- c('INS', 'INSo', 'DEL', 'DELo')
  CONTROL$Indel <- rowSums(CONTROL[,Indel])
  
  # Reference bases
  CONTROL_T <- CONTROL[CONTROL$REF == 'T',] 
  CONTROL_A <- CONTROL[CONTROL$REF == 'A',]
  
  CONTROL_G <- CONTROL[CONTROL$REF == 'G',] 
  CONTROL_C <- CONTROL[CONTROL$REF == 'C',] 
  
  COLNAMES <- c('DP_HQ','ALT_COUNT')
  
  # T>G | A>C
  TG <- CONTROL_T[,c('DP_HQ','G')]
  colnames(TG) <- COLNAMES
  
  AC <- CONTROL_A[,c('DP_HQ','C')]
  colnames(AC) <- COLNAMES
  
  CONTROL_TG <- rbind(TG, AC)
  
  # T>A | A>T
  TA <- CONTROL_T[,c('DP_HQ','A')]
  colnames(TA) <- COLNAMES
  
  AT <- CONTROL_A[,c('DP_HQ','T')]
  colnames(AT) <- COLNAMES
  
  CONTROL_TA <- rbind(TA, AT)
  
  # T>C | A>G
  TC <- CONTROL_T[,c('DP_HQ','C')]
  colnames(TC) <- COLNAMES
  
  AG <- CONTROL_A[,c('DP_HQ','G')]
  colnames(AG) <- COLNAMES
  
  CONTROL_TC <- rbind(TC, AG) 
  
  # G>T | C>A
  GT <- CONTROL_G[,c('DP_HQ','T')]
  colnames(GT) <- COLNAMES
  
  CA <- CONTROL_C[,c('DP_HQ','A')]
  colnames(CA) <- COLNAMES
  
  CONTROL_GT <- rbind(GT, CA) 
  
  # G>A | C>T
  GA <- CONTROL_G[,c('DP_HQ','A')]
  colnames(GA) <- COLNAMES
  
  CT <- CONTROL_C[,c('DP_HQ','T')]
  colnames(CT) <- COLNAMES
  
  CONTROL_GA <- rbind(GA, CT) 
  
  # G>C | C>G
  GC <- CONTROL_G[,c('DP_HQ','C')]
  colnames(GC) <- COLNAMES
  
  CG <- CONTROL_C[,c('DP_HQ','G')]
  colnames(CG) <- COLNAMES
  
  CONTROL_GC <- rbind(GC, CG)   
  
  ## For indels
  # A>Indel | T>Indel
  AIndel <- CONTROL_A[,c('DP_HQ','Indel')]
  colnames(AIndel) <- COLNAMES
  
  TIndel <- CONTROL_T[,c('DP_HQ','Indel')]
  colnames(TIndel) <- COLNAMES
  
  CONTROL_Indel_A <- rbind(AIndel, TIndel)
  
  # G>Indel | C>Indel
  GIndel <- CONTROL_G[,c('DP_HQ','Indel')]
  colnames(GIndel) <- COLNAMES
  
  CIndel <- CONTROL_C[,c('DP_HQ','Indel')]
  colnames(CIndel) <- COLNAMES
  
  CONTROL_Indel_G <- rbind(GIndel, CIndel)  
  
  
  ## Estimating parameters for each nucleotide change when we only have 1 duplicate per barcode group
  FIT_DP.TG <- MODEL_BETABIN_DP_DP(CONTROL_TG)
  
  FIT_DP.TA <- MODEL_BETABIN_DP_DP(CONTROL_TA)
  
  FIT_DP.TC <- MODEL_BETABIN_DP_DP(CONTROL_TC)
  
  FIT_DP.GT <- MODEL_BETABIN_DP_DP(CONTROL_GT)
  
  FIT_DP.GA <- MODEL_BETABIN_DP_DP(CONTROL_GA)
  
  FIT_DP.GC <- MODEL_BETABIN_DP_DP(CONTROL_GC)
  
  FIT_DP.Indel_A <- MODEL_BETABIN_DP_DP(CONTROL_Indel_A)
  
  FIT_DP.Indel_G <- MODEL_BETABIN_DP_DP(CONTROL_Indel_G)
  
  TABLE <- data.frame(FIT_DP.TG, FIT_DP.TA, FIT_DP.TC, FIT_DP.GT, FIT_DP.GA, FIT_DP.GC, FIT_DP.Indel_A, FIT_DP.Indel_G)
  colnames(TABLE) <- c('TG', 'TA', 'TC', 'GT', 'GA', 'GC', 'Indel_A', 'Indel_G')
  TABLE$Param <- c('Alpha', 'Beta') 
  TABLE <- TABLE[,c('Param','TG', 'TA', 'TC', 'GT', 'GA', 'GC', 'Indel_A', 'Indel_G')]
  rownames(TABLE) <- TABLE$Param
  
  return(TABLE)
}

###################
# Arguments
###################

parser <- ArgumentParser()

# setting parameters
parser$add_argument("-t1", "--tumor_file1", type="character", help="Tumor - Read count input file for duplicates = 1", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-t2", "--tumor_file2", type="character", help="Tumor - Read count input file for duplicates = 2", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-t3", "--tumor_file3", type="character", help="Tumor - Read count input file for duplicates = 3", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-t4", "--tumor_file4", type="character", help="Tumor - Read count input file for duplicates = 4", metavar="file", nargs=1, required=TRUE)
parser$add_argument("-o", "--out_file", type="character", help="output_file", metavar="file", nargs=1, required=TRUE)

# reading parameters
args <- parser$parse_args()

TD1 <- args$tumor_file1
TD2 <- args$tumor_file2
TD3 <- args$tumor_file3
TD4 <- args$tumor_file4

OUTF <- args$out_file

###################
# Running
###################


DP1 <- as.data.frame(fread(TD1))
DP2 <- as.data.frame(fread(TD2))
DP3 <- as.data.frame(fread(TD3))
DP4 <- as.data.frame(fread(TD4))
ALL <- rbind(DP1, DP2, DP3, DP4)

# DP1
DP1 <- SAMPLE(DP1)
PARAMS1 <- PARAMS(DP1)
PARAMS1$BARCODE <- 'DP1'

# DP2
DP2 <- SAMPLE(DP2)
PARAMS2 <- PARAMS(DP2)
PARAMS2$BARCODE <- 'DP2'

# DP3
DP3 <- SAMPLE(DP3)
PARAMS3 <- PARAMS(DP3)
PARAMS3$BARCODE <- 'DP3'

# DP4
DP4 <- SAMPLE(DP4)
PARAMS4 <- PARAMS(DP4)
PARAMS4$BARCODE <- 'DP4'

# ALL
ALL <- SAMPLE(ALL)
PARAMETERS_ALL <- PARAMS(ALL)
PARAMETERS_ALL$BARCODE <- 'ALL'

# Parameters table
PARAMETERS <- rbind(PARAMS1, PARAMS2, PARAMS3, PARAMS4, PARAMETERS_ALL)

# Printing final table
write.table(x = PARAMETERS, file = OUTF, quote = F, sep = '\t', col.names = T, row.names = F)
