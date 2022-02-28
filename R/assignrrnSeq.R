#'
#' Classifies sequences against reference training dataset.
#'
#' assignTaxonomy implements the RDP Naive Bayesian Classifier algorithm described in
#' Wang et al. Applied and Environmental Microbiology 2007, with kmer size 8 and 100 bootstrap
#' replicates. Properly formatted reference files for several popular taxonomic databases
#' are available \url{http://benjjneb.github.io/dada2/training.html}
#'
#' @param seqs (Required). A character vector of the sequences to be assigned, or an object
#' coercible by \code{\link{getUniques}}.
#'
#' @param rrndb (Required). The path to the reference rrndb fasta file, or an
#' R connection Can be compressed.
#'
#' @param minBoot (Optional). Default 50.
#' The minimum bootstrap confidence for assigning a taxonomic level.
#'
#' @param tryRC (Optional). Default FALSE.
#' If TRUE, the reverse-complement of each sequences will be used for classification if it is a better match to the reference
#' sequences than the forward sequence.
#'
#' @param outputBootstraps (Optional). Default FALSE.
#'  If TRUE, bootstrap values will be retained in an integer matrix. A named list containing the assigned taxonomies (named "taxa")
#'  and the bootstrap values (named "boot") will be returned. Minimum bootstrap confidence filtering still takes place,
#'  to see full taxonomy set minBoot=0
#'
#' @param multithread (Optional). Default is FALSE.
#'  If TRUE, multithreading is enabled and the number of available threads is automatically determined.
#'  If an integer is provided, the number of threads to use is set by passing the argument on to
#'  \code{\link{setThreadOptions}}.
#'
#' @param verbose (Optional). Default FALSE.
#'  If TRUE, print status to standard output.
#'
#' @return A character matrix of assigned taxonomies exceeding the minBoot level of
#'   bootstrapping confidence. Rows correspond to the provided sequences, columns to the
#'   taxonomic levels. NA indicates that the sequence was not consistently classified at
#'   that level at the minBoot threshhold.
#'
#'   If outputBootstraps is TRUE, a named list containing the assigned taxonomies (named "taxa")
#'   and the bootstrap values (named "boot") will be returned.
#'
#' @export
#'
#' @importFrom ShortRead readFasta
#' @importFrom ShortRead sread
#' @importFrom ShortRead id
#'

assignrrnSeq <- function(seqs, rrndb, minBoot=50, tryRC=FALSE, outputBootstraps=TRUE,
                           multithread=FALSE, verbose=FALSE) {
  MIN_REF_LEN <- 20 # Enforced minimum length of reference seqs. Must be bigger than the kmer-size used (8).
  MIN_TAX_LEN <- 50 # Minimum length of input sequences to get a taxonomic assignment
  # Get character vector of sequences
  seqs <- getSequences(seqs)
  if(min(nchar(seqs)) < MIN_TAX_LEN) {
    warning("Some sequences were shorter than ", MIN_TAX_LEN, " nts and will not receive a taxonomic classification.")
  }
  # Read in the reference fasta
  refsr <- readFasta(rrndb)
  lens <- width(sread(refsr))
  if(any(lens<MIN_REF_LEN)) {
    refsr <- refsr[lens>=MIN_REF_LEN]
    warning(paste0("Some reference sequences were too short (<", MIN_REF_LEN, "nts) and were excluded."))
  }
  refs <- unique(as.character(sread(refsr)))
  tax <- refs
  tax <- sapply(tax, function(x) gsub("^\\s+|\\s+$", "", x)) # Remove leading/trailing whitespace

  # Parse the taxonomies from the id string
  tax.depth <- sapply(strsplit(tax, ";"), length)
  td <- max(tax.depth)
  for(i in seq(length(tax))) {
    if(tax.depth[[i]] < td) {
      for(j in seq(td - tax.depth[[i]])) {
        tax[[i]] <- paste0(tax[[i]], "_DADA2_UNSPECIFIED;")
      }
    }
  }
  # Create the integer maps from reference to type ("genus") and for each tax level
  genus.unq <- unique(tax)
  ref.to.genus <- match(tax, genus.unq)
  tax.mat <- matrix(unlist(strsplit(genus.unq, ";")), ncol=td, byrow=TRUE)
  tax.df <- as.data.frame(tax.mat)
  for(i in seq(ncol(tax.df))) {
    tax.df[,i] <- factor(tax.df[,i])
    tax.df[,i] <- as.integer(tax.df[,i])
  }
  tax.mat.int <- as.matrix(tax.df)
  ### Assign
  # Parse multithreading argument
  if(is.logical(multithread)) {
    if(multithread==TRUE) { RcppParallel::setThreadOptions(numThreads = "auto") }
    else { RcppParallel::setThreadOptions(numThreads = 1) }
  } else if(is.numeric(multithread)) {
    RcppParallel::setThreadOptions(numThreads = multithread)
  } else {
    warning("Invalid multithread parameter. Running as a single thread.")
    RcppParallel::setThreadOptions(numThreads = 1)
  }
  # Run C assignemnt code
  assignment <- C_assign_taxonomy2(seqs, rc(seqs), refs, ref.to.genus, tax.mat.int, tryRC, verbose)

  # Parse results and return tax consistent with minBoot
  bestHit <- genus.unq[assignment$tax]
  boots <- assignment$boot
  taxes <- strsplit(bestHit, ";")
  taxes <- lapply(seq_along(taxes), function(i) taxes[[i]][boots[i,]>=minBoot])
  probs <- assignment$prob
  # const <- min(-700-min(probs),700-max(probs))
  const <- -700-apply(probs,1,min)
  adjusted <- probs+const # add a constant (to be canceled out) to log-posterior to avoid numerical issues.
  adjusted[adjusted > 700] = 700
  probs <- (exp(adjusted))/(rowSums(exp(adjusted)))
  # probs <- (exp(probs+const))/(rowSums(exp(probs+const))) # add a constant (to be canceled out) to log-posterior to avoid numerical issues.
  boot_prob <- assignment$boot_prob

  colnames(probs) <- genus.unq
  rownames(probs) <- seqs

  # Convert to character matrix
  tax.out <- matrix(NA_character_, nrow=length(seqs), ncol=td)
  for(i in seq(length(seqs))) {
    if(length(taxes[[i]]) > 0) {
      tax.out[i,1:length(taxes[[i]])] <- taxes[[i]]
    }
  }
  rownames(tax.out) <- seqs
  colnames(tax.out) <- "RefSeq"
  # tax.out[tax.out=="_DADA2_UNSPECIFIED"] <- NA_character_
  if(outputBootstraps){
    # Convert boots to integer matrix
    boots.out <- matrix(boots, nrow=length(seqs), ncol=td)
    rownames(boots.out) <- seqs
    colnames(boots.out) <- "RefSeq"

    boot_prob <- array(boot_prob,dim=c(nrow(probs),ncol(probs),100),dimnames = list(seq = seqs, genus = genus.unq, BOOT = 1:100))
    avg_boot_prob <- matrix(0,nrow=nrow(probs),ncol=ncol(probs))
    for (b in 1:100){
      # const <- min(-700-min(boot_prob[,,b]),700-max(boot_prob[,,b]))
      const <- -700-apply(boot_prob[,,b],1,min)
      adjusted <- boot_prob[,,b]+const # add a constant (to be canceled out) to log-posterior to avoid numerical issues.
      adjusted[adjusted > 700] = 700
      avg_boot_prob <- avg_boot_prob + exp(adjusted)/rowSums(exp(adjusted))
      # avg_boot_prob <- avg_boot_prob + exp(boot_prob[,,b]+const)/rowSums(exp(boot_prob[,,b]+const))
    }
    avg_boot_prob <- avg_boot_prob/100
    # boot_prob <- apply(boot_prob,"BOOT",function(x) apply(x,1,function(y) exp(y)/rowSums(exp(y))))
    # rownames(boot_prob) <- seqs
    # colnames(boot_prob) <- as.vector(sapply(1:100,function(x) paste(x,genus.unq,sep="|")))
    list(tax=tax.out, boot=boots.out, boot_tax=assignment$boot_tax, prob = probs, boot_prob = avg_boot_prob)
  } else {
    list(tax=tax.out, prob = probs)
  }
}
