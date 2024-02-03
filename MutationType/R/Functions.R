
#' This function determines the mutation types by returning and constructing a vector of SNV from VCF file, reference genome
#' and given context length
#'
#' @usage determine_mutation_types(vcf, refgen, context_len)
#' @param vcf The VCF file
#' @param refgen The reference genome (e.g. BSgenome.Hsapiens.UCSC.hg19)
#' @param context_len Positive, odd integer specifying the context length
#' @return A vector with mutation types for each SNV in the format UP REF>ALT DOWN
#' @author Aigerim Nugmanova \cr Politecnico di Milano \cr E-Mail: <aigerim.nugmanova@mail.polimi.it>
#' @examples
#' library("VariantAnnotation")
#' library("BSgenome")
#' library("BSgenome.Hsapiens.UCSC.hg19")
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(vcffile, "hg19")
#' determine_mutation_types(vcf, Hsapiens, 3)
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import IRanges
#' @import BSgenome
#' @import Biostrings
#' @export
determine_mutation_types <- function(vcf, refgen, context_len) {

  if(is(refgen, 'character')){
    stop("This function cannot accept string format for reference genome (e.g hg19), please send it in the format of, for example, Hsapiens from BSgenome.Hsapiens.UCSC.hg19 library")
  }

  if(context_len <= 0 | context_len %% 2 == 0){
    stop("The context length should be odd positive integer")
  }

  seqnames <- seqnames(vcf)

  REF = ref(vcf)
  ALT = unlist(alt(vcf))

  start <- start(vcf) - (context_len-1)/2
  end <- end(vcf) + (context_len-1)/2

  if (all(startsWith(levels(seqnames), 'chr')) == FALSE){
    levels(seqnames) <- paste('chr', levels(seqnames), sep = '')
  }

  ran <- GRanges(seqnames = seqnames, IRanges(start = start, end = end), REF = as.character(REF), ALT = as.character(ALT))

  ran <- ran[nchar(ran$REF) == 1 & nchar(ran$ALT) == 1]

  seqn <- getSeq(refgen, ran)

  rev_seq <- reverseComplement(seqn)

  mut <- c()

  for(i in seq_len(length(seqn))) {
    if(ran$REF[i] == 'A' || ran$REF[i] == 'G') {
      mut[i] <- paste0(substr(as.character(rev_seq[i]), 1, (context_len-1)/2), "[", as.character(reverseComplement(DNAStringSet(ran$REF[i]))), ">", as.character(reverseComplement(DNAStringSet(ran$ALT[i]))), "]", substr(as.character(rev_seq[i]), ((context_len+1)/2)+1, context_len), sep = '')
    }
    else{
      mut[i] <- paste0(substr(as.character(seqn[i]), 1, (context_len-1)/2), "[", ran$REF[i], ">", ran$ALT[i], "]", substr(as.character(seqn[i]), ((context_len+1)/2)+1, context_len), sep = '')
    }
  }


  return(mut)

}

#' This function returns count table of the SNV mutation types and their frequency, also shows it in the graph
#'
#' @usage count_SNVmutations(vcf, refgen)
#' @param vcf The VCF file
#' @param refgen The reference genome (e.g. BSgenome.Hsapiens.UCSC.hg19)
#' @return A count table with mutation types for each SNV and their frequancy
#' @author Aigerim Nugmanova \cr Politecnico di Milano \cr E-Mail: <aigerim.nugmanova@mail.polimi.it>
#' @examples
#' library("VariantAnnotation")
#' library("BSgenome")
#' library("BSgenome.Hsapiens.UCSC.hg19")
#' vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#' vcf <- readVcf(vcffile, "hg19")
#' count_SNVmutations(vcf, Hsapiens)
#' @import VariantAnnotation
#' @import GenomicRanges
#' @import IRanges
#' @import BSgenome
#' @import Biostrings
#' @import ggplot2
#' @export
count_SNVmutations <- function(vcf, refgen) {

  mutations <- determine_mutation_types(vcf, refgen, 1)

  count <- as.data.frame(table(mutations), responseName = "Frequency")

  graph <- ggplot(count, aes(x=mutations,y=Frequency)) + geom_bar(stat='identity') + labs(title = "Frequency of each mutation type", x = "Mutation type", y = "Frequency")

  print(graph)

  return(count)
}
