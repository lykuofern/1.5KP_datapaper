library("devtools")
library("dada2")

#learn error rate miseq, customized
path_el <- "/home/lykuo/lab_data/NGS_data/miseq/error_learn/SuperRed_35"
path_result <-"./"
fnFs <- sort(list.files(path_el, pattern=".r1.fq", full.names = TRUE))
fnRs <- sort(list.files(path_el, pattern=".r2.fq", full.names = TRUE))
errF <- learnErrors(fnFs, multithread=TRUE)
errR <- learnErrors(fnRs, multithread=TRUE)

numbers = c("01", "02", "03", "04", "05", "06", "07", "08", "09", 10:999)

#following the format like in "multiplex_table.txt" third element is the mimum ASV length, forth element is mimum r1 r2 ASVs overlapping length
rbcLN = c("fVGF", "rECL", "500", "4")
rbcLC = c("fNYG", "rVVG", "400", "20")
trnLF = c("L5675", "F4121", "180", "20")   

AP<-data.frame(rbcLN, rbcLC, trnLF)

multiplex <- read.table(
  "multiplex_cpDNAbarcode_clean.txt",
  sep="\t", header = TRUE)

for (a in 1:ncol(AP)){
  colnames(AP[a]) -> region
  AP[,a][1] -> Fp
  AP[,a][2] -> Rp
  as.numeric(AP[,a][3]) -> AP_minlength
  as.numeric(AP[,a][4]) -> minoverlap
  multiplex[!multiplex[,AP[,a][2]] %in% "",] -> amplicon
  miss = c()
  seqtable = c()
  dadamergfail = c()
  merg = c() #Kuo_modified
  nonmerg = c() #Kuo_modified
  path_demultiplex = paste0(path_result, region, "_demultiplex")
  path_trim <- paste0(path_demultiplex, "/trimmed")
  sort(list.files(path_trim, pattern=".r1.fq", full.names = FALSE))-> R1.names
  sort(list.files(path_trim, pattern=".r2.fq", full.names = FALSE))-> R2.names
  sort(list.files(path_trim, pattern=".r1.fq", full.names = TRUE))-> R1
  sort(list.files(path_trim, pattern=".r2.fq", full.names = TRUE))-> R2
  paste0(path_trim, "/filtered_", R1.names)-> filtFs
  paste0(path_trim, "/filtered_", R2.names)-> filtRs
  for(i in seq_along(R1)) {
    fastqPairedFilter(c(R1[i], R2[i]), c(filtFs[i], filtRs[i]),
                      verbose=TRUE, matchIDs = TRUE)
  }
  for (s in 1:nrow(amplicon)){
    s1 = paste0("filtered_trim_", region, "_", amplicon[s,Fp],"_", amplicon[s,Rp], "_r1.fq")
    r1 = paste0(path_trim, "/filtered_trim_", region, "_", amplicon[s,Fp],"_", amplicon[s,Rp], "_r1.fq")
    r2 = paste0(path_trim, "/filtered_trim_", region, "_", amplicon[s,Fp],"_", amplicon[s,Rp], "_r2.fq")
    r0 = paste0(amplicon[s,Fp],"_", amplicon[s,Rp])
    filename = paste(amplicon[s,1], amplicon[s,2], amplicon[s,3], amplicon[s,4], ".fas", sep = "_")  
    seqname = paste(amplicon[s,3], amplicon[s,4], amplicon[s,2], amplicon[s,1], sep = "_")       
    header = paste0(">",seqname)
    if (purrr::has_element(list.files(path= path_trim),s1)==TRUE){ 
      dadaFs <- dada(r1, err=errF, multithread=TRUE)
      dadaRs <- dada(r2, err=errR, multithread=TRUE)
      paste0(rep(header, length(dadaFs[["clustering"]][["abundance"]])), "_", numbers[1:length(dadaFs[["clustering"]][["abundance"]])], rep("_r1_", length(dadaFs[["clustering"]][["abundance"]])), sprintf(dadaFs[["clustering"]][["abundance"]]/sum(dadaFs[["clustering"]][["abundance"]]), fmt = '%#.3f'), rep("_abundance_", length(dadaFs[["clustering"]][["abundance"]])), dadaFs[["clustering"]][["abundance"]])->r1list
      cbind(r1list,dadaFs[["clustering"]][["sequence"]],dadaFs[["clustering"]][["abundance"]])-> r1fas
      r1fas[order(as.numeric(r1fas[,3]), decreasing = TRUE),]-> r1fas
      matrix(r1fas, ncol = 3)-> r1fas
      paste0(rep(header, length(dadaRs[["clustering"]][["abundance"]])), "_", numbers[1:length(dadaRs[["clustering"]][["abundance"]])], rep("_r2_", length(dadaRs[["clustering"]][["abundance"]])), sprintf(dadaRs[["clustering"]][["abundance"]]/sum(dadaRs[["clustering"]][["abundance"]]), fmt = '%#.3f'), rep("_abundance_", length(dadaRs[["clustering"]][["abundance"]])), dadaRs[["clustering"]][["abundance"]])->r2list
      cbind(r2list,dadaRs[["clustering"]][["sequence"]],dadaRs[["clustering"]][["abundance"]])-> r2fas
      r2fas[order(as.numeric(r2fas[,3]), decreasing = TRUE),]-> r2fas
      matrix(r2fas, ncol = 3)-> r2fas
      write.table(r1fas[,1:2], file = paste0(path_demultiplex, "/denoice/r1/",filename), append = FALSE, sep = "\n", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      write.table(r1fas[,1:2], file = paste0(path_demultiplex, "/denoice/r1/",r0), append = FALSE, sep = "\n", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      write.table(r2fas[,1:2], file = paste0(path_demultiplex, "/denoice/r2/",filename), append = FALSE, sep = "\n", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      write.table(r2fas[,1:2], file = paste0(path_demultiplex, "/denoice/r2/",r0), append = FALSE, sep = "\n", quote = FALSE,
                  row.names = FALSE, col.names = FALSE)
      try(mergers <- mergePairs(dadaFs, r1, dadaRs, r2, minOverlap = minoverlap, verbose=TRUE))
      if (nrow(mergers)>0 & max(nchar(mergers$sequence))>AP_minlength & max(mergers$abundance)>20){
        sum(mergers$abundance)->clustersum
        paste0(rep(header, length(mergers$abundance)), "_", numbers[1:length(mergers$abundance)], rep("_", length(mergers$abundance)), sprintf(mergers$abundance/clustersum, fmt = '%#.3f'), rep("_abundance_", length(mergers$abundance)), mergers$abundance)->merglist
        as.matrix(mergers)->mergers.table
        nchar(mergers.table[,1]) -> seqlen #Kuo_modified
        cbind(merglist,mergers.table,seqlen)-> fas0 #Kuo_modified
        matrix(fas0, ncol = 11)->fas0
        fas0[order(as.numeric(fas0[,3]), decreasing = TRUE),] -> fas0
        matrix(fas0, ncol = 11)->fas0 #Kuo_modified may write out this matrix
        rbind(merg, fas0)-> merg #Kuo_modified could be a table for dada2merge
        fas0[as.numeric(fas0[,11])>AP_minlength,] -> fas #Kuo_modified
        matrix(fas, ncol = 11)->fas
        if (nrow(fas)==0){
          cbind(r1, r2, header)->fail
          rbind(dadamergfail, fail)-> dadamergfail
        }
        write.table(fas[,1:2], file = paste0(path_demultiplex, "/denoice/", filename), append = FALSE, sep = "\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
        write.table(fas[1,1:2], file = paste0(path_demultiplex, "/denoice_best/",filename), append = FALSE, sep = "\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
      }else{
        cbind(r1, r2, header)->fail
        rbind(dadamergfail, fail)-> dadamergfail
      }
      
      try(nonmergers <- mergePairs(dadaFs, r1, dadaRs, r2, verbose=TRUE, justConcatenate=TRUE))#Kuo_modified in following lines else
      if (nrow(nonmergers)>0 ){
        sum(nonmergers$abundance)->clustersum
        paste0(rep(header, length(nonmergers$abundance)), "_", numbers[1:length(nonmergers$abundance)], rep("_", length(nonmergers$abundance)), sprintf(nonmergers$abundance/clustersum, fmt = '%#.3f'), rep("_abundance_", length(nonmergers$abundance)), nonmergers$abundance, rep("_10Ncat", length(nonmergers$abundance)))->nonmerglist
        as.matrix(nonmergers)->nonmergers.table
        cbind(nonmerglist,nonmergers.table)-> fascat #Kuo_modified
        matrix(fascat, ncol = 10)->fascat
        fascat[order(as.numeric(fascat[,3]), decreasing = TRUE),] -> fascat #Kuo_modified may write out this matrix
        matrix(fascat, ncol = 10)->fascat
        rbind(nonmerg, fascat)-> nonmerg #Kuo_modified could be a table for dada2nonmerge
        write.table(fascat[,1:2], file = paste0(path_demultiplex, "/denoice/nonmerged/",filename), append = FALSE, sep = "\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
        write.table(fascat[1,1:2], file = paste0(path_demultiplex, "/denoice_best/nonmerged/",filename), append = FALSE, sep = "\n", quote = FALSE,
                    row.names = FALSE, col.names = FALSE)
        }
      
      cbind(rep(r1, length(dadaFs[["clustering"]][["abundance"]])),r1list)-> seqtable01
      cbind(rep(r2, length(dadaRs[["clustering"]][["abundance"]])),r2list)-> seqtable02
      rbind(seqtable, seqtable01, seqtable02)-> seqtable
    }
    if (purrr::has_element(list.files(path= path_trim),s1)==FALSE){
      cbind(r1, header)->mis
      rbind(miss, mis)-> miss}
    }
  write.table(miss, file = paste0(path_demultiplex, "/denoice/missing_samples.txt"), append = FALSE, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  write.table(seqtable, file = paste0(path_demultiplex, "/denoice/sequence_table.txt"), append = FALSE, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  write.table(dadamergfail, file = paste0(path_demultiplex, "/denoice/dadamerge_fail.txt"), append = FALSE, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
}