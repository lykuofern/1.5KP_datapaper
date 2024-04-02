#starting miseq trim read files 'trim_pool-1_R1.fastq.gz' 'trim_pool-1_R2.fastq.gz' the trim PE read files of BioSample accession SAMN40650549
#using (universal) primer sequences to demultiplex the amplicon
#rbcL C terminal
#"TAGGTCTGTCTGCYAARAATTATGG" and "GTTCCCCYTCTAGTTTRCCTACTAC" are VARIABLEs (sequences of rbcLC primers)
#"rbcLC" the name of the amplicon, which is a VARIABLE too
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --discard-untrimmed --pair-filter=both --minimum-length 70 --pair-adapters -g TAGGTCTGTCTGCYAARAATTATGG -G GTTCCCCYTCTAGTTTRCCTACTAC --action=none -o rbcLC.1_R1 -p rbcLC.1_R2 trim_pool-1_R1.fastq.gz trim_pool-1_R2.fastq.gz -j 30
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --discard-untrimmed --pair-filter=both --minimum-length 70 --pair-adapters -g GTTCCCCYTCTAGTTTRCCTACTAC -G TAGGTCTGTCTGCYAARAATTATGG --action=none -o rbcLC.2_R1 -p rbcLC.2_R2 trim_pool-1_R1.fastq.gz trim_pool-1_R2.fastq.gz -j 30
cat rbcLC.1_R1 rbcLC.2_R2 > rbcLC_amplicon_r1.fq
cat rbcLC.1_R2 rbcLC.2_R1 > rbcLC_amplicon_r2.fq
rm rbcLC.*_R*
mkdir rbcLC_demultiplex
cd rbcLC_demultiplex
#"../barcodes_rbcLC_start_0.fasta" and "../barcodes_rbcLC_start2_0.fasta" are VARIABLEs (containing indexed sequences of demultiplex)
#{name1} and {name2} are header names of fasta
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0 --no-indels --pair-filter=both --discard-untrimmed -g file:../barcodes_rbcLC_start_0.fasta -G file:../barcodes_rbcLC_start2_0.fasta --action=none -o rbcLC_{name1}_{name2}_r1.fq -p rbcLC_{name1}_{name2}_r2.fq ../rbcLC_amplicon_r1.fq ../rbcLC_amplicon_r2.fq -j 30


for File in *r1.fq
	do
	/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --minimum-length 150 --no-indels -g TAGGTCTGTCTGCYAARAATTATGG  -o trim_${File} ${File} -j 30
	done

for file in *r2.fq
	do
	/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --minimum-length 150 --no-indels -g GTTCCCCYTCTAGTTTRCCTACTAC  -o trim_${file} ${file} -j 30
	done


mkdir trimmed
mv trim_* ./trimmed/
mkdir denoice
mkdir denoice_best
mkdir denoice/r1
mkdir denoice/r2
mkdir denoice/nonmerged
mkdir denoice_best/nonmerged
cd ..

# using (universal) primer sequences to demultiplex
#rbcL N terminal, explanation of the other variables are similar as above
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --discard-untrimmed --pair-filter=both --minimum-length 70 --pair-adapters -g GAGACTAAAGCAGGTGTTGGATTCA -G TCAAGTCCACCRCGAAGRCATTC --action=none -o rbcLN.1_R1 -p rbcLN.1_R2 trim_pool-1_R1.fastq.gz trim_pool-1_R2.fastq.gz -j 30
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --discard-untrimmed --pair-filter=both --minimum-length 70 --pair-adapters -g TCAAGTCCACCRCGAAGRCATTC -G GAGACTAAAGCAGGTGTTGGATTCA --action=none -o rbcLN.2_R1 -p rbcLN.2_R2 trim_pool-1_R1.fastq.gz trim_pool-1_R2.fastq.gz -j 30
cat rbcLN.1_R1 rbcLN.2_R2 > rbcLN_amplicon_r1.fq
cat rbcLN.1_R2 rbcLN.2_R1 > rbcLN_amplicon_r2.fq
rm rbcLN.*_R*
mkdir rbcLN_demultiplex
cd rbcLN_demultiplex
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0 --no-indels --pair-filter=both --discard-untrimmed -g file:../barcodes_rbcL_start_0.fasta -G file:../barcodes_rbcLN_start2_0.fasta --action=none -o rbcLN_{name1}_{name2}_r1.fq -p rbcLN_{name1}_{name2}_r2.fq ../rbcLN_amplicon_r1.fq ../rbcLN_amplicon_r2.fq -j 30


for File in *r1.fq
	do
	/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --minimum-length 200 -g GAGACTAAAGCAGGTGTTGGATTCA  -o trim_${File} ${File} -j 30
	done

for file in *r2.fq
	do
	/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --minimum-length 200 -g TCAAGTCCACCRCGAAGRCATTC -o trim_${file} ${file} -j 30
	done
	

mkdir trimmed
mv trim_* ./trimmed/
mkdir denoice
mkdir denoice_best
mkdir denoice/r1
mkdir denoice/r2
mkdir denoice/nonmerged
mkdir denoice_best/nonmerged
cd ..


# using (universal) primer sequences to demultiplex
#trnLF, explanation of the other variables are similar as above
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --discard-untrimmed --pair-filter=both --minimum-length 10 --pair-adapters -g TGAGGGTTCGANTCCCTCTA -G GGATTTTCAGTCCYCTGCTCT --action=none -o trnLF.1_R1 -p trnLF.1_R2 trim_pool-1_R1.fastq.gz trim_pool-1_R2.fastq.gz -j 30
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --discard-untrimmed --pair-filter=both --minimum-length 10 --pair-adapters -g GGATTTTCAGTCCYCTGCTCT -G TGAGGGTTCGANTCCCTCTA --action=none -o trnLF.2_R1 -p trnLF.2_R2 trim_pool-1_R1.fastq.gz trim_pool-1_R2.fastq.gz -j 30
cat trnLF.1_R1 trnLF.2_R2 > trnLF_amplicon_r1.fq
cat trnLF.1_R2 trnLF.2_R1 > trnLF_amplicon_r2.fq
rm trnLF.*_R*
mkdir trnLF_demultiplex
cd trnLF_demultiplex
/home/lykuo/cutadapt-venv/bin/cutadapt -e 0 --no-indels --pair-filter=both --discard-untrimmed -g file:../barcodes_trnL_3exonSTART_0.fasta -G file:../barcodes_trnF_0.fasta --action=none -o trnLF_{name1}_{name2}_r1.fq -p trnLF_{name1}_{name2}_r2.fq ../trnLF_amplicon_r1.fq ../trnLF_amplicon_r2.fq -j 30

for File in *r1.fq
	do
	/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --minimum-length 50 -g TGAGGGTTCGANTCCCTCTA -o trim_${File} ${File} -j 30
	done

for file in *r2.fq
	do
	/home/lykuo/cutadapt-venv/bin/cutadapt -e 0.125 --no-indels --minimum-length 50 -g GGATTTTCAGTCCYCTGCTCT -o trim_${file} ${file} -j 30
	done

mkdir trimmed
mv trim_* ./trimmed/
mkdir denoice
mkdir denoice_best
mkdir denoice/r1
mkdir denoice/r2
mkdir denoice/nonmerged
mkdir denoice_best/nonmerged
cd ..

