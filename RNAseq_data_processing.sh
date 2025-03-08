#!/bin/sh

#  Vaginal_gene_expression_analysis.sh
#  
#
#  Created by Rachel Petersen on 2/16/23.
#  
1. Trim reads

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=trim
#SBATCH --array=1-114


module load Trimmomatic/0.36

cd /nobackup/lea_lab/petersrm/Vaginal_gene_expression

data_path=$PWD/raw_data
trim_path=$PWD/trimmed

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
$data_path/$sampleID/*_L001_R1_001.fastq.gz $data_path/$sampleID/*_L001_R2_001.fastq.gz \
$trim_path/$sampleID.trim.L1.R1.fq.gz $trim_path/$sampleID.trim.L1.R1.unpaired.fq.gz \
$trim_path/$sampleID.trim.L1.R2.fq.gz $trim_path/$sampleID.trim.L1.R2.unpaired.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > \
$trim_path/$sampleID.L1.trimmomaticlog.txt 2>&1

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE -phred33 \
$data_path/$sampleID/*_L002_R1_001.fastq.gz $data_path/$sampleID/*_L002_R2_001.fastq.gz \
$trim_path/$sampleID.trim.L2.R1.fq.gz $trim_path/$sampleID.trim.L2.R1.unpaired.fq.gz \
$trim_path/$sampleID.trim.L2.R2.fq.gz $trim_path/$sampleID.trim.L2.R2.unpaired.fq.gz \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 > \
$trim_path/$sampleID.L2.trimmomaticlog.txt 2>&1

2. Index genome using STAR

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --job-name=star_index
#SBATCH --output=star_index.out

#module load GCC/6.4.0-2.28 Intel/2017.4.196 STAR/2.5.4b


cd /nobackup/lea_lab/petersrm/Vaginal_gene_expression

genome_path=$PWD/genomes

$PWD/STAR-2.7.10b/source/STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir $genome_path \
--genomeFastaFiles Panubis1.0.fna \
--sjdbGTFfile genomic.gtf \
--sjdbOverhang 99 \
--outFileNamePrefix Panubis1.0


3. Splice-aware alignment tool: STAR

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --job-name=star_align
#SBATCH --output=star_align_%A_%a.out
#SBATCH --array=1-114

module load GCC/6.4.0-2.28

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`

cd /nobackup/lea_lab/petersrm/Vaginal_gene_expression
trim_path=$PWD/trimmed
map_path=$PWD/mapped
genome_path=$PWD/genomes

$PWD/STAR-2.7.10b/source/STAR --genomeDir $genome_path \
--runThreadN 15 \
--readFilesIn $trim_path/$sampleID.trim.L1.R1.fq $trim_path/$sampleID.trim.L1.R2.fq \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 40712737493 \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--outFilterMatchNmin 0 \
--outFileNamePrefix $map_path/$sampleID.L1.

$PWD/STAR-2.7.10b/source/STAR --genomeDir $genome_path \
--runThreadN 15 \
--readFilesIn $trim_path/$sampleID.trim.L2.R1.fq $trim_path/$sampleID.trim.L2.R2.fq \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 40712737493 \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outFilterMismatchNoverReadLmax 0.04 \
--outFilterScoreMinOverLread 0.3 \
--outFilterMatchNminOverLread 0.3 \
--outFilterMatchNmin 0 \
--outFileNamePrefix $map_path/$sampleID.L2.

#filter for mapped reads

samtools view -b -F 4 $map_path/$sampleID.L1.Aligned.sortedByCoord.out.bam > $map_path/$sampleID.L1.mapped.bam
samtools view -b -F 4 $map_path/$sampleID.L2.Aligned.sortedByCoord.out.bam > $map_path/$sampleID.L2.mapped.bam


4. Build kraken database

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=load_libraries
#SBATCH --output=build_database.out

module purge
module load GCC/10.2.0  OpenMPI/4.0.5
module load Kraken2/2.1.1
module load BLAST+/2.11.0

cd /nobackup/lea_lab/petersrm/Vaginal_gene_expression/kraken_custom_db

kraken2-build --download-library bacteria --db /nobackup/lea_lab/petersrm/Vaginal_gene_expression/kraken_custom_db/microbe_anubis
kraken2-build --download-library archaea --db /nobackup/lea_lab/petersrm/Vaginal_gene_expression/kraken_custom_db/microbe_anubis
kraken2-build --download-library fungi --db microbe_anubis
kraken2-build --download-library protozoa --db microbe_anubis
kraken2-build --download-library viral --db microbe_anubis

#add anubis genome
kraken2-build --add-to-library /nobackup/lea_lab/petersrm/Vaginal_gene_expression/genomes/Panubis1.0.fna --db microbe_anubis

#build and clean db
kraken2-build --build --threads 15 --db microbe_anubis
kraken2-build --clean --db microbe_anubis
kraken2-inspect --db $PWD/microbe_anubis


5. Classify with Kraken

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --mem=190GB
#SBATCH --job-name=kraken-classify
#SBATCH --array=1-114

module load Kraken2/2.1.1

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`

cd /nobackup/lea_lab/petersrm/Vaginal_gene_expression
class_path=$PWD/kraken_classify/classified
unclass_path=$PWD/kraken_classify/unclassified
db_path=$PWD/kraken_custom_db
genome_path=$PWD/genomes
trim_path=$PWD/trimmed

kraken2 --threads 15 --use-names --classified-out $class_path/$sampleID.L1.cseqs#.fq \
--unclassified-out $unclass_path/$sampleID.L1.useqs#.fq \
--paired --report-zero-counts \
--db $db_path/microbe_anubis \
$trim_path/$sampleID.trim.L1.R1.fq \
$trim_path/$sampleID.trim.L1.R2.fq > classification_output

kraken2 --threads 15 --use-names --classified-out $class_path/$sampleID.L2.cseqs#.fq \
--unclassified-out $unclass_path/$sampleID.L2.useqs#.fq \
--paired --report-zero-counts \
--db $db_path/microbe_anubis \
$trim_path/$sampleID.trim.L2.R1.fq \
$trim_path/$sampleID.trim.L2.R2.fq > classification_output


7. Filter mapped reads for those classified as anubis or unclassified

grep_anubis.slurm

#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem=10GB
#SBATCH --job-name=grep_anubis
#SBATCH --array=1-114

sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p samples`

#Grep sequences classified as anubis
grep -A 3 'taxid|9555' $sampleID.L1.cseqs_1.fq | sed '/--/d' > $sampleID.L1.anubis_1.fq
grep -A 3 'taxid|9555' $sampleID.L1.cseqs_2.fq | sed '/--/d' > $sampleID.L1.anubis_2.fq
grep -A 3 'taxid|9555' $sampleID.L2.cseqs_1.fq | sed '/--/d' > $sampleID.L2.anubis_1.fq
grep -A 3 'taxid|9555' $sampleID.L2.cseqs_2.fq | sed '/--/d' > $sampleID.L2.anubis_2.fq

#Turn anubis fq files into read.ids.txt files
grep '@' $sampleID.L1.anubis_1.fq > $sampleID.anubis.read.ids.txt
grep '@' $sampleID.L1.anubis_2.fq >> $sampleID.anubis.read.ids.txt
grep '@' $sampleID.L2.anubis_1.fq >> $sampleID.anubis.read.ids.txt
grep '@' $sampleID.L2.anubis_2.fq >> $sampleID.anubis.read.ids.txt

#Turn unclassified fq files into read.ids.txt
grep '@' $sampleID.L1.useqs_1.fq > $sampleID.useqs.read.ids.txt
grep '@' $sampleID.L1.useqs_2.fq >> $sampleID.useqs.read.ids.txt
grep '@' $sampleID.L2.useqs_1.fq >> $sampleID.useqs.read.ids.txt
grep '@' $sampleID.L2.useqs_2.fq >> $sampleID.useqs.read.ids.txt

#concatonate anubis and useq read ids
cat $sampleID.anubis.read.ids.txt $sampleID.useqs.read.ids.txt > $sampleID.anubis.useqs.txt

#format read ids.
sed -i -e 's/ kraken:taxid|9555//g' $sampleID.anubis.useqs.txt
sed -i -e 's/@//' $sampleID.anubis.useqs.txt
sed -i -e 's/ .*$//' $sampleID.anubis.useqs.txt

#filter mapped reads for anubis_useq read names
module load picard/2.23.8

java -jar $PICARD_JAR FilterSamReads \
-INPUT $map_path/$sampleID.L1.mapped.bam \
-OUTPUT $map_path/$sampleID.L1.anubis.useqs.bam \
-READ_LIST_FILE $sampleID.anubis.useqs.txt
-FILTER includeReadList

java -jar $PICARD_JAR FilterSamReads \
-INPUT $map_path/$sampleID.L2.mapped.bam \
-OUTPUT $map_path/$sampleID.L2.anubis.useqs.bam \
-READ_LIST_FILE $sampleID.anubis.useqs.txt
-FILTER includeReadList

#merge L1 and L2
module load samtools/intel/1.12

samtools merge -@ 15 -f merged/$sampleID.anubis.useqs.merged.bam \
$sampleID.L1.anubis.useqs.bam \
$sampleID.L2.anubis.useqs.bam



8. Fix gene annotation file: Remove empty gene_id values

grep -v 'gene_id ""' genomic.gtf > FIXED_genomic.gtf


9. Run rsubread.R -- get read count matrix


#!/bin/sh
#
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=20GB

module load GCC/10.2.0  OpenMPI/4.0.5
module load R/4.0.5

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread", force=TRUE)

library(Rsubread)

bams = list.files(path = "/nobackup/lea_lab/petersrm/Vaginal_gene_expression/mapped/merged/", pattern = "*.anubis.useqs.merged.bam"
full.names=TRUE)
gtf.file = "/nobackup/lea_lab/petersrm/Vaginal_gene_expression/genomes/FIXED_genomic.gtf"
anubisfc = featureCounts(bams, annot.ext=gtf.file,
    isGTFAnnotationFile=TRUE,
    isPairedEnd=TRUE,
    nthreads=32,
    allowMultiOverlap=TRUE,
    tmpDir = "/nobackup/lea_lab/petersrm/Vaginal_gene_expression/Rsubread_tmp")

save.image("anubis_featurecounts.Rdata")
