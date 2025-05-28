
###### path to tools and directory ######

java="./bin/jdk-17.0.12/bin/java"
gatk="./bin/gatk-4.6.0.0/gatk-package-4.6.0.0-local.jar"

SAMPLE="Sample_name"
ref_dir="./ref/hg38"
work_dir="./DLBCL_CH"



###### reference genome alignment and processing ######

# load other required tools available on HPC

module load bwa/0.7.17-gcc-11.2.0
module load samtools/1.13-gcc-11.2.0

echo bwa mem ${SAMPLE} start at `date`

bwa mem ${ref_dir}/Homo_sapiens_assembly38.fasta FASTQ1 FASTQ2 | samtools addreplacerg --input-fmt SAM -r 'ID:${SAMPLE}' -r 'LB:${SAMPLE}' -r 'PL:ILLU' -r 'PU:ILLU01' -r 'SM:${SAMPLE}' -o ${work_dir}/bwa/${SAMPLE}.addRG.bam - 

samtools sort -t 20 ${work_dir}/bwa/${SAMPLE}.addRG.bam -o ${work_dir}/bwa/${SAMPLE}.sorted.bam
samtools index -@ 12 ${work_dir}/bwa/${SAMPLE}.sorted.bam

rm ${work_dir}/bwa/${SAMPLE}.addRG.bam

echo alignment ${SAMPLE} end at `date`



###### duplication marking ######

# duplication marking of the alignment results

module load python/3.11.7-gcc-12.3.0

echo gatk ${SAMPLE} dedupend start at `date`

${java} -jar ${gatk} MarkDuplicates \
    -ASO coordinate \
    -I ${work_dir}/bwa/${SAMPLE}.sorted.bam \
    -O ${work_dir}/dedup/${SAMPLE}.dedup.bam \
    -M ${work_dir}/dedup/${SAMPLE}.dedup.metrics.txt \
    --VALIDATION_STRINGENCY SILENT \
    --TMP_DIR ${work_dir}/tmp \

samtools index -@ 12 ${work_dir}/dedup/${SAMPLE}.dedup.bam

echo gatk ${SAMPLE} dedupend end at `date`



###### gatk BaseRecalibration ######

module load python/3.11.7-gcc-12.3.0

echo gatk ${SAMPLE} BaseRecalibrator start at `date`

${java} -jar ${gatk} BaseRecalibrator \
   -I ${work_dir}/dedup/${SAMPLE}.dedup.bam  \
   -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
   --known-sites ${ref_dir}/dbsnp_146.hg38.vcf.gz \
   --known-sites ${ref_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
   --known-sites ${ref_dir}/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
   -L ${ref_dir}/wgs_calling_regions.hg38.interval_list \
   -O ${work_dir}/cleanbam/${SAMPLE}.bqsr.table

echo gatk ${SAMPLE} BaseRecalibrator end at `date`


echo gatk ${SAMPLE} ApplyBQSR start at `date`

${java} -jar ${gatk} ApplyBQSR \
  -I ${work_dir}/dedup/${SAMPLE}.dedup.bam  \
  -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
  -L ${ref_dir}/wgs_calling_regions.hg38.interval_list \
  --bqsr-recal-file ${work_dir}/cleanbam/${SAMPLE}.bqsr.table \
  -O ${work_dir}/cleanbam/${SAMPLE}.bqsr.bam

rm ${work_dir}/dedup/${SAMPLE}.dedup.bam 

echo gatk ${SAMPLE} ApplyBQSR end at `date`



###### gatk DepthOfCoverage ######

module load python/3.11.7-gcc-12.3.0

echo DepthOfCoverage ${SAMPLE} start at `date`

${java} -jar ${gatk} DepthOfCoverage \
   -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
   -I ${work_dir}/cleanbam/${SAMPLE}.bqsr.bam \
   -L ${ref_dir}/SeqCap_EZ_Exome_v3_primary_targets.hg38.bed \
   -O ${work_dir}/coverage/coverage.${SAMPLE}.out.txt

echo DepthOfCoverage ${SAMPLE} done at `date`



###### mutation calling by Mutect2 ######

module load python/3.11.7-gcc-12.3.0

# mutation calling with mutect2

echo Mutect2 start at `date`

${java} -jar ${gatk} Mutect2 \
   -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
   -I ${work_dir}/cleanbam/${SAMPLE}.bqsr.bam \
   --germline-resource ${ref_dir}/af-only-gnomad.hg38.noCHIP.vcf.gz \
   --panel-of-normals ${ref_dir}/1000g_pon.hg38.noCHIP.vcf.gz \
   -O ${work_dir}/mutect2/mutect2.${SAMPLE}.raw.vcf.gz

echo Mutect2 done at `date`


# apply default filter by mutect2

echo FilterMutectCalls start at `date`

${java} -jar ${gatk} FilterMutectCalls \
   -R ${ref_dir}/Homo_sapiens_assembly38.fasta \
   -V ${work_dir}/mutect2/mutect2.${SAMPLE}.raw.vcf.gz \
   -O ${work_dir}/mutect2/mutect2.${SAMPLE}.filtered.vcf.gz

echo FilterMutectCalls finish at `date`



###### variant filtering ######

# additional filter by coverage, vaf, etc
bcftools annotate -x ^FORMAT/GT,^FORMAT/DP,^FORMAT/AD,^FORMAT/AF,^FORMAT/SB ${work_dir}/mutect2/mutect2.${SAMPLE}.filtered.vcf.gz | bcftools filter -e ' FORMAT/DP < 20 ||  FORMAT/AD<3 ||  FORMAT/SB<1 || FORMAT/AF<0.02 | FORMAT/AF>0.3' --SnpGap 10 --IndelGap 20 - -O z | bcftools view  -f 'PASS' - -O z > ${work_dir}/mutations/mutect2.${SAMPLE}.filtered.pass.vcf.gz && bcftools index -t ${work_dir}/mutations/mutect2.${SAMPLE}.filtered.pass.vcf.gz

# retain only CHIP mutations
bedtools intersect -a ${work_dir}/mutations/mutect2.${SAMPLE}.filtered.pass.vcf.gz -b CH_mutation_regions_MSK2023_NM2021.bed -header > ${work_dir}/mutations/mutect2.${SAMPLE}.filtered.pass.CHIP.vcf

## with VEP results
bcftools +split-vep -f '%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t[%${SAMPLE}=%AF:%AD:%SB\t]\n' -s worst -A tab ${work_dir}/mutations/mutect2.${SAMPLE}.filtered.pass.CHIP.vep.vcf | les > ${work_dir}/mutations/mutect2.${SAMPLE}.filtered.pass.CHIP.vep.txt
