#!/bin/bash
: '
  Copyright (c) 2015, Fonleap Ltd
  All rights reserved.

  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:

  1. Redistributions of source code must retain the above copyright notice, this
     list of conditions and the following disclaimer.

  2. Redistributions in binary form must reproduce the above copyright notice, this
     list of conditions and the following disclaimer in the documentation and/or
     other materials provided with the distribution.

  3. Neither the name of the copyright holder nor the names of its contributors may
     be used to endorse or promote products derived from this software without
     specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
  INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
  NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
  OF SUCH DAMAGE.
'

TOOLSPATH=/enter/your/path/to/tools

PERL=$TOOLSPATH/perl/perl
BWA=$TOOLSPATH/bwa-0.7.12/bwa
SAMTOOLS=$TOOLSPATH/samtools/samtools
BCFTOOLS=$TOOLSPATH/bcftools/bcftools
VCFUTILS=$TOOLSPATH/vcfutils/vcfutils.pl
GATKPATH=$TOOLSPATH/gatk/GenomeAnalysisTK.jar

AMREF=$TOOLSPATH/1000genomes/hs37d5.seq
REF=$TOOLSPATH/1000genomes/hs37d5.fa
REFNAME=$TOOLSPATH/1000genomes/hs37d5

INDELMILLS=$TOOLSPATH/variants/Mills_and_1000G_gold_standard.indels.b37.vcf
DBSNP=$TOOLSPATH/variants/dbsnp_138.b37.vcf

if [ "$#" -lt 5 ]
then
    echo "A script to execute Samtools WGS/WES mapping to variant calls (version 1.0)"
    echo "Usage:"
    echo $0 /path/to/pair1.fastq /path/to/pair2.fastq /path/to/unpaired.fastq /output/dir/ iddesc [testid] [seqid]
    echo "    If you have no unpaired data, enter a path to an empty file"
    echo "    iddesc is used to identify reads (e.g. SRR622461)"
    echo "    seqid is used to identify sequence sample (e.g. NA12878)"
    echo "    testid can be used to identify particular instances of pipeline execution"
    echo "    use /output/dir/ to point to the path where the intermediate and output files will be stored"
    exit
fi

function checkrc {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "error with $1" >&2
        exit 1
    fi
    return $status
}

PAIR1=$1
PAIR2=$2
PAIR0=$3
IDDESC=$5
if [ "$#" -ge 6 ]
then
    TESTID=$6
    OUTPUT="${4}/${IDDESC}.${TESTID}.${FILEPOSTFIX}"
else
    OUTPUT="${4}/${IDDESC}.${FILEPOSTFIX}"
fi
SEQID=NA12878
if [ "$#" -ge 7 ]
then
    SEQID=$7
fi


if [ ! -e "${REF}.fai" ]
then
    checkrc $SAMTOOLS faidx $REF
    date
    echo "Reference indexing completed!"
fi

FILEPOSTFIX=samtools.workflow

BAMFILE=${OUTPUT}.bam

date
echo "Running BWA mem alignment..."
$BWA mem -t8 -R"@RG\tID:${IDDESC}\tSM:NA12878\tPL:ILLUMINA" $REF ${PAIR0} | $SAMTOOLS view -b - -o ${BAMFILE}.up.bam
$BWA mem -t8 -R"@RG\tID:${IDDESC}\tSM:NA12878\tPL:ILLUMINA" $REF ${PAIR1} ${PAIR2} | $SAMTOOLS view -b - -o ${BAMFILE}.p.bam

date
echo "BWA mem alignment completed. Running cleanup (samtools fixmate)..."
checkrc $SAMTOOLS fixmate -O bam ${BAMFILE}.up.bam ${BAMFILE}.up.fm.bam
checkrc $SAMTOOLS fixmate -O bam ${BAMFILE}.p.bam ${BAMFILE}.p.fm.bam

date
echo "Merging..."
checkrc $SAMTOOLS merge $BAMFILE ${BAMFILE}.up.fm.bam ${BAMFILE}.p.fm.bam

#tmp file cleanup
checkrc rm ${BAMFILE}.up.bam ${BAMFILE}.p.bam ${BAMFILE}.up.fm.bam ${BAMFILE}.p.fm.bam

date
echo "Sorting..."
$SAMTOOLS sort -@8 -m 1073741824 -O bam -o ${BAMFILE}.sorted -T ${BAMFILE}.sort.tmp $BAMFILE
checkrc mv ${BAMFILE}.sorted $BAMFILE
checkrc $SAMTOOLS index $BAMFILE

date
echo "Mapping, conversion, cleaning, sorting and indexing completed!"

date
echo "Starting improvement phase..."
# realign indels
checkrc $JAVA -jar $GATKPATH -T RealignerTargetCreator -R $REF -I $BAMFILE -o $BAMFILE.intervals --known $INDELMILLS -nt 8
checkrc $JAVA -jar $GATKPATH -T IndelRealigner -R $REF -I $BAMFILE -targetIntervals $BAMFILE.intervals -known $INDELMILLS -o $BAMFILE.ra.bam
checkrc mv $BAMFILE.ra.bam $BAMFILE
checkrc mv $BAMFILE.ra.bai $BAMFILE.bai
checkrc rm $BAMFILE.intervals
checkrc $SAMTOOLS index $BAMFILE

# BQSR
checkrc $JAVA -jar $GATKPATH -T BaseRecalibrator -R $REF -knownSites $DBSNP -I $BAMFILE -o $BAMFILE.recal.table
checkrc $JAVA -jar $GATKPATH -T PrintReads -R $REF -I $BAMFILE -BQSR $BAMFILE.recal.table -o $BAMFILE.recal.bam
checkrc mv $BAMFILE.recal.bam $BAMFILE
checkrc mv $BAMFILE.recal.bai $BAMFILE.bai
checkrc rm $BAMFILE.recal.table
checkrc $SAMTOOLS index $BAMFILE
# realign indels again!
checkrc $JAVA -jar $GATKPATH -T RealignerTargetCreator -R $REF -I $BAMFILE -o $BAMFILE.intervals --known $INDELMILLS -nt 8
checkrc $JAVA -jar $GATKPATH -T IndelRealigner -R $REF -I $BAMFILE -targetIntervals $BAMFILE.intervals -known $INDELMILLS -o $BAMFILE.ra.bam
checkrc mv $BAMFILE.ra.bam $BAMFILE
checkrc mv $BAMFILE.ra.bai $BAMFILE.bai
checkrc rm $BAMFILE.intervals
checkrc $SAMTOOLS index $BAMFILE

date
echo "Improvement phase completed!"

date
# run mpile in parallel (map)
echo "Running samtools mpileup..."
$SAMTOOLS view -H $BAMFILE | grep '\@SQ' | sed 's/^.*SN://g' | cut -f 1 | xargs -I {} -n 1 -P 8 sh -c "${SAMTOOLS} mpileup -d 100000 -g -f ${REF} -r {} ${BAMFILE} | ${BCFTOOLS} call -vmO v | ${VCFUTILS} varFilter -2 0 > ${BAMFILE}.flt2_0.tmp.{}.vcf"
# merge results (reduce)
ITEMLIST=$($SAMTOOLS view -H $BAMFILE | grep '\@SQ' | sed 's/^.*SN://g' | awk '{print $1}')
ITEM0=$($SAMTOOLS view -H $BAMFILE | grep '\@SQ' | sed 's/^.*SN://g' | awk '{print $1}' | head -n 1)
grep "#" ${BAMFILE}.flt2_0.tmp.$ITEM0.vcf > ${BAMFILE}.flt2_0.vcf
for item in $ITEMLIST;
do
  cat ${BAMFILE}.flt2_0.tmp.$item.vcf | grep -v "#" >> ${BAMFILE}.flt2_0.vcf
  rm ${BAMFILE}.flt2_0.tmp.$item.vcf
done

date
echo "Samtools mpileup completed"
echo "Completed calling variants"
