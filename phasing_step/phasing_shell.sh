ls *.fasta > namelist.txt

for i in `ls *.fasta`;
do mv $i temp;
echo $i > namefile
Rscript extract_hybrid.R;
mv temp.fa $i;
rm -rf temp;
rm -rf namefile;
done;

gatk=`tail -n+1 phasing_settings | head -n1`
picard=`tail -n+2 phasing_settings | head -n1`
sequencing=`tail -n+5 phasing_settings | head -n1`
forward=`tail -n+6 phasing_settings | head -n1`

sed -i 's/\?/N/g' hybrid_ref.fa;
sed -i 's/-//g' hybrid_ref.fa;
bwa index -a is hybrid_ref.fa;
samtools faidx hybrid_ref.fa;
java -jar $picard CreateSequenceDictionary R=hybrid_ref.fa O=hybrid_ref.dict;

 if [ $sequencing == paired ]
 then
reverse=`tail -n+7 phasing_settings | head -n1`
bwa mem hybrid_ref.fa $forward $reverse > temp.sam;
fi

 if [ $sequencing == single ]
 then
bwa mem hybrid_ref.fa $forward > temp.sam;
fi

java -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=pixme SM=hybrid;
java -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
java -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
java -jar $gatk -T RealignerTargetCreator -R hybrid_ref.fa -I tempsortmarked.bam -o tempintervals.list;
java -jar $gatk -T IndelRealigner -R hybrid_ref.fa -I  tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
java -jar $gatk -T HaplotypeCaller -R hybrid_ref.fa -I temp_realigned_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 -o temp_raw_variants.vcf;
java -jar $gatk -T ReadBackedPhasing -R hybrid_ref.fa -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -o temp_phased_SNPs.vcf;
java -jar $gatk -T FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R hybrid_ref.fa -o temp_alt.fa;

Rscript onelining.R;

rm -rf hybrid_ref.*;
mv hybrid_ref2.fa hybrid_ref.fa;
rm -rf temp*;

sed -i 's/\?/N/g' hybrid_ref.fa;
sed -i 's/-//g' hybrid_ref.fa;
bwa index -a is hybrid_ref.fa;
samtools faidx hybrid_ref.fa;
java -jar $picard CreateSequenceDictionary R=hybrid_ref.fa O=hybrid_ref.dict;

 if [ $sequencing == paired ]
 then
reverse=`tail -n+7 phasing_settings | head -n1`
bwa mem hybrid_ref.fa $forward $reverse > temp.sam;
fi

 if [ $sequencing == single ]
 then
bwa mem hybrid_ref.fa $forward > temp.sam;
fi

java -jar $picard AddOrReplaceReadGroups I=temp.sam O=tempsort.sam SORT_ORDER=coordinate LB=rglib PL=illumina PU=pixme SM=hybrid;
java -jar $picard MarkDuplicates MAX_FILE_HANDLES=1000 I=tempsort.sam O=tempsortmarked.sam M=temp.metrics AS=TRUE;
java -jar $picard SamFormatConverter I=tempsortmarked.sam O=tempsortmarked.bam;
samtools index tempsortmarked.bam;
java -jar $gatk -T RealignerTargetCreator -R hybrid_ref.fa -I tempsortmarked.bam -o tempintervals.list;
java -jar $gatk -T IndelRealigner -R hybrid_ref.fa -I  tempsortmarked.bam -targetIntervals tempintervals.list -o temp_realigned_reads.bam;
java -jar $gatk -T HaplotypeCaller -R hybrid_ref.fa -I temp_realigned_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 30 -stand_call_conf 30 -o temp_raw_variants.vcf;
java -jar $gatk -T ReadBackedPhasing -R hybrid_ref.fa -I temp_realigned_reads.bam  --variant temp_raw_variants.vcf -o temp_phased_SNPs.vcf;
java -jar $gatk -T FastaAlternateReferenceMaker -V temp_phased_SNPs.vcf -R hybrid_ref.fa -o temp_alt.fa;

Rscript onelining.R;

mv hybrid_ref.fa safehybrid_ref.fa
rm -rf hybrid_ref.*;
mv hybrid_ref2.fa hybrid_ref.fa;
mv safehybrid_ref.fa hybrid_ref2.fa
rm -rf temp*;

for i in `ls *.fasta`;
name1="ONE_$i";
name2="TWO_$i";
cp $i $name1;
mv $i $name2;
done;



Rscript allelelifying.R;
