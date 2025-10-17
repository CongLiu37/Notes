fastq="/bucket/BourguignonU/RAW_DATA/AzentaSQC_06062025/CCS/fastX/POR-A.hifireads.fastq.gz"
threads=50
host_genome.fna=""
#rmHost_minimap=""
nr="/bucket/BourguignonU/Cong/public_db/ncbi.blastdb_20250610/nr.dmdb.dmnd"
megan.db="/bucket/BourguignonU/Cong/public_db/megan.db/megan-nr-r2.mdb"
viralVerify.db="/bucket/BourguignonU/Cong/public_db/virusVerify/nbc_hmms.hmm"
out_dir="/flash/BourguignonU/Cong/gut_bacteria_termite/assembly"
out_prefix="PoroAdam.gut"
Step=1

fastq="/bucket/BourguignonU/Cong/gut_bacteria_termite/assembly/fastq.gz/m84168_251003_081605_s4_fastq.zip.fq.gz"
threads=32
host_genome.fna="/bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/Ncas/Ncas_maskedGenome.fna"
#rmHost_minimap=""
nr="/bucket/BourguignonU/Cong/public_db/ncbi.blastdb_20250610/nr.dmdb.dmnd"
megan.db="/bucket/BourguignonU/Cong/public_db/megan.db/megan-nr-r2.mdb"
viralVerify.db="/bucket/BourguignonU/Cong/public_db/virusVerify/nbc_hmms.hmm"
out_dir="/flash/BourguignonU/Cong/gut_bacteria_termite/assembly"
out_prefix="NeoSug.gut"
Step=1

# PoroAdam.gut	/bucket/BourguignonU/RAW_DATA/AzentaSQC_06062025/CCS/fastX/POR-A.hifireads.fastq.gz
# StoloVict.gut	/bucket/BourguignonU/RAW_DATA/AzentaSQC_06062025/CCS/fastX/STO-A.hifireads.fastq.gz

#find /bucket/BourguignonU/Cong/termite_pca/seqs.HOG/seqinr.kaks/ -maxdepth 1 -type f -name '*_seqinr.kaks.tsv' | xargs grep "Znev" | grep "Hsjo" > Hsjo_vs_Znev.kaks

cd /flash/BourguignonU/Cong/gut_bacteria_termite/assembly/find_termites
for i in `ls /bucket/BourguignonU/Cong/termite_genome_annotation/protein_3/*/*_maskedGenome.fna`
do
sp=`basename ${i} | sed 's/_maskedGenome.fna$//'`
echo ${sp}
seqkit replace -p '^' -r "${sp}." ${i} >>  47genome.fna
done
makeblastdb -dbtype nucl -in 47genome.fna
