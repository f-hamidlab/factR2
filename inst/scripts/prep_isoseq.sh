wget https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/P14_M5_
wget https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/P28_M1_
wget https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/P14_M5_
wget https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/P14_M5_
wget https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/P14_M5_
wget https://data.nemoarchive.org/biccn/grant/rf1_tilgner/tilgner/transcriptome/scell/pacbio/mouse/raw/P14_M5_

wd=/media/cdn-bc/RAID/Datasets/Tilgner_Biorxiv_2023_BrainLR/

mkdir originals
while IFS=$'\t' read -r NAME CLASS LINK CLASSNUM;do
    wget $LINK -P FASTQ_sc
    tar -xvf "FASTQ_sc/$NAME.fastq.tar" -C FASTQ_sc
    rm FASTQ_sc/$NAME.fastq.tar
done < sc_fastq_curated.txt

for file in *.tar;do
    folder=${file%.fastq.tar}
    tar -xvf $file 
    mv $folder/${folder}.fastq.gz . 
    #rm $file
done


for file in *.fastq.gz;do
    echo $file
done


conda activate isoquant
files=(*.gz)
isoquant.py --reference /media/cdn-bc/RAID/Genomes/GRCm39_mm39/Gencode_M28/GRCm39.genome.fa.gz \
 --fastq ${files[@]} \
  --data_type pacbio_ccs \
   -o /media/cdn-bc/RAID/Datasets/Tilgner_Biorxiv_2023_BrainLR/outputs/isoquant_out

cd /media/cdn-bc/RAID/Datasets/Tilgner_Biorxiv_2023_BrainLR/outputs/isoquant_out/OUT
gzip -c OUT.transcript_models.gtf > ~/Github_repo/factR2/inst/extdata/pb_custom.gtf.gz

echo feature_id	P14_M5_1	P14_M5_2	P14_M5_3	P14_M6_1	P14_M6_2	P14_M6_3	P28_M1_1	P28_M1_2	P28_M1_3	P28_M2_1	P28_M2_2	P28_M2_3 > ~/Github_repo/factR2/inst/extdata/pb_expression.tsv

tail -n +2 OUT.transcript_model_grouped_counts.tsv >> ~/Github_repo/factR2/inst/extdata/pb_expression.tsv
gzip ~/Github_repo/factR2/inst/extdata/pb_expression.tsv