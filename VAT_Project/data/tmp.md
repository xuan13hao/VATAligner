./kraken2 --paired --threads 16 --report /home/h392x566/data/16S_UNC/kraken_results/kraken_taxonomy.txt --output /home/h392x566/data/16S_UNC/kraken_results/kraken_output.txt  --db .     /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/Gr5m2_2_F_filt.fastq.gz  /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/Gr5m2_2_R_filt.fastq.gz



for i in *R*;
do
  j=$(echo $i | sed s/_R_filt.fastq.gz//)_F_filt.fastq.gz;
  k=$(echo $i | sed s/_L002_R\._001.fastq.gz//);

  ./home/h392x566/data/16S_UNC/kraken2/kraken2 --paired --threads 16 --report /home/h392x566/data/16S_UNC/kraken_results/$i_taxonomy.txt --output /home/h392x566/data/16S_UNC/kraken_results/$i_output.txt  --db   /home/h392x566/data/16S_UNC/kraken2   /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/$i  /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/$j

done

  ./home/h392x566/data/16S_UNC/kraken2/kraken2 --paired --threads 16 --report /home/h392x566/data/16S_UNC/kraken_results/$i_taxonomy.txt --output /home/h392x566/data/16S_UNC/kraken_results/$i_output.txt  --db   /home/h392x566/data/16S_UNC/kraken2   /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/$i  /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/$j

  /home/h392x566/data/16S_UNC/kraken2/kraken2 --paired --threads 16 --report /home/h392x566/data/16S_UNC/kraken_results/$i_taxonomy.txt --output /home/h392x566/data/16S_UNC/kraken_results/$i_output.txt  --db   /home/h392x566/data/16S_UNC/kraken2   /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/$j  /home/h392x566/data/16S_UNC/UNC_Data/alldata1/filtered/$i