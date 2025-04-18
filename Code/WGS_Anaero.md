---
output:
  html_document: default
  pdf_document: default
---
# *Anaerostipes* whole genome assembly and non-R related code

### Author: DB

This is how all the programs that generated the data for R visualization were executed. All the work was done on the Clemson University HPCC: Palmetto with PBS initially and then later with SLURM. Some of the programs were utilised through Clemson University Genomics and Bioinformatics Facility (CUGBF) and some were downloaded on personal accounts.

## WGS ASSEMBLY:

### Organizing data :

Make necessary directories (Windows 'folder' == Mac/Linux 'directory'); I have it organized by `WGS_datedatareceived` directory separate from actual working directory called `data`

#### this directory also contains unzipped fastq files both forward and reverse

```         
mkdir data  ## makes directory data 
cd data  ## goes to data
mkdir data/wgs_daterecieved  ## makes a new directory inside data
cp ./wgs_data/wgs_daterecieved/*.fastq.gz ./data/wgs_daterecieved/  ##copies over all the fastq.gz files into new folder
cd ./data/wgs_daterecieved/  ## goes to working folder 
gunzip *.fastq.gz  #unzips .gz files and removes the original from the folder so keep a copy of the original in wgs_data 
```

## Running the whoe genome assembly pipeline

Input files: \*.fastq files for whole genome assembly\
Output files and folders: each step produces multiple files

Create a `sample_id` file in the directory where you will run the pipeline which lists the names of all the sample file names

```         
ls *R1_001.fastq |awk -F '_' '{print $1 "_" $2 "_" $3}' > ./sample_id
## for those filenames that have long names:
ls *contigs.db | sed 's|\(.*\)_.*|\1|' > sample2 
```

#### to check if the for loop will work

```         
    for i in `cat sample_id`; do echo ${i}; done
```

#### this will print the file names on screen

I use CUGBF softwares for easier run, as I dont have to install all the softwares eg trimgalore\
I have clubbed all the assembly steps together. You can separate them if you need only portions of the assembly.\
open nano/yourfavoriteeditor and paste the following:

```         
#!/bin/sh
#PBS -N pipeline_final 
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=1g,walltime=100:00:00 

source activate seekatz_wgs_assembly
module add samtools/1.13
module add prokka/1.14.5
module add bowtie2/2.5.4
module add spades/3.15.5

cd /scratch1/dbhatta/data/WGS_daterecieved/  


for i in `cat ./sample_id`; do
        ~/anaconda3/bin/trim_galore --fastqc --trim-n --paired ${i}_R1_001.fastq ${i}_R2_001.fastq --path_to_cutadapt ~/anaconda3/bin/cutadapt --output_dir ./${i}_trim/
        spades.py -o ./${i}_spades/ -1 ./${i}_trim/${i}_R1_001_val_1.fq -2 ./${i}_trim/${i}_R2_001_val_2.fq -k 55,77,127 --careful -t 8
        ~/anaconda3/bin/bowtie2-build ./${i}_spades/contigs.fasta ./${i}_spades/${i}_index #### I had bowtie2; you can use the one Rooksie has
        ~/anaconda3/bin/bowtie2 -x ./${i}_spades/${i}_index -1 ${i}_R1_001.fastq -2 ${i}_R2_001.fastq -S ./${i}_spades/sam_${i}.sam
        samtools view -b ./${i}_spades/sam_${i}.sam -o ./${i}_spades/bam_${i}.bam
        samtools sort ./${i}_spades/bam_${i}.bam -o ./${i}_spades/sort_${i}.bam
        samtools coverage ./${i}_spades/sort_${i}.bam | awk '{sum += $6}END{ print "Average = ", sum/NR}' > ./${i}_spades/${i}_average.txt
        prokka --outdir ./${i}_spades/prokka_${i} --prefix ${i} ./${i}_spades/contigs.fasta
        awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ./${i}_spades/prokka_${i}/${i}.ffn > ./${i}_spades/prokka_${i}/${i}_sl.ffn
        grep -A1 "16S ribosomal RNA" ./${i}_spades/prokka_${i}/${i}_sl.ffn > ${i}_16S.txt
done

qstat -xf $PBS_JOBID
```

## QUAST:

For quality metrics per genome assembly

```         
#!/bin/sh
#PBS -N quast_pipe 
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=1g,walltime=72:00:00 

module load quast/5.0.2

cd /scratch1/dbhatta/data/WGS_daterecieved/  

for i in `cat sample_id`; do
        quast.py --output-dir ./${i}_spades/quast_${i}/ ./${i}_spades/contigs.fasta
done

qstat -xf $PBS_JOBID
```

## MULTIQC:

To see all Quast output in a single txt or html file, use this otherwise you will be pulling all info from separate html files from Quast for 30-50 samples and it will take ages\
Using program from my personal account

```         
#!/bin/sh
#PBS -N multiqc_batchno
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=94gb:interconnect=1g,walltime=72:00:00

module load multiqc/1.27.1

cd /scratch1/dbhatta/wgs/WGS_batchno/

multiqc ./*_spades/quast_*/ --dirs --outdir multiqc_batchno

qstat -xf $PBS_JOBID
```

Once the multiqc runs and creates the directory multiqc_batchno, copy the multiqc_quast.txt from the multiqc_results directory into your desktop and open it in Excel. You will have all the values Quast gives for all of your strains in one sheet, make a note of that in WGS_id_spades.xlsx

### Additionally, concatenate all the `*_16S.txt` files into a single text file and run the 16S one by one in both NCBI and EZBiocloud to get the 16S ids. Note them down in WGS_id_spades.xlsx.

```         
for i in `cat sample_id3`; do sed -i '1{x;s/.*/fich=$(ps -p $PPID -o args=);fich=${fich##*\\} };echo ${fich%%.*}/e;G}' ${i}_16S.txt ; done # this adds the filename to the first line in the text
cat *_16S.txt > all_16S.txt
```

#### Now move `all_16S.txt` on your desktop to identify the 16S

### Do this for samplename_average.txt so you can get the depth of average in a single file as well and move it to your desktop

```         
cp ./*_spades/*_average.txt . #*_average.txt was in the spades folder, so I put it in the same folder as sample_id3
for i in `cat sample_id3`; do sed -i '1{x;s/.*/fich=$(ps -p $PPID -o args=);fich=${fich##*\\} };echo ${fich%%.*}/e;G}' ${i}_average.txt ; done
cat *_average.txt > all_average.txt
```

## DBCAN4:

the folder contains `dbcan_wgs` contains .fna files from prokka, so copy those over first along with the copy of `sample_id`

there is command for db_dir with run_dbcan it is to indicate where the db is. I installed the db from the run_dbcan site in the /scratch1/dbhatta for easier path

```         
cp ./WGS_20121/*_spades/prokka_*/*.fna ./dbcan_wgs
```

```         
#!/bin/sh
#SBATCH --job-name=dbcan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

module load anaconda3/2023.09-0
source activate /home/dbhatta/anaconda3/envs/run_dbcan4 ## personal dbcan4 envt

cd /scratch/dbhatta/dbcan

for i in `cat ./sample_anaero126`; do
    run_dbcan ./fna_anaero/${i}.fna prok -c cluster --dbcan_thread 36 --hmm_cpu 36 --dia_cpu 36 --db_dir /scratch/dbhatta/dbcan/db_dbcan4/ --out_dir dbcan_${i}
done
```

## ANVIO-8:

Basic steps for pangenome analysis; explanations of all commands are at: <https://merenlab.org/software/anvio/>

#### Create contigs database:

```         
#!/bin/sh
#SBATCH --job-name=anvi_init
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

source activate anvio-8 ##personal anvio8 envt

cd /scratch/dbhatta/anvio2/

#anvi-setup-kegg-data
#anvi-setup-ncbi-cogs --just-do-it --reset
for i in `cat ./sample_anaero126_3`; do
        anvi-script-reformat-fasta ./contigs/${i}_contigs.fasta -o ./contigs/${i}_contigs.fa -l 1000 --simplify-names --report-file ./contigs/${i}_report.txt
        anvi-gen-contigs-database -f ./contigs/${i}_contigs.fa -o ./contigs/${i}_contigs.db -n ${i}
        anvi-run-hmms -c ./contigs/${i}_contigs.db -T 36
        anvi-run-ncbi-cogs -c ./contigs/${i}_contigs.db -T 36
        anvi-run-kegg-kofams -c ./contigs/${i}_contigs.db -T 36
done
```

### dereplication of genomes

create extgen_anaero126.txt with 2 tab separated columns with 'name' and 'contigs_db_path' like this:

```         
name    contigs_db_path
CM1_51A_S205    ./cinn_contigs/CM1_51A_S205_contigs.db
CM1_51B_S206    ./cinn_contigs/CM1_51B_S206_contigs.db
CM1_51C_S207    ./cinn_contigs/CM1_51C_S207_contigs.db
CM1_52_S208 ./cinn_contigs/CM1_52_S208_contigs.db
CM1_53_S216 ./cinn_contigs/CM1_53_S216_contigs.db
```

Now dereplication was run

```         
#!/bin/sh
#SBATCH --job-name=anvi_derep
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

source activate anvio-8

cd /scratch/dbhatta/anvio2/

anvi-dereplicate-genomes -e ./extgen_anaero126.txt --output-dir ./dereplicate_a8_90 --program pyANI --similarity-threshold 0.9 -T 36
anvi-dereplicate-genomes -e ./extgen_anaero126.txt --output-dir ./dereplicate_a8_95 --program pyANI --similarity-threshold 0.95 -T 36
anvi-dereplicate-genomes -e ./extgen_anaero126.txt --output-dir ./dereplicate_a8_98 --program pyANI --similarity-threshold 0.98 -T 36
anvi-dereplicate-genomes -e ./extgen_anaero126.txt --output-dir ./dereplicate_a8_99 --program pyANI --similarity-threshold 0.99 -T 36
anvi-dereplicate-genomes -e ./extgen_anaero126.txt --output-dir ./dereplicate_a8_997 --program pyANI --similarity-threshold 0.997 -T 36
anvi-dereplicate-genomes -e ./extgen_anaero126.txt --output-dir ./dereplicate_a8_100 --program pyANI --similarity-threshold 1.0 -T 36
```

### Estimate metabolism

```         
#!/bin/sh
#SBATCH --job-name=anvi_met
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

source activate anvio-8

cd /scratch/dbhatta/anvio2/

anvi-estimate-metabolism -e ./extgen_anaero126.txt --matrix-format -O final_a8_anaero_mat
anvi-estimate-metabolism -e ./extgen_anaero126.txt -O final_a8_anaero
```

### pangenome creation

```         
#!/bin/sh
#SBATCH --job-name=anvi_pan
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

source activate anvio-8

cd /scratch/dbhatta/anvio2/

anvi-gen-genomes-storage -e ./extgen_anaero126.txt -o ./AN8-GENOMES.db
anvi-pan-genome -g ./AN8-GENOMES.db -n AN8ANPAN --output-dir ./PANGA8 --num-threads 36 --mcl-inflation 5 --use-ncbi-blast
anvi-compute-genome-similarity --external-genomes ./extgen_anaero126.txt --program pyANI --output-dir AN8ANI --num-threads 10 --pan-db ./PANGA8/AN8ANPAN-PAN.db
```

### pangenome visualization and functional enrichment across genomes

```         
## when in anvio-8 environment
anvi-display-pan -p ./cinn_contigs/CINN/CINNWGS-PAN.db -g ./cinn_contigs/CINN_GENOMES.db

anvi-compute-functional-enrichment-across-genomes -e extgen_anaero126.txt -o functional-enrichment.txt -G groups.txt --annotation-source FUNCTION_SOURCE
```

groups.txt contains 2 tab separated columns `isolate` and `characteristic`. The headers can be anything but it has to be 2 columns. FUNCTION_SOURCES as defined by anvio: KEGG_Module, KOfam, COG_Category etc.

## GTDB-TK:

Collected the fna file from prokka and put them in a file together

```         
#!/bin/sh
#SBATCH --job-name=gtdbtk
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=1000G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

source activate gtdbtk-2.4.0 ## personal gtdbtk envt

cd /scratch1/dbhatta/gtdbtk_db/

gtdbtk classify_wf --genome_dir ./fna_anaero/ --out_dir classify_wf_out_anaero --skip_ani_screen --cpus 36
```

## PATHOFACT:

Two files to be modded when running: config.yaml and pbs script that you create\
Put all .fna from prokka into a folder called fna_genus\
This will be where the results will be generated

### Example of config.yaml

```         
pathofact:
  sample: ["CM1_51A_S205", "CM1_51B_S206", "CM1_51C_S207", "CM1_52_S208", "CM1_53_S216"] # requires user input
  project: PathoFact_results_CM1 # requires user input
  datadir:  /scratch1/dbhatta/virulence_factors/pathofact/PathoFact/fna_cinn/ # requires user input
  workflow: "complete" #options: "complete", "AMR", "Tox", "Vir"
  size_fasta: 10000 #Adjustable to preference
  scripts: "scripts"
  signalp: "/scratch1/dbhatta/virulence_factors/pathofact/PathoFact/signalp-5.0b/bin/" # requires user input
  deepvirfinder: "submodules/DeepVirFinder/dvf.py"
  tox_hmm: "databases/toxins/combined_Toxin.hmm"
  tox_lib: "databases/library_HMM_Toxins.csv"
  tox_threshold: 40 #Bitscore threshold of the toxin prediction, adjustable by user to preference
  vir_hmm: "databases/virulence/Virulence_factor.hmm"
  vir_domains: "databases/models_and_domains"
  plasflow_threshold: 0.7
  plasflow_minlen: 1000
  runtime:
    short: "00:10:00"
    medium: "01:00:00"
    long: "02:00:00"
  mem:
    normal_mem_per_core_gb: "4G"
    big_mem_cores: 12 # change if required 
    big_mem_per_core_gb: "200G" # change if required
    
```

DeepVirFinder can't handle sequences longer than 2.1Mb. I tried removing sequences longer than 2.1Mb in my input, and then DeepVirFinder works. I also tried inputting only one sequence with 2.8Mb length, DeepVirFinder got stuck for a few hours and can't predict it. When I cut this same sequence into two smaller subsequences (2.1Mb and 0.7Mb), DeepVirFinder predicted it successfully in a few minutes.

### Exmaple of PBS script that will use the above config file

```         
#!/bin/sh
#PBS -N pathofact_cm1
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=251gb:interconnect=1g,walltime=100:00:00

source activate PathoFact ## personal pathofact envt

cd /scratch1/dbhatta/virulence_factors/pathofact/PathoFact/

snakemake -s Snakefile --configfile configcm1.yaml --use-conda --reason --cores 12 -p

qstat -xf $PBS_JOBID
```

## PHYLOGENETIC TREES:

### Core SNP tree

#### ROARY

Collect all gff files. Make sure prokka has run on the samples or get the gff files from ncbi (I did not do that for the dataset from gtdb, I put it through prokka)

```         
mkdir gff_cinn_all
cp *.gff ./gff_cinn_all
```

Paste the following in nano. The explanations of the actual commands are in: <https://sanger-pathogens.github.io/Roary/>

```         
#!/bin/sh
#PBS -N roary_cinn
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=12:mem=94gb:interconnect=1g,walltime=72:00:00

source activate roary1 ## personal roary envt

cd /scratch1/dbhatta/roary/gff_files/

roary -e --mafft -b blastp -n -v -p 8 *.gff

qstat -xf $PBS_JOBID
```

#### SNP-sites

```         
snp-sites -mvp -o snp_sites core_gene_alignment.aln 
```

#### RAXML

```         
#!/bin/sh
#PBS -N raxml_16S
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=12:mem=94gb:interconnect=1g,walltime=100:00:00

source activate raxml

cd /scratch1/dbhatta/raxml/raxml_cinn_final/

raxmlHPC-PTHREADS -s snp_sites_anaero.phylip -f a -k -m GTRGAMMA -p 12345 -x 1234 -N 500 -n snp_anaero


qstat -xf $PBS_JOBID
```

### EzAAI

This program uses fna files, converts to protein sequences, calculates AAI and clusters using UPGMA\
The following converts to protein sequences:

```         
#!/bin/sh
#SBATCH --job-name=ezaai
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=24:00:00

module load anaconda3/2023.09-0
source activate /home/dbhatta/.conda/envs/ezaai

cd /scratch/dbhatta/ezaai_anaero/

input_dir="./input_fasta"
output_dir="./ezaai_out"

for fna_file in $input_dir/*.fna; do
    base_name=$(basename "$fna_file" .fna)
    output_db="$output_dir/$base_name.db"
    ezaai extract -i "$fna_file" -o "$output_db" -l "$base_name"
done
```

The following creates the calculations and the tree:

```         
#!/bin/sh
#SBATCH --job-name=ezaai
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=24:00:00

module load anaconda3/2023.09-0
source activate /home/dbhatta/.conda/envs/ezaai

cd /scratch/dbhatta/ezaai_anaero/

ezaai calculate -i ./ezaai_out/ -j ./ezaai_out/ -o aai_anaero.tsv
ezaai cluster -i aai_anaero.tsv -o aai_anaero.nwk
```

## SCFA, GERMINATION, TOXIN GENES AND TREES:

Upload the genes fasta from uniprot\
Copy all faa files from genomes

### DIAMOND

The following is done in interactive:

```         
diamond makedb --in ./uniprotkb_ygfh_AND_taxonomy_id_2_2024_08_07.fasta -d ygfH_prop
diamond makedb --in ./uniprotkb_pduP_AND_taxonomy_id_2_2024_08_07.fasta -d pduP_prop
diamond makedb --in ./uniprotkb_pduC_AND_taxonomy_id_2_2024_08_07.fasta -d pduC_prop
diamond makedb --in ./uniprotkb_mmdA_AND_taxonomy_id_2_2024_08_07.fasta -d mmdA_prop
diamond makedb --in ./uniprotkb_methylmalonyl_CoA_mutase_AND_2024_08_07.fasta -d mut_prop
etc
```

This is the blast script:

```         
#!/bin/sh
#SBATCH --job-name=diamond
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=350G
#SBATCH --time=48:00:00
#SBATCH --constraint=interconnect_hdr

module load anaconda3/2023.09-0
source activate /home/dbhatta/anaconda3/envs/run_dbcan4

cd /scratch/dbhatta/scfa/

for g in `cat ./genes_acetate`; do
    for i in `cat ./sample_anaero126`; do
        diamond blastp -d ./acetate/${g}.dmnd -q ./faa_files/${i}.faa --threads 36 --id 50 --query-cover 80 --top 1 -o ./acetate/${i}_${g}_diamond.tsv
    done
done

for g in `cat ./genes_bukyrate`; do
    for i in `cat ./sample_anaero126`; do
        diamond blastp -d ./bukyrate/${g}.dmnd -q ./faa_files/${i}.faa --threads 36 --id 50 --query-cover 80 --top 1 -o ./bukyrate/${i}_${g}_diamond.tsv
    done
done

for g in `cat ./genes_propionate`; do
    for i in `cat ./sample_anaero126`; do
        diamond blastp -d ./propionate/${g}.dmnd -q ./faa_files/${i}.faa --threads 36 --id 50 --query-cover 80 --top 1 -o ./propionate/${i}_${g}_diamond.tsv
    done
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./gerAA_germ.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_gerAA_diamond.tsv
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./gerAB_germ.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_gerAB_diamond.tsv
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./spo0A_spor.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_spo0A_diamond.tsv
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./zona_tox.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_zona_diamond.tsv
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./cspB_germ.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_cspB_diamond.tsv
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./cspC_germ.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_cspC_diamond.tsv
done
for i in `cat ./sample_id`; do
    diamond blastp -d ./sleC_germ.dmnd -q ./faa_files/${i}.faa --threads 36 --id 70 --query-cover 80 --top 1 -o ./${i}_sleC_diamond.tsv
done
```

Move all the results in different folders based on the gene then R

### TREES

To create the trees, we need an aligned file of the target sequences. I created a list of target protein names from RStudio using the diamond.tsv file. I captured the target protein sequences in a single fasta file called `target_protein_fasta.fasta` for each gene

#### CLUSTALO

To align the sequences:

```         
clustalo -i ./target_protein_fasta.fasta --threads 4 -o out_gene.fa -v
```

#### RAXML

To create the trees:

```         
#!/bin/sh
#PBS -N raxml_bsh
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=250gb:interconnect=1g,walltime=100:00:00

source activate raxml

cd /scratch/dbhatta/theriot/clustal_raxml_tree/

raxmlHPC-PTHREADS -s out_gene.fa -f a -k -m PROTCATLG -p 12345 -x 1234 -N 500 -n raxml_bsh

qstat -xf $PBS_JOBID
```
