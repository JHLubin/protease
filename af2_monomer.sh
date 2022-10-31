#!/bin/bash
module purge
module load singularity/3.6.4
module use /projects/community/modulefiles
module load alphafold

input_fasta_file=$1
output_data_dir=$2

echo "--------------------------------" 
date
echo "--begin time---" 
bdate=$(date +%s)
echo $bdate 

echo "input fasta file"
echo $input_fasta_file
echo "Output data directory"
echo $output_data_dir

singularity run -B $ALPHAFOLD_DATA_PATH:/data -B .:/etc --pwd /app/alphafold --nv $CONTAINERDIR/alphafoldmm.sif \
    --fasta_paths=${input_fasta_file} \
    --output_dir=${output_data_dir} \
    --data_dir=/data \
    --uniref90_database_path=/data/uniref90/uniref90.fasta \
    --mgnify_database_path=/data/mgnify/mgy_clusters_2018_12.fa \
    --template_mmcif_dir=/data/pdb_mmcif/mmcif_files/ \
    --obsolete_pdbs_path=/data/pdb_mmcif/obsolete.dat \
    --model_preset=monomer \
    --db_preset=reduced_dbs \
    --small_bfd_database_path=/data/small_bfd/bfd-first_non_consensus_sequences.fasta \
    --pdb70_database_path=/data/pdb70/pdb70 \
    --max_template_date=2017-11-01 \
    --use_gpu_relax=True

echo "--end Time---" 
edate=$(date +%s)
echo $edate 
echo "--time spend---" 
echo $(( $bdate - $edate ))
echo "--------------------------------"