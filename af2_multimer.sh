#!/bin/bash
export TF_FORCE_UNIFIED_MEMORY=1
export XLA_PYTHON_CLIENT_MEM_FRACTION=4.0
#export DB_DIR="/projects/datasets/alphafold/dbs2022.03"
export DB_DIR="/projects/community/alphafold/dbs"
export BIN_DIR="/projects/community/ai-fold/2021/bd387/envs/af2.2/bin/"

export PATH=/usr/lpp/mmfs/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin
module purge
module load gcc/11.2 openmpi
module use /projects/community/modulefiles
module load cuda/11.4.1
module load ai-fold
source activate af2.2

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

#CUDA_VISIBLE_DEVICES=0
python /projects/community/ai-fold/2021/bd387/envs/af2.2/alphafold/run_alphafold.py \
    --data_dir=$DB_DIR \
    --use_gpu_relax=FALSE \
    --uniref90_database_path=$DB_DIR/uniref90/uniref90.fasta \
    --uniclust30_database_path=$DB_DIR/uniclust30/uniclust30_2018_08/uniclust30_2018_08 \
    --mgnify_database_path=$DB_DIR/mgnify/mgy_clusters_2018_12.fa \
    --bfd_database_path=$DB_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt  \
    --template_mmcif_dir=$DB_DIR/pdb_mmcif/mmcif_files/ \
    --obsolete_pdbs_path=$DB_DIR/pdb_mmcif/obsolete.dat \
    --kalign_binary_path=$BIN_DIR/kalign \
    --jackhmmer_binary_path=$BIN_DIR/jackhmmer \
    --hhblits_binary_path=$BIN_DIR/hhblits \
    --hhsearch_binary_path=$BIN_DIR/hhsearch \
    --pdb_seqres_database_path=$DB_DIR/pdb_seqres/pdb_seqres.txt \
    --uniprot_database_path=$DB_DIR/uniprot/uniprot.fasta \
    --model_preset=multimer \
    --num_multimer_predictions_per_model=2 \
    --db_preset=full_dbs \
    --fasta_paths=$input_fasta_file \
    --output_dir=$output_data_dir \
    --max_template_date=2020-12-01 

echo "--end Time---" 
edate=$(date +%s)
echo $edate 
echo "--time spend---" 
echo $(( $bdate - $edate ))
echo "--------------------------------"