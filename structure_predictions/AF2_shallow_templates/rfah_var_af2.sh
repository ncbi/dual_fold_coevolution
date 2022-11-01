#!/bin/bash
module load alphafold2/2.1.2
run_singularity \
    --use_precomputed_msas \
    --model_preset=monomer \
    --fasta_paths=/data/porterll/AlphaFold/Schafer_coevolution/shallow_msas2/rfah_var.fasta \
    --max_template_date=2022-04-20 \
    --output_dir=/data/porterll/AlphaFold/Schafer_coevolution/shallow_msas2
