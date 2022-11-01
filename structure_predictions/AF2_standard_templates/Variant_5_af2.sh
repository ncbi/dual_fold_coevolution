#!/bin/bash
module load alphafold2/2.1.2
run_singularity \
    --use_precomputed_msas \
    --model_preset=monomer \
    --fasta_paths=/data/porterll/AlphaFold/RfaH_variants_2.1.1/Variant_5.fasta \
    --max_template_date=2022-04-20 \
    --output_dir=/data/porterll/AlphaFold/RfaH_variants_2.1.1
