#!/bin/bash
module load alphafold2/2.1.2
run_singularity \
    --use_precomputed_msas \
    --model_preset=monomer \
    --fasta_paths=/data/porterll/AlphaFold_new/RfaH_variants_2.1.1/delete_beta_coev_notemplates/Variant_5.fasta \
    --max_template_date=1960-04-20 \
    --output_dir=/data/porterll/AlphaFold_new/RfaH_variants_2.1.1/delete_beta_coev_notemplates
