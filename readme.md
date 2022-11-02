# Evolutionary selection of proteins with two folds

This program is intended to evaluate coevolutionary structural predictions of fold-switching proteins (also known as metamorphic proteins). Fold-switching proteins have the ability to transition between two sets of stable secondary and tertiary structure. This approach successfully revealed coevolution of amino acid pairs uniquely corresponding to both conformations of 58 fold-switching proteins from distinct families. 

This program utilizes both GREMLIN and MSA Transformer applied to multiple sequence alignments to generate structural predictions. Signal is enhanced by generating predictions from subfamily alignments using HHSUITE's QID filter.

For more information: [Porter lab](https://www.nlm.nih.gov/research/researchstaff/labs/porter/index.html)


    - For convenience, coevolution.yml is provided with all python dependencies to run fsc.py
    
    - Example run:
        python fsc.py --msa XXXX.msa --pdb1 XXXX.pdb  --pdb2 YYYY.pdb --Extra n/y

        Note: Extra flag defaults to n if set to y the scripts will generate an additional 6-panel image that compares
              output and runs the hypergeometric test

    - HMMER3.3.2 was used to generate the multiple sequence alignments to generate an alignment
      An example run is provided (update E-values as needed)
        jackhmmer -N 10 --cpu 16 -o XXXX.out -A XXXX.msa --tblout XXXX.tbl --incdomE 10E-20 --incE 10E-20 XXXX.fa uniref90.fasta
