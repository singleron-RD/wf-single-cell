nextflow run /workspace/wf-single-cell \
    --fastq ./mt_10K.fastq.gz \
    --kit "3prime:GEXSCOPE-V2" \
    --expected_cells 100 \
    --ref_genome_dir '/workspace/test/genome/human_mt/' \
    --matrix_min_genes 1 \
    --matrix_max_mito 100 \
    --barcode_max_ed 3 \
    -profile standard \
    -resume \