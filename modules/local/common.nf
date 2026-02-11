// Merge TSVs and sum the specified column
// Currently support only headerless inputs and summing of the second column
process merge_and_publish_tsv {
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    label "wf_common"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta),
              path("inputs/input*.tsv")
        val(output_fname)
    output:
        tuple val(meta),
              path("${meta.alias}.${output_fname}")
    script:
    """
    find inputs -name "*.tsv" \
        -exec cat {} + \
        | csvtk -t summary -H -f 2:sum -g 1 \
        > "${meta.alias}.${output_fname}"
    """
}

process build_minimap_index {
    /*
    Build minimap index from reference genome
    */
    label "singlecell"
    cpus params.threads
    memory '16 GB'
    input:
        path "reference.fa"
    output:
        path "genome_index.mmi", emit: index
    script:
    """
    minimap2 -t ${task.cpus} -I 16G -d "genome_index.mmi" "reference.fa"
    """
}

process make_fasta_index {
    label "wf_common"
    memory "2 GB"
    cpus 1
    input:
        path "fasta.fa"
    output:
        path "fasta.fa.fai"
    script:
    """
    samtools faidx fasta.fa
    """
}

process call_paftools {
    label "singlecell"
    memory "2 GB"
    cpus 1
    input:
        path "ref_genes.gtf"
    output:
        path "ref_genes.bed", emit: ref_genes_bed
    script:
    """
    paftools.js gff2bed -j ref_genes.gtf > ref_genes.bed
    """
}

process cat_tags_by_chrom {
    // Merge per-chunk tags to create per-chromosome tags
    label "wf_common"
    cpus 4
    memory "12 GB"
    input:
        tuple val(meta),
              path('tags/*tags.tsv.zst')
    output:
        tuple val(meta),
              path("chr_tags/*"),
              emit: tags
        
    script:
        def tmpdir = "tmpdir"
        def buffer_size = "8192M"
        def sort_threads = 4
    """
    mkdir -p chr_tags
    mkdir -p ${tmpdir}
    
    first_file=\$(ls tags/*.tsv.zst | head -n1)
    [ -z "\$first_file" ] && { echo "No input files" >&2; exit 1; }
    
    # Get header from first line of first zstd file
    # head -n1 stops after first line and closes its stdin; zstdcat then receives SIGPIPE (exit 141).
    # So treat 141 as normal; any other non‑zero exit code is a real error.
    header=\$( (zstdcat "\$first_file" | head -n1) 2>&1 ) || {
        status=\$?
        if [ "\$status" -ne 141 ]; then
            echo "Error reading header: exit \$status" >&2
            exit "\$status"
        fi
    }
    [ -z "\$header" ] && { echo "Empty header" >&2; exit 1; }
    
    # Find chr column
    CHR_COL=\$(awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="chr"){print i; exit}}' <<< "\$header")
    CB_COL=\$(awk -F'\\t' '{for(i=1;i<=NF;i++) if(\$i=="CB"){print i; exit}}' <<< "\$header")
    [ -z "\$CHR_COL" ] && { echo "No chr column" >&2; exit 1; }
    [ -z "\$CB_COL" ] && { echo "No CB column" >&2; exit 1; }
    
    # Split by chromosome
    find tags -name '*.tsv.zst' -exec zstdcat {} + | \\
        awk -F'\\t' -v CHR_COL="\$CHR_COL" -v hdr="\$header" '
            BEGIN { split(hdr, h, FS); first_col = h[1] }
            \$1 == first_col { next }  # Skip headers
            {
                chr = \$CHR_COL
                out = "chr_tags/" chr ".tsv"
                if (!(seen[chr]++)) print hdr > out
                print > out
            }
        '
    
    # Sort by corrected cell barcode so downstream processess can read contiguous CB blocks.
    find -L chr_tags -name "*.tsv" -type f | while read -r file; do
        head -n 1 "\$file" | zstd -3 > "\${file}.zst"
        tail -n +2 "\$file" \
            | sort \
                --buffer-size=${buffer_size} \
                --temporary-directory=${tmpdir} \
                --parallel=${sort_threads} \
                --field-separator=\$'\\t' \
                --key=\${CB_COL},\${CB_COL} \
            | zstd -3 >> "\${file}.zst"
            rm \${file}
    done
    rm -rf ${tmpdir}
    """

}

process split_gtf_by_chroms {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        path(ref_gtf)
    output:
        path("*.gtf"), emit: chrom_gtf
    script:
    """
    if [[ "${ref_gtf}" == *.gz ]]; then
        cat_cmd="zcat" 
    else
        cat_cmd="cat"
        fi
    \${cat_cmd} ${ref_gtf} | awk '/^[^#]/ {print>\$1".gtf"}'
    """
}