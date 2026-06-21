process LABEL_TRANSFER {
    tag "$meta.id"
    label 'process_high_memory'
    
    // Container specified in conf/containers.config
    
    input:
    tuple val(meta), path(rds)
    path(reference_rds)

    output:
    // path "output", emit: results_dir
    tuple val(meta), path("**_labeled.rds") ,emit: rds
    tuple val(meta), path("**.csv*"), emit: data_csv
    tuple val(meta), path("**.png")  , emit: plots
    path "versions.yml", emit: versions
    script:
    """
    # Create the output directory
    mkdir label_transfer
     
    # Run the R script with explicit output directory
    label_transfer.r \\
        --inrds $rds\\
        --reference_rds $reference_rds \\
        --label_col $params.single_cell_label_col \\
        --outdir label_transfer \\
        --sample_name $meta.id\\
        --min_celltype_probability $params.min_celltype_probability \\
        --future_mem_limit $params.future_mem_limit \\
        --seed $params.seed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/R version //')
        seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
