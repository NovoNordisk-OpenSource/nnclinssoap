process SEURAT_XENIUM {
    tag "$meta.id"
    label 'process_medium'
    
    // Container specified in conf/containers.config
    
    input:
    tuple val(meta), path(xenium_path)


    output:
    // path "output", emit: results_dir
    tuple val(meta), path("**_processed.rds") ,emit: rds
    tuple val(meta), path("**.csv*") , emit: data_csv
    tuple val(meta), path("**.html") , emit: html
    tuple val(meta), path("**.png")  , emit: plots
    tuple val(meta), path("**.txt")  , emit: session_info
    path("**_average_qv_by_nCount_log10_bygene_transcript_type.png") , emit: qc_plots_to_report
    path "versions.yml", emit: versions
    script:
    def generic_name = params.generic_name != null ? "--genericName " + (params.generic_name ? "TRUE" : "FALSE") : ""
    """
    # Create the output directory
    mkdir output
       
    # Run the R script with explicit output directory
    Process_Individual_Xenium_Sample.r \\
        --indir ${xenium_path} \\
        --outdir output \\
        --sample_name $meta.id \\
        --min_reads_per_cell $params.min_reads_per_cell \\
        --max_reads_per_cell $params.max_reads_per_cell \\
        --max_control_prop_per_cell $params.max_control_prop_per_cell \\
        --n_highly_var_genes $params.n_highly_var_genes \\
        --seed $params.seed \\
        $generic_name \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | grep "R version" | sed 's/R version //')
        seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
    END_VERSIONS
    """
}
