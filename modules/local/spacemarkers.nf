process SPACEMARKERS {
  tag "$meta.id"
  label 'process_high_memory'
  container 'ghcr.io/fertiglab/spacemarkers:1.1.2.3'

  input:
    tuple val(meta), path(cogapsResult), path(data)
  output:
    tuple val(meta), path("${prefix}/spPatterns.rds"),   emit: spPatterns
    tuple val(meta), path("${prefix}/optParams.rds"),    emit: optParams
    tuple val(meta), path("${prefix}/spaceMarkers.rds"), emit: spaceMarkers
    path  "versions.yml",                                emit: versions

  stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}"
    touch "${prefix}/spPatterns.rds"
    touch "${prefix}/optParams.rds"
    touch "${prefix}/spaceMarkers.rds"
    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
  """
  script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    dir.create("${prefix}", showWarnings = FALSE)
    library("SpaceMarkers")
    
    #load spatial coords from tissue positions, deconvolved patterns, and expression
    coords <- load10XCoords("$data")
    features <- getSpatialFeatures("$cogapsResult")
    dataMatrix <- load10XExpr("$data")

    #add spatial coordinates to deconvolved data, only use barcodes present in data
    spPatterns <- merge(coords, features, by.x = "barcode", by.y = "row.names")
    spPatterns <- spPatterns[which(spPatterns[,"barcode"] %in% colnames(dataMatrix)),]
    saveRDS(spPatterns, file = "${prefix}/spPatterns.rds")

    #remove genes with low expression, only barcodes present in spatial data
    keepGenes <- which(apply(dataMatrix, 1, sum) > 10)
    keepBarcodes <- which(colnames(dataMatrix) %in% spPatterns[,"barcode"])
    dataMatrix <- dataMatrix[keepGenes, keepBarcodes]

    #compute optimal parameters for spatial patterns
    optParams <- getSpatialParameters(spPatterns);
    saveRDS(optParams, file = "${prefix}/optParams.rds")

    #find genes that are differentially expressed in spatial patterns
    spaceMarkers <- getPairwiseInteractingGenes(data = dataMatrix,
                                                  optParams = optParams,
                                                  spPatterns = spPatterns,
                                                  mode = "DE",
                                                  analysis="enrichment")

    saveRDS(spaceMarkers, file = "${prefix}/spaceMarkers.rds")

    # Get the versions of the packages
    spaceMarkersVersion <- packageVersion("SpaceMarkers")
    rVersion <- packageVersion("base")
    cat(sprintf('"%s":\n  SpaceMarkers: %s\n  R: %s\n', 
            "${task.process}", spaceMarkersVersion, rVersion), 
        file = "versions.yml")
    """
}

process SPACEMARKERS_MQC {
  tag "$meta.id"
  label 'process_low'
  container 'ghcr.io/fertiglab/spacemarkers:1.1.2.3'

  input:
    tuple val(meta), path(spaceMarkers)
  output:
    tuple val(meta), path("${prefix}/spacemarkers_mqc.json"), emit: spacemarkers_mqc
    path  "versions.yml",                                    emit: versions

  script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    dir.create("${prefix}", showWarnings = FALSE)

    #[['']] notation needed to allow nextflow var susbtitution

    sm <- readRDS("$spaceMarkers")
    smi <- sm[which(sapply(sm, function(x) length(x[['interacting_genes']]))>0)]

    #interacting patterns stats
    n_pairs_total <- length(sm)
    n_pairs_interact <- length(smi)

    #spacemarker metric
    max_spacemarker_metric <- max(sapply(smi, function(x) {
      max(x[['interacting_genes']][[1]][['SpaceMarkersMetric']])
    }))
    min_spacemarker_metric <- min(sapply(smi, function(x) {
      min(x[['interacting_genes']][[1]][['SpaceMarkersMetric']])
    }))

    #average number of genes in each pair
    avg_genes_in_pair <- mean(sapply(smi, function(x) {
      length(x[['interacting_genes']][[1]][['Gene']])
    }))

    #average percent overlap across interacting patterns
    avg_hotspot_area <- mean(sapply(smi, function(x) {
      sum(!is.na(x[['hotspots']]))/length(x[['hotspots']][,1])
    }))

    #report
    report_data <- list(
      "${prefix}" = list(
        n_pairs_total = n_pairs_total,
        n_pairs_interact = n_pairs_interact,
        max_spacemarker_metric = max_spacemarker_metric,
        min_spacemarker_metric = min_spacemarker_metric,
        avg_genes_in_pair = avg_genes_in_pair,
        avg_hotspot_area = avg_hotspot_area
      )
    )

    report <- list(
        id = "spacemarkers_mqc",
        section_name = "SpaceMarkers",
        description = "Stats for the spacemarkers run.",
        plot_type = "table",
        pconfig = list(
            id = "custom_data_table",
            title = "SpacemMarkers Stats"
            ),
        data = report_data
    )
    jsonlite::write_json(
              x=report, 
              path = "${prefix}/spacemarkers_mqc.json", 
              auto_unbox = TRUE, 
              pretty = TRUE)
    
    # Get the versions of the packages
    spaceMarkersVersion <- packageVersion("SpaceMarkers")
    rVersion <- packageVersion("base")
    cat(sprintf('"%s":\n  SpaceMarkers: %s\n  R: %s\n', 
            "${task.process}", spaceMarkersVersion, rVersion), 
        file = "versions.yml")
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir "${prefix}"
    touch "${prefix}/spacemarkers_mqc.json"
    cat <<-END_VERSIONS > versions.yml
      "${task.process}":
          SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
          R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}