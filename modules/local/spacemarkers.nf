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
    mkdir "${prefix}"
    Rscript -e 'library("SpaceMarkers");
      #load spatial coords from tissue positions, deconvolved patterns, and expression
      coords <- load10XCoords("$data")
      features <- getSpatialFeatures("$cogapsResult")
      dataMatrix <- load10XExpr("$data")

      #add spatial coordinates to deconvolved data, only use barcodes present in data
      spPatterns <- merge(coords, features, by.x = "barcode", by.y = "row.names")
      spPatterns <- spPatterns[which(spPatterns[,"barcode"] %in% colnames(dataMatrix)),]
      saveRDS(spPatterns, file = "${prefix}/spPatterns.rds");

      #remove genes with low expression, only barcodes present in spatial data
      keepGenes <- which(apply(dataMatrix, 1, sum) > 10)
      keepBarcodes <- which(colnames(dataMatrix) %in% spPatterns[,"barcode"])
      dataMatrix <- dataMatrix[keepGenes, keepBarcodes]

      #compute optimal parameters for spatial patterns
      optParams <- getSpatialParameters(spPatterns);
      saveRDS(optParams, file = "${prefix}/optParams.rds");

      #find genes that are differentially expressed in spatial patterns
      spaceMarkers <- getPairwiseInteractingGenes(data = dataMatrix, \
                                                  optParams = optParams, \
                                                  spPatterns = spPatterns, \
                                                  mode = "DE", \
                                                  analysis="enrichment");

      saveRDS(spaceMarkers, file = "${prefix}/spaceMarkers.rds");
                '
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}
