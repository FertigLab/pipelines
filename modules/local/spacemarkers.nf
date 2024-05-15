process SPACEMARKERS {
  tag "$meta.id"
  label 'process_medium'
  container 'ghcr.io/fertiglab/spacemarkers:0.99.8'

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
      dataMatrix <- load10XExpr("$data");
      coords <- load10XCoords("$data");
      features <- getSpatialFeatures("$cogapsResult");
      spPatterns <- cbind(coords, features);
      saveRDS(spPatterns, file = "${prefix}/spPatterns.rds");

      #temp fix to remove barcodes with no spatial data
      barcodes <- intersect(rownames(spPatterns), colnames(dataMatrix))
      dataMatrix <- dataMatrix[,barcodes]
      spPatterns <- spPatterns[barcodes,]

      optParams <- getSpatialParameters(spPatterns);
      saveRDS(optParams, file = "${prefix}/optParams.rds");

      spaceMarkers <- getInteractingGenes(data = dataMatrix, \
                                          optParams = optParams, \
                                          spPatterns = spPatterns, \
                                          refPattern = "Pattern_1", \
                                          mode = "DE", analysis="enrichment");
      saveRDS(spaceMarkers, file = "${prefix}/spaceMarkers.rds");
                '
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SpaceMarkers: \$(Rscript -e 'print(packageVersion("SpaceMarkers"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
    END_VERSIONS
    """
}