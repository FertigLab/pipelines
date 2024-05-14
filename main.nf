// usage
// nextflow run main.nf --input ./[samplesheet] -w wd -resume -profile docker
// make sure wd has really low access requirements for docker to write there 

// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true, to ensure ENTRYPOINT is 
// the /bin/bash command and not anything else that the author has specified
// to debug inside docker run e.g. 
// docker run -it \
//            --rm -v $(pwd):/spacemarkers \
//            --entrypoint /bin/bash ghcr.io/fertiglab/spacemarkers:x.y.z

nextflow.enable.dsl=2

process PREPROCESS {
  tag "$meta.id"
  label 'process_single'
  container 'docker.io/satijalab/seurat:5.0.0'

  input:
      tuple val(meta), path(data) 
  output:
      tuple val(meta), path("${prefix}/dgCMatrix.rds"), emit: dgCMatrix
      path "versions.yml"                             , emit: versions

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"

  """
  mkdir "${prefix}"
  touch "${prefix}/dgCMatrix.rds"
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"

  Rscript -e 'res <- Seurat::Read10X("$data/raw_feature_bc_matrix/");
              res <- Seurat::NormalizeData(res);
              saveRDS(res, file="${prefix}/dgCMatrix.rds")';

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seurat: \$(Rscript -e 'print(packageVersion("Seurat"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}

process COGAPS {
  tag "$meta.id"
  label 'process_long'
  container 'ghcr.io/fertiglab/cogaps:3.21.5'
  cpus = params.nsets

  input:
    tuple val(meta), path(dgCMatrix)
  output:
    tuple val(meta), path("${prefix}/cogapsResult.rds"), emit: cogapsResult
    path  "versions.yml",                                emit: versions

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  touch "${prefix}/cogapsResult.rds"
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  Rscript -e 'library("CoGAPS");
      sparse <- readRDS("$dgCMatrix");
      data <- as.matrix(sparse);
      params <- CogapsParams(seed=42,
                             nIterations = $params.niterations,
                             nPatterns = $params.npatterns,
                             sparseOptimization = as.logical($params.sparse),
                             distributed=$params.distributed);
      params <- setDistributedParams(params, nSets = 7);
      cogapsResult <- CoGAPS(data = data, params = params)
      saveRDS(cogapsResult, file = "${prefix}/cogapsResult.rds")'

  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CoGAPS: \$(Rscript -e 'print(packageVersion("CoGAPS"))' | awk '{print \$2}')
        R: \$(Rscript -e 'print(packageVersion("base"))' | awk '{print \$2}')
  END_VERSIONS
  """
}

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

workflow COSPACE {

  samplesheet = Channel.fromPath(params.input)
    .splitCsv(header:true, sep: ",")
    .map { row-> tuple(meta=[id:row.sample], data=file(row.data_directory)) }

  PREPROCESS(samplesheet)

  ch_gaps = PREPROCESS.out.dgCMatrix.map { tuple(it[0], it[1]) }
  COGAPS(ch_gaps)

  ch_spacemarkers = COGAPS.out.cogapsResult.map { tuple(it[0], it[1]) }
  .join(samplesheet)

  SPACEMARKERS(ch_spacemarkers)

  emit:
    dgCMatrix       = PREPROCESS.out.dgCMatrix
    cogapsResult    = COGAPS.out.cogapsResult
    spPatterns      = SPACEMARKERS.out.spPatterns
    optParams       = SPACEMARKERS.out.optParams
    spaceMarkers    = SPACEMARKERS.out.spaceMarkers
    versions        = SPACEMARKERS.out.versions

}

workflow {
    COSPACE()
}