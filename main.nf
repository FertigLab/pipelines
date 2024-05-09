// usage
// nextflow run main.nf --input ./samplesheet -w wd -resume -profile docker
// make sure wd has really low access requirements for docker to write there 

// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true, trouble is ep is /bin/bash
// to debug inside docker run e.g. 
// docker run -it \
//            --rm -v $(pwd):/spacemarkers \
//            --entrypoint /bin/bash ghcr.io/fertiglab/spacemarkers:x.y.z

nextflow.enable.dsl=2

// Script parameters
params.input = ''
params.npatterns = 2
params.nsets = 2
params.niterations = 10
params.sparse = 1
params.seed = 42
params.distributed = '"genome-wide"'


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

  stub:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"
  touch "${prefix}/cogapsResult.rds"
  """

  script:
  def args = task.ext.args ?: ''
  prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir "${prefix}"

  Rscript -e 'library("CoGAPS");
      sparse <- readRDS("dgCMatrix.rds");
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
  label 'process_single'
  container 'ghcr.io/fertiglab/spacemarkers:0.99.8'
  
  input:
    tuple val(sample), path(data)
    path 'cogapsResult.rds'
  output:
    path 'spPatterns.rds'
    path 'optParams.rds'
    path 'spaceMarkers.rds'

  stub:
  """
  touch spPatterns.rds
  touch optParams.rds
  touch spaceMarkers.rds
  """
  script:
  """
  Rscript -e 'library("SpaceMarkers");
    dataMatrix <- load10XExpr("$data");
    coords <- load10XCoords("$data");
    features <- getSpatialFeatures("cogapsResult.rds");
    spPatterns <- cbind(coords, features);
    saveRDS(spPatterns, file = "spPatterns.rds");

    #temp fix to remove barcodes with no spatial data
    barcodes <- intersect(rownames(spPatterns), colnames(dataMatrix))
    dataMatrix <- dataMatrix[,barcodes]
    spPatterns <- spPatterns[barcodes,]
    message("dim check:", dim(coords), ",", dim(dataMatrix))

    optParams <- getSpatialParameters(spPatterns);
    saveRDS(optParams, file = "optParams.rds");

    spaceMarkers <- getInteractingGenes(data = dataMatrix, \
                                        optParams = optParams, \
                                        spPatterns = spPatterns, \
                                        refPattern = "Pattern_1", \
                                        mode = "DE", analysis="enrichment");
    saveRDS(spaceMarkers, file = "spaceMarkers.rds");
              '
  """
}

// workflow {
//   def ss=Channel.fromPath(params.input)
//     | splitCsv(header:true, sep: ",")
//     | map { row-> tuple(sample=row.sample, data=file(row.data_directory)) }
//     PREPROCESS(ss)
//     COGAPS(PREPROCESS.out)
//     SPACEMARKERS(ss, COGAPS.out)
// }

// workflow {
// Channel.fromList([tuple([id: "sample"],
//         file(params.input))])
//         | map { tuple(it[0], it[1]) }
// | PREPROCESS

// }

workflow COSPACE {
  def samplesheet = Channel.fromPath(params.input)
    | splitCsv(header:true, sep: ",")
    | map { row-> tuple(meta=[id:row.sample], data=file(row.data_directory)) }
  
  PREPROCESS(samplesheet)

  ch_pre = PREPROCESS.out.dgCMatrix.map { tuple(it[0], it[1]) }
  COGAPS(ch_pre)

}

workflow {
  COSPACE()
}