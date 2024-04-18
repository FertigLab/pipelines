// usage
// nextflow run main.nf --data ~/Documents/data/breastcancer -w wd -resume -profile docker
// make sure wd has really low access requirements for docker to write there 

// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true, trouble is ep is /bin/bash
// to debug inside docker run e.g. 
// docker run -it \
//            --rm -v $(pwd):/spacemarkers \
//            --entrypoint /bin/bash ghcr.io/fertiglab/spacemarkers:x.y.z

nextflow.enable.dsl=2

// Script parameters
params.data = ''
params.npatterns = 8
params.nsets = 7
params.niterations = 100
params.sparse = 1
params.seed = 42
params.distributed = '"genome-wide"'


process PREPROCESS {
  container 'docker.io/satijalab/seurat:5.0.0'

  input:
      path data 
  output:
      path 'dgCMatrix.rds'

  stub:
  """
  touch dgCMatrix.rds
  """

  script:
  """
  Rscript -e 'res <- Seurat::Read10X("$data/raw_feature_bc_matrix/");
              res <- Seurat::NormalizeData(res);
              saveRDS(res, file="dgCMatrix.rds")';
  """
}

process COGAPS {
  container 'ghcr.io/fertiglab/cogaps:3.21.5'
  cpus = params.nsets
  input:
    path 'dgCMatrix.rds'
  output:
    path  'cogapsResult.rds'

  stub:
  """
  touch cogapsResult.rds
  """

  script:
  """
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
      saveRDS(cogapsResult, file = "cogapsResult.rds")'
  """
}

process SPACEMARKERS {
  container 'ghcr.io/fertiglab/spacemarkers:0.99.8'
  input:
    path data
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

workflow COSPACE {
  def input = Channel.fromPath(params.data)
  PREPROCESS(input)
  COGAPS(PREPROCESS.out)
  SPACEMARKERS(input, COGAPS.out)
  view
}

workflow {
  COSPACE()
}