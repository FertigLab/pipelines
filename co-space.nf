// usage
// nextflow run co-space.nf --data ~/Documents/data/breastcancer -w wd -resume
// make sure wd has really low access requirements for docker to write there 

// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true, trouble is ep is /bin/bash
// to debug inside docker run e.g. 
// docker run -it \
//            --rm -v $(pwd):/spacemarkers \
//            --entrypoint /bin/bash ghcr.io/fertiglab/spacemarkers:x.y.z

nextflow.enable.dsl=2

// Script parameters
params.data = ''
params.npatterns = 5
params.niterations = 100
params.sparse = 1


process PREPROCESS {
  container 'docker.io/satijalab/seurat:5.0.0'

  input:
      path data 
  output:
      path 'dgCMatrix.rds'

  """
  Rscript -e 'res <- Seurat::Read10X("$data/raw_feature_bc_matrix/");
              res <- Seurat::NormalizeData(res);
              saveRDS(res, file="dgCMatrix.rds")';
  """
}

process COGAPS {
  container 'ghcr.io/fertiglab/cogaps:3.21.5'
  input:
    path 'dgCMatrix.rds'
  output:
    path  'cogapsResult.rds'

  """
  Rscript -e 'library("CoGAPS");
      sparse <- readRDS("dgCMatrix.rds");
      data <- as.matrix(sparse) #this converts to a dense matrix unfortunately
      cogapsResult <- CoGAPS(data = data, nPatterns = $params.npatterns, \
                              nIterations = $params.niterations, \
                              sparseOptimization = as.logical($params.sparse))
      saveRDS(cogapsResult, file = "cogapsResult.rds")'
  """
}

process SPACEMARKERS {
  container 'ghcr.io/fertiglab/spacemarkers:0.99.5'
  input:
    path data
    path 'cogapsResult.rds'
  output:
    path 'spPatterns.rds'
    path 'spaceMarkers.rds'

  """
  Rscript -e 'library("SpaceMarkers");
    dataMatrix <- load10XExpr("$data");
    coords <- load10XCoords("$data");
    features <- getSpatialFeatures("cogapsResult.rds");
    spPatterns <- cbind(coords, features);
    saveRDS(spPatterns, file = "spPatterns.rds");

    optParams <- getSpatialParameters(spPatterns);
    saveRDS(optParams, file = "optParams.rds");
    spaceMarkers <- getInteractingGenes(data = dataMatrix, \
                                        optParams = optParams, \
                                        spPatterns = spPatterns, \
                                        refPattern = "Pattern_1", \
                                        mode = "DE", analysis="enrichment");
              '
  """
}


workflow {
  def input = Channel.fromPath(params.data)
  PREPROCESS(input)
  COGAPS(PREPROCESS.out)
  SPACEMARKERS(input, COGAPS.out)
  view
}