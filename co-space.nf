// usage
// nextflow run co-space.nf --data ~/Documents/data/breastcancer/filtered_feature_bc_matrix -w wd
// make sure wd directory has really low access requirements for docker to be able to populate it
// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true in case entrypoint of container is not /bin/bash
// to debug inside docker e.g. docker run -it --rm -v $(pwd):/spacemarkers ghcr.io/fertiglab/spacemarkers:0.0.0

nextflow.enable.dsl=2

// Script parameters
params.data = ''
params.npatterns = 2
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
    path  'cogapsresult.rds'

  """
  Rscript -e 'library("CoGAPS");
      sparse <- readRDS("dgCMatrix.rds");
      data <- as.matrix(sparse) #this converts to a dense matrix unfortunately
      cogapsresult <- CoGAPS(data = data, nPatterns = $params.npatterns, \
                              nIterations = $params.niterations, \
                              sparseOptimization = as.logical($params.sparse))
      saveRDS(cogapsresult, file = "cogapsresult.rds")'
  """
}

process SPACEMARKERS {
  container 'ghcr.io/fertiglab/spacemarkers:0.0.0'
  input:
    path data
    path 'cogapsresult.rds'
  output:
    path 'sppatterns.rds'
    path 'hotspots.rds'

  """
  Rscript -e 'library("SpaceMarkers");
    coords <- load10XCoords("$data");
    features <- getSpatialFeatures("cogapsresult.rds");
    patterns <- cbind(coords, features);
    saveRDS(patterns, file = "sppatterns.rds");

    hotspots <- getSpatialParameters(patterns);
    saveRDS(hotspots, file = "hotspots.rds");'
  """
}


workflow {
  def input = Channel.fromPath(params.data)
  PREPROCESS(input)
  COGAPS(PREPROCESS.out)
  SPACEMARKERS(input, COGAPS.out)
  view
}