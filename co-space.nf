// usage
// nextflow run co-space.nf --data ~/Documents/data/breastcancer/filtered_feature_bc_matrix -w wd
// make sure wd directory has really low access requirements for docker to be able to populate it
// export NXF_CONTAINER_ENTRYPOINT_OVERRIDE=true in case entrypoint of container is not /bin/bash

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
    Rscript -e 'res <- Seurat::Read10X("$data");
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

process SPM_PATTERNS {
    container 'ghcr.io/fertiglab/spacemarkers:0.0.0'
    input:
      path data
      path 'cogapsresult.rds'
    output:
      path 'sppatterns.rds'
    """
    Rscript -e 'library("SpaceMarkers");
      coords <- load10Xcoords("$data");
      features <- getSpatialFeatures("cogapsresult.rds");
      patterns <- cbind(coord, features);
      saveRDS(patterns, file = "sppatterns.rds")';
    """
}


workflow {
  def input = Channel.fromPath(params.data)
  PREPROCESS(input)
  COGAPS(PREPROCESS.out)
  //SPM_PATTERNS(input, COGAPS.out)
  view
}