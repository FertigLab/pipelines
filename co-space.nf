// usage
// sudo nextflow run co-space.nf --data ~/Documents/data/breastcancer/filtered_feature_bc_matrix
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.data = ''
params.npatterns = 2
params.niterations = 100
params.sparse = 1


process PREPROCESS {
    container 'satijalab/seurat:5.0.0'
    input:
        path data 
    output:
        path 'dgCMatrix.rds'
    """
    Rscript -e 'res <- Seurat::Read10X("$data");
                res <- Seurat::NormalizeData(res);
                saveRDS(res, file="dgCMatrix.rds")'
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


workflow {
  def input = Channel.fromPath(params.data)
  PREPROCESS(input)
  COGAPS(PREPROCESS.out)
  view
}