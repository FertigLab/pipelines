// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.data = ''
params.npatterns = 2
params.niterations = 100
params.sparse = 1


process preprocess {
    container 'satijalab/seurat:5.0.0'
    input:
        path data 
    output:
        path 'data.rds'
    """
    Rscript -e 'res <- Seurat::Read10X("$data");
                res <- Seurat::NormalizeData(res);
                saveRDS(res, file="data.rds")'
    """
}

process cogaps {
    container 'ghcr.io/fertiglab/cogaps:3.21.5'
    input:
      path 'matrix.mtx'
    output:
      path  'cogapsresult.rds'
    """
    Rscript -e 'library("CoGAPS");
        params <- new("CogapsParams");
        params <- setParam(params, "nPatterns", $params.npatterns);
        params <- setParam(params, "nIterations", $params.niterations);
        params <- setParam(params, "sparseOptimization", as.logical($params.sparse));
        res <- CoGAPS("matrix.mtx", params, outputFrequency = 10);
        saveRDS(res, file="cogapsresult.rds")'
    """
}



workflow {
  def input = Channel.fromPath(params.data)
  preprocess(input) | view
}