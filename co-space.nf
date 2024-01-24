//usage nextflow run --data ~/Documents/data/breastcancer/
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
params.data = ''
params.npatterns = 2
params.niterations = 100
params.sparse = 1

// add log transform if not yet transformed for cogaps, in Seurat
// propose adding Seutrat to cogaps inputs
// think about images and other human output for cogaps, markers
// think about parallel hotspots 

process preprocess {
    input:
        path data 
    output:
        path 'matrix.mtx'
    """
    gunzip -c $data/filtered_feature_bc_matrix/matrix.mtx.gz > matrix.mtx
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

process markers-preprocess {
  //1. call spacemarkers::load10Xcoords(data dir) -> 2 columns with x and y coords (needs cbind)
  //2. call spacemarkers::getSpatialFeatures(path cogaps results, methods) -> n cols with barcodes, 
  //3. cbind 1&2 -> spPatterns -> save 
}

process spacemarkers {
  //4. getSpatialParameters(spPatterns) -> hotspots (save)
}

process spacemarkers2 {
  //5. getInteractingGenes(hotposts, ref_patterns?) -> save as result
}

//visualizations

workflow {
  def input = Channel.fromPath(params.data)
  preprocess(input) | cogaps | view
}