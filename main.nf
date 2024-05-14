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

include { PREPROCESS } from './modules/local/preprocess'
include { SPACEMARKERS } from './modules/local/spacemarkers'
include { COGAPS } from './modules/local/cogaps'

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