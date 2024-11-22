nextflow.enable.dsl=2

include { PREPROCESS } from './modules/local/preprocess.nf'
include { COGAPS } from './modules/local/cogaps/nextflow/'
include { SPACEMARKERS; 
          SPACEMARKERS_MQC;
          SPACEMARKERS_IMSCORES } from './modules/local/spacemarkers/nextflow/'

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

  ch_spacemarkers_mqc = SPACEMARKERS.out.spaceMarkers.map { tuple(it[0], it[1]) }
  SPACEMARKERS_MQC(ch_spacemarkers_mqc)

  ch_spacemarkers_imscores = SPACEMARKERS.out.spaceMarkers.map { tuple(it[0], it[1]) }
  SPACEMARKERS_IMSCORES(ch_spacemarkers_imscores)

  emit:
    dgCMatrix       = PREPROCESS.out.dgCMatrix
    cogapsResult    = COGAPS.out.cogapsResult
    spPatterns      = SPACEMARKERS.out.spPatterns
    optParams       = SPACEMARKERS.out.optParams
    spaceMarkers    = SPACEMARKERS.out.spaceMarkers
    versions        = SPACEMARKERS.out.versions
    spacemarkers_mqc = SPACEMARKERS_MQC.out.spacemarkers_mqc
    spacemarkers_imscores = SPACEMARKERS_IMSCORES.out.spacemarkers_imscores
}

workflow {
    COSPACE()
}