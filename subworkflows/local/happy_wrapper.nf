//
// Run Happy preparation and benchmarking modules
//

//
// MODULE: Installed directly from nf-core/modules
//
include { HAPPY_HAPPY                 } from '../../modules/nf-core/modules/happy/happy/main'
include { HAPPY_PREPY                 } from '../../modules/nf-core/modules/happy/prepy/main'

workflow HAPPY_WRAP {
  take:

  inputs_channel

  main:
  
  //
  // CODE: Prepare input channels
  //
  ch_refgen_fasta = Channel.value([params.fasta, params.fasta_fai])
  //ch_refgen_fasta.view()
  ch_preppy_runs = inputs_channel.map{giab_id, sample_id, sample_vcf, sample_bam, ref_vcf, bed_file, bed_id -> tuple(['id': sample_id], sample_vcf, bed_file).flatten()}
  //ch_preppy_runs.view()

  //
  // MODULE: Prepare happy
  //
  HAPPY_PREPY (
    ch_preppy_runs,
    ch_refgen_fasta
  )
  
  //
  // MODULE: Run happy
  //

  //versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}