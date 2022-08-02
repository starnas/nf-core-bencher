//
// Run Happy preparation and benchmarking modules
//

//
// MODULE: Installed directly from nf-core/modules
//
include { HAPPY_PREPY                 } from '../../modules/nf-core/modules/happy/prepy/main'
include { HAPPY_HAPPY                 } from '../../modules/nf-core/modules/happy/happy/main'

workflow HAPPY_WRAP {
  take:

  ch_formatted_inputs

  main:
  
  //
  // CODE: Prepare preppy input channel
  //
  ch_refgen_fasta = Channel.value([params.fasta, params.fasta_fai])
  //ch_refgen_fasta.view()
  ch_preppy_runs = ch_formatted_inputs.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed -> tuple(meta, sample_vcf, final_region_bed).flatten()}                
  //ch_preppy_runs.view()

  //
  // MODULE: Prepare happy
  //
  HAPPY_PREPY (
    ch_preppy_runs,
    ch_refgen_fasta
  )
  
  //
  // CODE: Prepare happy input channel
  //
  ch_happy_runs = ch_formatted_inputs.combine(HAPPY_PREPY.out.preprocessed_vcf, by: 0)
  ch_happy_runs = ch_happy_runs.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed, sample_preppy_vcf -> tuple(meta, ref_vcf, sample_preppy_vcf, final_region_bed).flatten()}         
  //ch_happy_runs.view()

  //
  // MODULE: Run happy
  //
  HAPPY_HAPPY (
    ch_happy_runs,
    ch_refgen_fasta
  )

  emit:
  HAPPY_HAPPY.out.metrics                
  //versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}