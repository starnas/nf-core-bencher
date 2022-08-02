//
// Run collectmultiplemetrics

//
// MODULE: Installed directly from nf-core/modules
//
include { MOSDEPTH } from '../../modules/nf-core/modules/mosdepth/main'
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main'
//include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'

workflow BAM_METRICS {
  take: 

  inputs_channel

  main:
  
  //
  // CODE: Prepare picard input channel
  //
  ch_multimetrics_runs = inputs_channel.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed -> tuple(meta, sample_bam).flatten()}                     
  //ch_multimetrics_runs.view()

  //
  // MODULE: collectmultiplemetrics
  //
  // PICARD_COLLECTMULTIPLEMETRICS (
  //   ch_multimetrics_runs,
  //   params.fasta, 
  //   params.fasta_fai
  // )
  
  //
  // MODULE: index from samtools
  //
  SAMTOOLS_INDEX (
    ch_multimetrics_runs
  )

  //
  // CODE: Prepare mosdepth input channel
  //  
  ch_mosdepth_bams = inputs_channel.combine(SAMTOOLS_INDEX.out.bai, by: 0)
  ch_mosdepth_bams = ch_mosdepth_bams.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed, sample_bai -> tuple(meta, sample_bam, sample_bai).flatten()}         
  ch_mosdepth_beds = inputs_channel.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed -> tuple(final_region_bed).flatten()}   

  //
  // MODULE: mosdepth
  //
  MOSDEPTH (
    ch_mosdepth_bams,
    ch_mosdepth_beds,
    params.fasta
  )

  emit:
  MOSDEPTH.out.summary_txt                
  //versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}