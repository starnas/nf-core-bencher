//
// Run collectmultiplemetrics

//
// MODULE: Installed directly from nf-core/modules
//
include { SAMTOOLS_INDEX } from '../../modules/nf-core/modules/samtools/index/main'
include { MOSDEPTH } from '../../modules/nf-core/modules/mosdepth/main'
include { QUALIMAP_BAMQC } from '../../modules/nf-core/modules/qualimap/bamqc/main'
include { PICARD_COLLECTMULTIPLEMETRICS } from '../../modules/nf-core/modules/picard/collectmultiplemetrics/main'

workflow BAM_METRICS {
  take: 

  inputs_channel

  main:
  
  //
  // CODE: Prepare picard input channel
  //
  ch_multimetrics_runs = inputs_channel.map{meta, ref_vcf, sample_vcf, bed_file, sample_bam -> tuple(meta, sample_bam).flatten()}                     
  //ch_multimetrics_runs.view()

  //
  // MODULE: collectmultiplemetrics
  //
  PICARD_COLLECTMULTIPLEMETRICS (
    ch_multimetrics_runs,
    params.fasta, 
    params.fasta_fai
  )

  //
  // MODULE: samtools index
  //
  SAMTOOLS_INDEX (
    ch_multimetrics_runs
  )

  tmp_list = SAMTOOLS_INDEX.out.bai

  //
  // CODE: Prepare mosdepth input channel
  //
  ch_mosdepth_runs = inputs_channel.combine(SAMTOOLS_INDEX.out.bai, by: 0)
  ch_mosdepth_runs = ch_mosdepth_runs.map{meta, ref_vcf, sample_vcf, bed_file, sample_bam, sample_bai -> tuple(meta, sample_bam).flatten()}         
  ch_mosdepth_beds = inputs_channel.map{meta, ref_vcf, sample_vcf, bed_file, sample_bam -> tuple(bed_file).flatten()}   
  //
  // MODULE: mosdepth
  //
  QUALIMAP_BAMQC (
    ch_mosdepth_runs,
    ch_mosdepth_beds
  )

  //versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}