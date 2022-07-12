//
// Check input samplesheet and get read channels
//

import groovy.json.JsonSlurper

workflow INPUT_PREP {
  take:
  
  inputsheet // file: /path/to/inputsheet.csv

  main:
  
  // slurp using groovy standard library
  def jSlurp = new JsonSlurper()
  def inputFile = new File(inputsheet.toString())
  def inputs = jSlurp.parseText(inputFile.text)
  print inputs.input_files

  // placeholders for the 2 channels
  l_bams = []
  l_vcfs = []

  // iterate over inputs, check them and split to bams and vcfs
  inputs.input_files.each {

    // check for valid identifier
    if (it[0].startsWith('NA')){

      // check if vcf supplied
      if (it[1].endsWith('vcf.gz')){

        // check if it exists
        File file = new File(it[1])
        if (file.exists()){

          // append to vcf channel
          l_vcfs << [it[0], it[1]]

        }

        // check if bam supplied
        if (it.size() > 2){

          // check if it is bam :)
          if (it[2].endsWith('.bam')){

            // check if it exists
            if (file.exists()){

              // append to vcf channel
              l_bams << [it[0], it[2]]

            }
          }
        }
      }
    }
  }

  //----------------------------------------------------------------------
  // high confidence channel preparations

  // channel with reference locations based on NIST ftp download
  ch_refgen_loc = Channel.from(
    ["NA12878", "${params.nist_release}/NA12878_HG001/${params.nist_version}/${params.genome_version}"],
    ["NA24143", "${params.nist_release}/AshkenazimTrio/HG004_NA24143_mother/${params.nist_version}/${params.genome_version}"],
    ["NA24149", "${params.nist_release}/AshkenazimTrio/HG003_NA24149_father/${params.nist_version}/${params.genome_version}"],
    ["NA24385", "${params.nist_release}/AshkenazimTrio/HG002_NA24385_son/${params.nist_version}/${params.genome_version}"],
    ["NA24631", "${params.nist_release}/ChineseTrio/HG005_NA24631_son/${params.nist_version}/${params.genome_version}"],
    ["NA24694", "${params.nist_release}/ChineseTrio/HG006_NA24694_father/${params.nist_version}/${params.genome_version}"],
    ["NA24695", "${params.nist_release}/ChineseTrio/HG007_NA24695_mother/${params.nist_version}/${params.genome_version}"]
  )

  // final high confidence channel with tuples of nist_id, ref_vcf, strata_bed and strata_name
  ch_refgen_hc = ch_refgen_loc.map{a,b -> tuple(a, get_file(b, "vcf.gz"), get_file(b, "bed")).flatten()}
  ch_refgen_hc = ch_refgen_hc.combine(Channel.from("${params.genome_version}_high_confidence"))

  // prepare channel for hc vcfs
  ch_inputs_hc = ch_vcfs.combine(ch_refgen_hc, by: 0)
  // ch_inputs_hc.view()



  emit:
  l_vcfs                                     // channel: [ val(meta), [ reads ] ]
  versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}
