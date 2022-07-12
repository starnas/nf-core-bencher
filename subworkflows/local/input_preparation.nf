//
// Check input samplesheet and get read channels
//

import groovy.json.JsonSlurper

workflow INPUT_PREP {
  take:
  
  inputsheet // file: /path/to/jsonfile

  main:
  
  //----------------------------------------------------------------------
  // prepare reference genome channeæ

  ch_refgen_fasta = Channel.from([[params.fasta, params.fasta_fai]])
  
  //----------------------------------------------------------------------
  // read the input json 

  // slurp input json using groovy standard library
  def jSlurp = new JsonSlurper()
  def inputFile = new File(inputsheet.toString())
  def inputs = jSlurp.parseText(inputFile.text)

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
  //ch_refgen_loc.view()

  // final high confidence channel with tuples of nist_id, ref_vcf, strata_bed and strata_name
  ch_refgen_hc = ch_refgen_loc.map{a,b -> tuple(a, get_file(b, "vcf.gz"), get_file(b, "bed")).flatten()}
  ch_refgen_hc = ch_refgen_hc.combine(Channel.from("${params.genome_version}_high_confidence"))
  //ch_refgen_hc.view()

  //----------------------------------------------------------------------
  // all strata channel preparations

  // all strata list based on folders from nist website
  def l_strata_beds = []
  new File("${params.nist_release}/genome-stratifications/${params.genome_strata_version}/${params.genome_version}").traverse(type: groovy.io.FileType.FILES, nameFilter: ~/.*.bed.gz/) {it -> l_strata_beds.add(it)}

  // filter out GenomeSpecificskoda greenline
  l_strata_beds.removeAll{it =~/GenomeSpecific/}

  // filter out anything not in subset if defined
  if (inputs.run_filter_strata){
    def l_strata_beds_filtered = []
    for (b in inputs.filter_strata){
      l_strata_beds_cp = l_strata_beds
      l_strata_beds_filtered.add(l_strata_beds_cp.findAll{it =~/${b}/})
    }
    l_strata_beds = l_strata_beds_filtered.flatten()
  }

  // generate the channel
  ch_strata_beds = Channel.from(l_strata_beds)
  //ch_strata_beds.view()

  //----------------------------------------------------------------------
  // custom strata channel preparations

  // generate the channel
  ch_custom_beds = Channel.from(inputs.custom_strata)
  ch_custom_beds.view()
  
  //----------------------------------------------------------------------
  // prepare the input files channel

  // placeholder for inputs as list
  l_input_files = []

  // populating the inputs list - loop over the nist id
  for (giab_id in inputs.input_files.keySet()){

    // check if it is a valid giab_id
    if (giab_id in ["NA12878", "NA24143", "NA24149", "NA24385", "NA24631", "NA24694", "NA24695"]){

      // loop over sample id
      for (sample_id in inputs.input_files[giab_id].keySet()){

        vcf_file = inputs.input_files[giab_id][sample_id][0]
        bam_file = inputs.input_files[giab_id][sample_id][1]

        // check if vcf supplied
        if (vcf_file.endsWith('vcf.gz')){

          // check if bam supplied
          if (bam_file.endsWith('bam')){

            // append giab_id, sample_id. vcf and bam to the input list
            l_input_files << [giab_id, sample_id, vcf_file, bam_file]

          }
        }
      }
    }
  }
  ch_input_files = Channel.from(l_input_files)
  //ch_input_files.view()

  //----------------------------------------------------------------------
  // final mixing of channels

  // prepare channel for hc = high confidence
  ch_high_confidence = ch_input_files.combine(ch_refgen_hc, by: 0)
  //ch_high_confidence.view()

  // final strata channel with tuples of nist_id, ref_vcf, strata_bed and strata_name
  //ch_refgen_strata = ch_refgen_loc.map{a,b -> tuple(a, get_file(b, "vcf.gz")).flatten()}
  //ch_refgen_strata = ch_refgen_strata.combine(ch_strata_beds)
  //ch_refgen_strata = ch_refgen_strata.map{a,b,c -> tuple(a, b, c, base_name(c)).flatten()}

  // prepare channel for strata vcfs
  //ch_inputs_strata = ch_vcfs.combine(ch_refgen_strata, by: 0)
  //ch_inputs_strata.view()

  // final strata channel with tuples of nist_id, ref_vcf, strata_bed and strata_name
  //ch_refgen_custom = ch_refgen_loc.map{a,b -> tuple(a, get_file(b, "vcf.gz")).flatten()}
  //ch_refgen_custom = ch_refgen_custom.combine(ch_custom_beds)
  //ch_refgen_custom = ch_refgen_custom.map{a,b,c -> tuple(a, b, c, base_name(c)).flatten()}

  // prepare channel for strata vcfs
  //ch_inputs_custom = ch_vcfs.combine(ch_refgen_custom, by: 0)
  //ch_inputs_custom.view()

  // combined channel with all runs
  ch_preppy_runs = ch_high_confidence.map{giab_id, sample_id, sample_vcf, sample_bam, ref_vcf, bed_file, bed_id -> tuple(['id': sample_id], sample_vcf, bed_file).flatten()}
  //ch_preppy_runs.view()

  emit:
  ch_preppy_runs                                     
  ch_refgen_fasta
  //versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}


// function to find a file based on extension
def get_file (ref_folder, ext){
  file(ref_folder + "/*" + ext)
}

// function to get simple name
def base_name (file_path){
  file(file_path).getBaseName()
}