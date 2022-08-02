//
// Check input samplesheet and get the input channel ready
//

//
// MODULE: Installed directly from nf-core/modules
//
include { BEDTOOLS_SORT } from '../../modules/nf-core/modules/bedtools/sort/main'
include { BEDTOOLS_MERGE } from '../../modules/nf-core/modules/bedtools/merge/main'
include { BEDTOOLS_INTERSECT } from '../../modules/nf-core/modules/bedtools/intersect/main'

// required for handling input json
import groovy.json.JsonSlurper

workflow INPUT_PREP {
  take:
  
  inputsheet // file: /path/to/jsonfile

  main:
  
  //
  // CODE: Read json using groovy
  //

  // slurp input json using groovy standard library
  def jSlurp = new JsonSlurper()
  def inputFile = new File(inputsheet.toString())
  def inputs = jSlurp.parseText(inputFile.text)

  //
  // CODE: Prepare high confidence reference channel
  //

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

  //
  // CODE: Prepare all strata channel
  //

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
  } else {
    l_strata_beds = []
  }

  // generate the channel
  ch_strata_beds = Channel.from(l_strata_beds)
  //ch_strata_beds.view()

  //
  // CODE: Prepare custom stratifications channel
  //
  if (inputs.run_custom_strata){
    
    // generate the channel
    ch_custom_beds = Channel.from(inputs.custom_strata)

    // add it to the strata channel, since treatment is the same
    ch_strata_beds = ch_strata_beds.concat(ch_custom_beds) 
  }
  
  //
  // CODE: Prepare input files channel
  //

  // placeholder for inputs as list
  l_input_genomes = []
  l_input_targetted = []

  // populating the inputs list - loop over the nist id
  for (giab_id in inputs.input_files.keySet()){

    // check if it is a valid giab_id
    if (giab_id in ["NA12878", "NA24143", "NA24149", "NA24385", "NA24631", "NA24694", "NA24695"]){

      // loop over sample id
      for (sample_id in inputs.input_files[giab_id].keySet()){

        vcf_file = inputs.input_files[giab_id][sample_id][0]
        bam_file = inputs.input_files[giab_id][sample_id][1]
        bed_file = inputs.input_files[giab_id][sample_id][2]

        // check if vcf supplied
        if (vcf_file.endsWith('vcf.gz')){

          // check if bam supplied
          if (bam_file.endsWith('bam')){
            
            // check if bed supplied - if yes, targetted analysis
            if (bed_file && bed_file.endsWith('bed')){
              
              //  giab_id, sample_id. vcf, bam and target bed to the input list
              l_input_targetted << [giab_id, sample_id, vcf_file, bam_file, bed_file]

            }

            // if not, non-targetted analysis
            else {

              // append giab_id, sample_id. vcf and bam to the input list
              l_input_genomes << [giab_id, sample_id, vcf_file, bam_file]

            }
          }
        }
      }
    }
  }

  // make into channels
  ch_input_genomes = Channel.from(l_input_genomes)
  //ch_input_genomes.view()
  ch_input_targetted = Channel.from(l_input_targetted)
  //ch_input_targetted.view()

  //
  // CODE: Final mixing of the channels
  //

  // prepare channel for hc = high confidence
  ch_hc_genomes = ch_input_genomes.combine(ch_refgen_hc, by: 0)
  ch_hc_targetted = ch_input_targetted.combine(ch_refgen_hc, by: 0)
  //ch_hc_targetted.view()
  //ch_hc_genomes.view()

  // final strata channel with tuples of nist_id, ref_vcf, strata_bed and strata_name
  ch_refgen_strata = ch_refgen_loc.map{a,b -> tuple(a, get_file(b, "vcf.gz")).flatten()}
  ch_refgen_strata = ch_refgen_strata.combine(ch_strata_beds)
  ch_refgen_strata = ch_refgen_strata.map{a,b,c -> tuple(a, b, c, base_name(c)).flatten()}

  // prepare channel for strata vcfs
  ch_st_genomes = ch_input_genomes.combine(ch_refgen_strata, by: 0)
  ch_st_targetted = ch_input_targetted.combine(ch_refgen_strata, by: 0)
  //ch_st_targetted.view()
  //ch_st_genomes.view()

  // concat the genome and targetted channles
  ch_genomes = ch_hc_genomes.concat(ch_st_genomes)
  ch_targetted = ch_hc_targetted.concat(ch_st_targetted)
  //ch_targetted.view()
  //ch_genomes.view()

  // prepare final structure for channels
  ch_genomes = ch_genomes.map{giab_id, sample_id, sample_vcf, sample_bam, ref_vcf, ref_bed_file, bed_id -> tuple(['id': sample_id+'_'+bed_id, 'sample_id': sample_id, 'ref_id': giab_id, 'target_bed': false, 'strata_bed_id': bed_id, 'target_bed_id': 'None'], ref_vcf, sample_vcf, ref_bed_file, sample_bam, 'None', ref_bed_file).flatten()}
  ch_targetted = ch_targetted.map{giab_id, sample_id, sample_vcf, sample_bam, sample_bed, ref_vcf, ref_bed_file, bed_id -> tuple(['id': sample_id+'_'+base_name(sample_bed)+'_vs_'+bed_id, 'sample_id': sample_id, 'ref_id': giab_id, 'target_bed': true, 'strata_bed_id': bed_id, 'target_bed_id': base_name(sample_bed)], ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed).flatten()}
  
  //
  // CODE: Intersecting the target bed with reference bed for the targetted channels
  // MODULE: bedtools_intersect
  //
  ch_bedtools = ch_targetted.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed -> tuple(meta, ref_bed_file, sample_bed).flatten()}

  BEDTOOLS_INTERSECT (
    ch_bedtools,
    Channel.value('bed')
  )
  
  //
  // MODULE: bedtools_sort
  //
  BEDTOOLS_SORT (
    BEDTOOLS_INTERSECT.out.intersect,
    Channel.value('sorted.bed')
  )

  //
  // MODULE: bedtools_merge
  //
  BEDTOOLS_MERGE (
    BEDTOOLS_SORT.out.sorted,
  )
  
  ch_targetted = ch_targetted.combine(BEDTOOLS_MERGE.out.bed, by: 0)
  ch_targetted = ch_targetted.map{meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed -> tuple(meta, ref_vcf, sample_vcf, ref_bed_file, sample_bam, sample_bed, final_region_bed).flatten()}

  //
  // CODE: Joining the targetted and genomes into one output channel
  //

  ch_formatted_inputs = ch_genomes.concat(ch_targetted) 
  //ch_formatted_inputs.view()

  emit:
  ch_formatted_inputs                
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