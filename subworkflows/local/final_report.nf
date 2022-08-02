//
// Prepare the final report using jupyter
//

//
// MODULE: Installed directly from nf-core/modules
//
include { JUPYTERNOTEBOOK } from '../../modules/nf-core/modules/jupyternotebook/main'

workflow FINAL_REPORT {
  take:

  happy_results
  mosdepth_results

  main:
  
  //
  // CODE: Prepare preppy input channel
  //
  ch_results = happy_results.combine(mosdepth_results, by: 0)
  ch_results = ch_results.map{meta, happy_csvs, happy_json, mosdepth_txt -> tuple(happy_csvs, mosdepth_txt).flatten()}
  ch_results = ch_results.collect()
  
  //
  // MODULE: jupyternotebook
  //
  
  JUPYTERNOTEBOOK (
    Channel.value([[id: 'test'], file("assets/220711_testjupyter.ipynb")]),
    [:],
    Channel.value(file("assets/220711_testdata.txt"))
  )
  //ch_results.view()
  
  //versions = INPUT_PREP.out.versions // channel: [ versions.yml ]
}