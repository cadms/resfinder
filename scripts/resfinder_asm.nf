#!/usr/bin/env nextflow

python3 = "python3"
resfinder = "/home/projects/cge/people/rkmo/resfinder4/src/resfinder/run_resfinder.py"

params.input = './*.fa'
// params.indir = './'
// params.ext = '.fa'
params.outdir = '.'
params.species

println("Search pattern: $params.input")

infile_ch = Channel
              .fromPath("$params.input", followLinks: true)
              .map{ file -> tuple(file.baseName, file) }

process resfinder{

    cpus 1
    time '30m'
    memory '1 GB'
    clusterOptions '-V -W group_list=cge -A cge'
    executor "PBS"

    input:
    set sampleID, file(datasetFile) from infile_ch

    output:
    stdout result

    """
    set +u
    module unload perl
    source /home/projects/cge/apps/env/rf4_env/bin/activate
    module load ncbi-blast/2.8.1+
    $python3 $resfinder -acq --point -ifa $datasetFile -o '$params.outdir/$sampleID' -s '$params.species'
    """
}

/*
result.subscribe {
    println it
}
*/
