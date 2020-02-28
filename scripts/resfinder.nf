#!/usr/bin/env nextflow

python3 = "python3"
resfinder = "/home/projects/cge/apps/resfinder/resfinder/run_resfinder.py"

params.indir = './'
params.ext = '.fq.gz'
params.outdir = '.'
params.species

println("Search pattern: $params.indir*{1,2}$params.ext")

Channel
    .fromFilePairs("$params.indir*{1,2}$params.ext", followLinks: true)
    .set{ infile_ch }

process resfinder{

    cpus 5
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
    $python3 $resfinder -acq --point -ifq $datasetFile -o '$params.outdir/$sampleID' -s '$params.species'
    """
}

/*
result.subscribe {
    println it
}
*/
