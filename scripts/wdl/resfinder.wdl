workflow Resistance {
    File inputSamplesFile
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)

		String outputDir

    scatter (sample in inputSamples) {
        call ResFinder {
            input: fq1=sample[0],
                fq2=sample[1],
                species=sample[2],
								outputRoot=outputDir
        }
    }
}

task ResFinder {

    String fq1
    String fq2
    String species
		String outputRoot
    String filename = basename(fq1)
    String out_dir_name = "${outputRoot}/${filename}.rf_out"

    command {
        mkdir ${out_dir_name}

				if [ ${species} = 'other' ]
				then
					$CGE_PYTHON $CGE_RESFINDER \
							-ifq ${fq1} ${fq2} \
							--blastPath $CGE_BLASTN \
							--kmaPath $CGE_KMA \
							--species "${species}" \
							--db_path_res $CGE_RESFINDER_RESGENE_DB \
							--acquired \
							--acq_overlap 30 \
							--min_cov $CGE_RESFINDER_GENE_COV \
							--threshold $CGE_RESFINDER_GENE_ID \
							-o ${out_dir_name}
				else
		        $CGE_PYTHON $CGE_RESFINDER \
		            -ifq ${fq1} ${fq2} \
		            --blastPath $CGE_BLASTN \
		            --kmaPath $CGE_KMA \
		            --species "${species}" \
		            --db_path_res $CGE_RESFINDER_RESGENE_DB \
		            --acquired \
		            --acq_overlap 30 \
		            --min_cov $CGE_RESFINDER_GENE_COV \
		            --threshold $CGE_RESFINDER_GENE_ID \
		            --point \
		            --db_path_point $CGE_RESFINDER_RESPOINT_DB \
		            --min_cov_point $CGE_RESFINDER_POINT_COV \
		            --threshold_point $CGE_RESFINDER_POINT_ID \
		            -o ${out_dir_name}
				fi
    }
    output {
        File rf_out = "${out_dir_name}/std_format_under_development.json"
    }
    runtime {
        walltime: "1:00:00"
        cpu: 2
        memory_mb = 3900.0
        queue = "cge"
    }
}
