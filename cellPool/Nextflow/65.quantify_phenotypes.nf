#!/usr/bin/env nextflow

import java.nio.file.Paths
nextflow.enable.dsl = 2

/*****************************************************************************************
/ Channels feeding the wells, fields and cycles
/****************************************************************************************/

wellsFile = file(Path.of("${params.analysisDir}", "${params.listOfWells}"))
wells = wellsFile.readLines()
wellCh = Channel.fromList(wells)

cycleCh = Channel.fromList(params.cycles)
cycWellCh = cycleCh.combine(wellCh)

/******************************************************************************************
/ Processes
/*****************************************************************************************/
process GET_PHENO {
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.featureDirSuffix}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each well

	output:
	path '*.txt'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.analysisDir}/!{params.phenotypeProfilingParams.script} \
		!{well} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	GET_PHENO(wellCh)
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




