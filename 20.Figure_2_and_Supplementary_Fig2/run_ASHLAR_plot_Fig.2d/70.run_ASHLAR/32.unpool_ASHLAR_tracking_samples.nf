#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*****************************************************************************************
/ Channels feeding the wells, fields and cycles
/****************************************************************************************/

wellsFile = file(Path.of("${params.analysisDir}", 
	"${params.ashlarAnalysisDir}", "00.config", "wellsToShowTracking.txt"))
wells = wellsFile.readLines()
wellCh = Channel.fromList(wells)

wellOmeTiffCh = Channel.fromPath(Path.of("${params.analysisDir}", 
	"${params.ashlarAnalysisDir}", "20.mosaics", "*_maxShift_15_sigma_3.ome.tiff"))

/******************************************************************************************
/ Processes
/*****************************************************************************************/
process MKOUTPUTDIR {
	executor 'local'
	cpus 1
	containerOptions "-B ${params.analysisDir}"

	input:

	output:
	val 'Ok to proceed with analysis ...', emit: greenFlag

	shell:
	'''
	mkdir -p !{params.analysisDir}/!{params.ashlarAnalysisDir}/50.showMaskTransfer/images
	'''
}

//--------------------------
process SHOWRNDMSMPL {
	publishDir "${params.analysisDir}/${params.ashlarAnalysisDir}/50.showMaskTransfer/images", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	path wellOmeTiff
	//val greenFlag

	output:
	path '*.png'
	val "Save samples for another well ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.analysisDir}/!{params.scripts}/40.sampleMaskTransfer.py \
		!{wellOmeTiff} \
		!{params.analysisDir}/!{params.ashlarAnalysisDir}/!{params.ashlarConfig}
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	//MKOUTPUTDIR()
	SHOWRNDMSMPL(wellOmeTiffCh)//, MKOUTPUTDIR.out.greenFlag)
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




