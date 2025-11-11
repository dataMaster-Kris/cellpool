#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

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
	mkdir -p !{params.analysisDir}/!{params.unpooledPhenoDir}
	mkdir -p !{params.analysisDir}/!{params.unpooledPhenoDir}/images
	'''
}

//--------------------------
process LIST2UNPOOL {
	executor 'local'
	publishDir "${params.analysisDir}/${params.unpooledPhenoDir}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	val greenFlag

	output:
	path '*.txt', emit: allLists

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail
	
	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/70.nominate_for_unpooling.R \
		!{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process UNPOOLPHENO {
	publishDir "${params.analysisDir}/${params.unpooledPhenoDir}/images", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	each objList

	output:
	path '*.png'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/72.unpool.py \
		!{objList} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	MKOUTPUTDIR()
	LIST2UNPOOL(MKOUTPUTDIR.out.greenFlag)
	UNPOOLPHENO(LIST2UNPOOL.out.allLists)
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




