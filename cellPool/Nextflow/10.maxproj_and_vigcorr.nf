#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*****************************************************************************************
/ Channels feeding the wells, fields and cycles
/****************************************************************************************/

wellsFile = file(Path.of("${params.analysisDir}", "${params.listOfWells}"))
wells = wellsFile.readLines()
wellCh = Channel.fromList(wells)

cycleCh = Channel.fromList(params.cycles)

cycFldCh = Channel.fromPath(Path.of("${params.analysisDir}", "${params.wellMapFile}"))
		.splitCsv(header : true, sep : '\t')
		.map { row -> tuple(row.cycle, row.FieldID) }

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
	mkdir -p !{params.analysisDir}/!{params.maxProjDir}
	mkdir -p !{params.analysisDir}/!{params.vigCorrectDir}
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/vignetting_correction
	'''
}

//----------------------------
process MAXPROJ {
	publishDir "${params.analysisDir}/${params.maxProjDir}"
	containerOptions "-B ${params.rawDataDir} -B ${params.analysisDir}"
	memory '1 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cycFldTpl
	val greenFlag

	output:
	path '*.tiff', emit: cycFldFiles

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/10.maxProj.py \
		!{cycFldTpl[0]} !{cycFldTpl[1]} \
		!{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//-----------------------------
process VIGCORRECTPY {
	publishDir "${params.analysisDir}/${params.vigCorrectDir}", pattern: "*.vigCorr.tiff"
	publishDir "${params.analysisDir}/${params.qualityControlDir}/vignetting_correction", pattern: "*-field.tiff"
	containerOptions "-B ${params.analysisDir}"
	memory '15 GB'
	time '8h'
	cpus 1

	input:
	path cycFldFiles

	output:
	path '*.tiff'
		
	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/20.vigCorrect.py !{cycFldFiles}
	'''
}

//-----------------------------
process VIGCORRECTM {
	publishDir "${params.analysisDir}/${params.vigCorrectDir}", pattern: "*.vigCorr.tiff"
	publishDir "${params.analysisDir}/${params.qualityControlDir}/vignetting_correction", pattern: "*-field.tiff"
	container false
	memory '15 GB'
	time '10h'
	cpus 1
	scratch true

	input:
	path cycFldFiles

	output:
	path '*.tiff'
		
	shell:
	'''
	if [[ !{params.matlabModuleName} != 'None' ]] ; then
		module load "!{params.matlabModuleName}"
	fi

	apptainer exec !{params.containerDir}/cellpool.sif \
		cp -R /app/BaSiC-master .
	apptainer exec !{params.containerDir}/cellpool.sif \
		cat /app/cellpool_procedures/20.vigCorrect.m | \
		matlab -nodesktop -nosplash \
		-r "files='!{cycFldFiles}';"
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	MKOUTPUTDIR()
	MAXPROJ(cycFldCh, MKOUTPUTDIR.out.greenFlag)
	
	if (params.BaSiCImplementation.MATLAB)
		VIGCORRECTM(MAXPROJ.out.cycFldFiles)
	else 
		VIGCORRECTPY(MAXPROJ.out.cycFldFiles)
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




