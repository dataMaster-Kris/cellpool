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
process MKOUTPUTDIR {
	executor 'local'
	cpus 1
	containerOptions "-B ${params.analysisDir}"

	input:

	output:
	val 'Ok to proceed with analysis ...', emit: greenFlag

	shell:
	'''
	for i in {!{params.segmentationParams.trainingDataDirSuffix},!{params.segmentationParams.trainedModelDirSuffix},!{params.segmentationParams.maskDirSuffix},!{params.segmentationParams.featureDirSuffix},!{params.segmentationParams.unpoolDirSuffix}}; do \
		mkdir -p !{params.analysisDir}/!{params.segmentationDir}/iter!{params.segmentationParams.segIter}${i}; done

	mkdir -p !{params.analysisDir}/!{params.segmentationDir}/iter!{params.segmentationParams.segIter}!{params.segmentationParams.unpoolDirSuffix}/!{params.segmentationParams.unpooledImagesDir}
	mkdir -p !{params.analysisDir}/!{params.segmentationDir}/iter!{params.segmentationParams.segIter}!{params.segmentationParams.unpoolDirSuffix}/!{params.segmentationParams.unpooledNominationsDir}
	'''
}

//--------------------------
process GET_TRAINING_DATA {
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.trainingDataDirSuffix}", pattern: '*.tiff'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag

	output:
	path '*.tiff'
	tuple val("${cycWellTpl[0]}"), val("${cycWellTpl[1]}"), emit: cycWellTpl

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	argsString="!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}"
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/40.step0.nominate_for_training.py ${argsString}
	'''
}

//--------------------------
process SEGMENT {
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.maskDirSuffix}", pattern: '*.tif'
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.maskDirSuffix}", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '30 GB'
	time '3h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag
	
	output:
	path '*'
	tuple val("${cycWellTpl[0]}"), val("${cycWellTpl[1]}"), emit: cycWellTpl

	shell:
	'''
	NUMBA_CACHE_DIR='\${TMPDIR:?}'
	export NUMBA_CACHE_DIR
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_cellpose
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/40.step2.segment_with_trained_model.C4.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process GET_FEATURES {
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.featureDirSuffix}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cycWellTpl

	output:
	path '*.txt'
	tuple val("${cycWellTpl[1]}"), val("${cycWellTpl[0]}"), emit: wellCycTpl

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/40.step3.makeFeatureTbl.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process SAMPLE4QC {
	executor 'local'
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.unpoolDirSuffix}/${params.segmentationParams.unpooledNominationsDir}", pattern: '*.nominations.txt'
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.unpoolDirSuffix}", pattern: 'nominations.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	val wells

	output:
	path '*.txt'
	val 'Sampled segmentation masks for QC. Proceed with unpooling ...', emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/40.step4.a.nominateForUnpooling.R \
		!{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process UNPOOL_SEGQC {
	publishDir "${params.analysisDir}/${params.segmentationDir}/iter${params.segmentationParams.segIter}${params.segmentationParams.unpoolDirSuffix}/${params.segmentationParams.unpooledImagesDir}", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag

	output:
	path '*.png'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/40.step4.b.saveUnpooled.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	MKOUTPUTDIR()

	segModel = Paths.get(params.analysisDir,
		params.segmentationDir, 
		"iter" + params.segmentationParams.segIter + 
			params.segmentationParams.trainedModelDirSuffix, 
		params.segmentationParams.cellposeParams.model)

	if (segModel.exists()) {
		println "Segmentation model found. Proceeding to segmentation."
		SEGMENT(cycWellCh, MKOUTPUTDIR.out.greenFlag)
		GET_FEATURES(SEGMENT.out.cycWellTpl)
		SAMPLE4QC(GET_FEATURES.out
				.wellCycTpl
				.groupTuple()
				.map{ it -> it[0] }
				.unique()
				.toList())
		UNPOOL_SEGQC(cycleCh.combine(GET_FEATURES.out
					.wellCycTpl
					.groupTuple()
					.map{ it -> it[0] }
					.unique()),
				SAMPLE4QC.out.greenFlag)
	} else {
		println "Segmentation model not found."
		println "Make sure that you added your model file in the segmentation directory." 
		println "I am saving image patches for you to train a Cellpose model."
		GET_TRAINING_DATA(
			cycWellCh.randomSample(params.segmentationParams.nWellsForTrainingData), 
			MKOUTPUTDIR.out.greenFlag)
		println "When you have a trained model, you can resume this pipeline."
	}	
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




