#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*****************************************************************************************
/ Channels feeding the wells, fields and cycles
/****************************************************************************************/

wellsFile = file(Path.of("${params.analysisDir}", "${params.listOfWells}"))
wells = wellsFile.readLines()
wellCh = Channel.fromList(wells)

cycleCh = Channel.fromList(params.cycles)

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
	mkdir -p !{params.analysisDir}/!{params.mosaicDir}
	mkdir -p !{params.analysisDir}/!{params.mosaicDir}/!{params.mosaickingParams.mosaicIntermediateFilesDir}
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/mosaics/showStitches/
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/mosaics/showMosaics/
	'''
}

//--------------------------
process MOSAIC {
	publishDir "${params.analysisDir}/${params.mosaicDir}/", pattern: '*.tileCentroidCoords.txt'
	publishDir "${params.analysisDir}/${params.mosaicDir}/${params.mosaickingParams.mosaicIntermediateFilesDir}", pattern: '*.intmdt.*'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag

	output:
	path '*.txt'
	path '*.graphml'
	tuple val("${cycWellTpl[1]}"), val("${cycWellTpl[0]}"), emit: wellCycTpl

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	argsString="!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}"
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/30.estRefStitchCoords.py ${argsString}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/31.estStitchCoordsByGraphFusionAndCameraModel.py ${argsString}
	'''
}

//--------------------------
process REGISTER {
	publishDir "${params.analysisDir}/${params.mosaicDir}/", pattern: '*.coords.txt'

	publishDir "${params.analysisDir}/${params.mosaicDir}/${params.mosaickingParams.mosaicIntermediateFilesDir}", pattern: '*.intmdt.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '3 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	each wellCycTpl
	
	output:
	path '*.txt'
	tuple val("${wellCycTpl[0]}"), val("${wellCycTpl[1]}"), emit: wellCycTpl

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/33.registerCycles.py \
		!{wellCycTpl[0]} !{wellCycTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process STITCH {
	publishDir "${params.analysisDir}/${params.mosaicDir}", pattern: '*.npy'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/mosaics/showMosaics", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each wellCycTpl

	output:
	path '*.npy'
	path '*.png'
	tuple val("${wellCycTpl[0]}"), val("${wellCycTpl[1]}"), emit: wellCycTpl

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/35.makeStitched.py \
		!{wellCycTpl[0]} !{wellCycTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process TRIMMOSAICANDQC {
	publishDir "${params.analysisDir}/${params.mosaicDir}", pattern: '*.npy'
	publishDir "${params.analysisDir}/${params.mosaicDir}", pattern: '*.txt'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/mosaics/showStitches", pattern: '*.showOverlap.png'
	publishDir "${params.analysisDir}/${params.mosaicDir}", pattern: '*.trimmed.tileOverlapZone.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '3h'
	cpus 1
	scratch true

	input:
	each well

	output:
	val "${well}", emit: well
	val 'Summarizing QC results ...', emit: greenFlag
	path '*.png'
	path '*.npy'
	path '*.txt'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	argsString="!{well} !{params.analysisDir}/!{params.cellPoolConfig}"	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/36.trimMosaics.py ${argsString}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/37.markTileOverlaps.py ${argsString}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/38.showEdgeAlignments.py ${argsString}
	'''
}

//--------------------------
process INTMOSAICQC {
	executor 'local'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/mosaics", pattern: '*.png'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/mosaics", pattern: '*.pdf'
	containerOptions "-B ${params.analysisDir}"
	cpus 1

	input:
	val greenFlag

	output:
	path '*.png'
	path '*.pdf'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/39.qc.mosaic_tileCoords.R !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process MKVIDEO_STITCHES {
	executor 'local'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/mosaics", pattern: '*.mp4'
	containerOptions "-B ${params.analysisDir}"
	cpus 1

	input:
	each cycle
	val greenFlag

	output:
	path '*.mp4'

	shell:
	'''
	ffmpeg -f image2 -r 1 \
		-i !{params.analysisDir}/!{params.qualityControlDir}/mosaics/showStitches/%02d.cycle!{cycle}.showOverlap.png \
		-vcodec mpeg4 -y cycle!{cycle}.mp4
	'''
}

//--------------------------
process GETNORMINFO {
	executor 'local'
	publishDir "${params.analysisDir}/${params.mosaicDir}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	cpus 1

	input:
	val greenFlag

	output:
	path '*.txt'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/39.saveNormalizers.R !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	MKOUTPUTDIR()
	MOSAIC(cycleCh.combine(wellCh), MKOUTPUTDIR.out.greenFlag)	
	REGISTER(MOSAIC.out.wellCycTpl
			.groupTuple()
			.map{ it -> it[0] }
			.combine(cycleCh.filter{ it != 1 }) )
	STITCH(REGISTER.out.wellCycTpl
			.groupTuple()
			.map{ it -> it[0]}
			.combine(cycleCh))
	TRIMMOSAICANDQC(STITCH.out.wellCycTpl
			.groupTuple()
			.map{ it -> it[0]})
	INTMOSAICQC(TRIMMOSAICANDQC.out.greenFlag.collect())
	MKVIDEO_STITCHES(cycleCh, TRIMMOSAICANDQC.out.greenFlag.collect())
	GETNORMINFO(TRIMMOSAICANDQC.out.greenFlag.collect())
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




