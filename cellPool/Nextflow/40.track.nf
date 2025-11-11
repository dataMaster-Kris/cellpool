#!/usr/bin/env nextflow

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
	mkdir -p !{params.analysisDir}/!{params.trackObjsDir}
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/tracking/stage1
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/tracking/stage2
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/tracking/stage3
	'''
}

//--------------------------
process TRACKSTG1 {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag

	output:
	path '*.txt'
	val "Done with tracking stage 1 for another well-cycle ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/50.trackNuclei.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process GETNGBRS {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag

	output:
	path '*.txt'
	val "Saved neighbor stats for another well-cycle ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/52.findNeighbors.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process FISHERRORANDQC {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage1", pattern: '*.png'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage1", pattern: '*.pdf'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2.5h'
	cpus 1
	scratch true

	input:
	each cycle
	val greenFlag1
	val greenFlag2

	output:
	path '*.png'
	path '*.txt'
	path '*.pdf'
	val "${cycle}", emit: cycle

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/54.fishErrors.R \
		!{cycle} !{params.analysisDir}/!{params.cellPoolConfig}
	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/55.qc.tracking.R \
		!{cycle} 1 !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process TRACKSTG2 {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt'
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
	val "Done with tracking stage 2 for another well-cycle ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/56.corrTrackNuclei.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process TRACKINGQC {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage2", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	each cycle
	val greenFlag

	output:
	path '*.png'
	path '*.txt'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/57.qc.corrTracking.R \
		!{cycle} 2 !{params.analysisDir}/!{params.cellPoolConfig}
	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/55.qc.tracking.R \
		!{cycle} 2 !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process SMPLTRKMATES {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage2/${params.unpoolTrackingQcParams.showUnpoolDir}", pattern: '*.png', failOnError: false
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage2/${params.unpoolTrackingQcParams.showUnpoolDir}", pattern: '*.allLabels.txt', failOnError: false
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.highConfValidMates.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each well

	output:
	path '*.png'
	path '*.allLabels.txt'
	path '*.highConfValidMates.txt'
	val "Save samples for another well ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.a.sampleTrackMates.py \
		!{well} Rndm stg2OutSuffix !{params.analysisDir}/!{params.cellPoolConfig}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.a.sampleTrackMates.py \
		!{well} LowConf stg2OutSuffix !{params.analysisDir}/!{params.cellPoolConfig}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.a.sampleTrackMates.py \
		!{well} HighConf stg2OutSuffix !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process MKMP4VIDEOS {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage2", pattern: '*.mp4'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	val greenFlag

	output:
	path '*.mp4'
	
	shell:
	'''
	now=`date '+%F_%H%M%S'`
	for smpl in {Rndm,LowConf,HighConf}; do \
		ffmpeg -f image2 -r 1 \
			-i !{params.analysisDir}/!{params.qualityControlDir}/tracking/stage2/!{params.unpoolTrackingQcParams.showUnpoolDir}/%02d.${smpl}.png \
			-vcodec mpeg4 \
			-y showing${smpl}Mates.${now}.mp4; \
	done
	'''
}

//--------------------------
process MKAVIVIDEOS {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage2", pattern: '*.avi'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	each batch
	val greenFlag

	output:
	path '*.avi'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_cellpose
	set -euo pipefail
	
	for smpl in {Rndm,LowConf,HighConf}
	do
		python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.b.videoTrackMates.py \
			!{batch} ${smpl} !{params.analysisDir}/!{params.cellPoolConfig}
	done
	'''
}

//--------------------------
process TRACKSTG3 {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt', failOnError: false
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cycWellTpl
	val greenFlag
	val greenFlag2

	output:
	path '*.txt'
	tuple val("${cycWellTpl[1]}"), val("${cycWellTpl[0]}"), emit: wellCycTpl
	val "Done with tracking stage 2 for another well-cycle ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.d.trackUsingLocalRefs.py \
		!{cycWellTpl[0]} !{cycWellTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process TRACKINGQC1 {
	publishDir "${params.analysisDir}/${params.trackObjsDir}", pattern: '*.txt'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage3", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	each cycle
	val greenFlag

	output:
	path '*.png'
	path '*.txt'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/57.qc.corrTracking.R \
		!{cycle} 3 !{params.analysisDir}/!{params.cellPoolConfig}
	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/55.qc.tracking.R \
		!{cycle} 3 !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process SMPLTRKMATES1 {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage3/${params.unpoolTrackingQcParams.showUnpoolDir}", pattern: '*.png', failOnError: false
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage3/${params.unpoolTrackingQcParams.showUnpoolDir}", pattern: '*.allLabels.txt', failOnError: false
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each well

	output:
	path '*.png'
	path '*.allLabels.txt'
	val "Save samples for another well ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.a.sampleTrackMates.py \
		!{well} Rndm stg3OutSuffix !{params.analysisDir}/!{params.cellPoolConfig}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.a.sampleTrackMates.py \
		!{well} LowConf stg3OutSuffix !{params.analysisDir}/!{params.cellPoolConfig}
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.a.sampleTrackMates.py \
		!{well} HighConf stg3OutSuffix !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process MKMP4VIDEOS1 {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage3", pattern: '*.mp4'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	val greenFlag

	output:
	path '*.mp4'
	
	shell:
	'''
	now=`date '+%F_%H%M%S'`
	for smpl in {Rndm,LowConf,HighConf}; do \
		ffmpeg -f image2 -r 1 \
			-i !{params.analysisDir}/!{params.qualityControlDir}/tracking/stage3/!{params.unpoolTrackingQcParams.showUnpoolDir}/%02d.${smpl}.png \
			-vcodec mpeg4 \
			-y showing${smpl}Mates.${now}.mp4; \
	done
	'''
}

//--------------------------
process MKAVIVIDEOS1 {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/tracking/stage3", pattern: '*.avi'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	each batch
	val greenFlag

	output:
	path '*.avi'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_cellpose
	set -euo pipefail
	
	for smpl in {Rndm,LowConf,HighConf}
	do
		python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/59.b.videoTrackMates.py \
			!{batch} ${smpl} !{params.analysisDir}/!{params.cellPoolConfig}
	done
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	MKOUTPUTDIR()
	if (params.trackingParams.doStg1AndStg2) {
		TRACKSTG1(cycleCh
				.filter{ it != 1 }
				.combine(wellCh), 
			MKOUTPUTDIR.out.greenFlag)

		GETNGBRS(cycWellCh, 
			MKOUTPUTDIR.out.greenFlag)
		FISHERRORANDQC(cycleCh
				.filter{ it != 1 },
			TRACKSTG1.out.greenFlag.collect(), 
			GETNGBRS.out.greenFlag.collect())
		TRACKSTG2(FISHERRORANDQC.out.cycle
				.combine(wellCh))
	}

	if (params.trackingParams.doStg2QC) {
		TRACKINGQC(cycleCh
				.filter{ it != 1 },
			TRACKSTG2.out.greenFlag.collect()) 
		SMPLTRKMATES(TRACKSTG2.out.wellCycTpl
				.groupTuple()
				.map{ it -> it[0]})
		if (params.unpoolTrackingQcParams.videoFormat.mp4) 
			MKMP4VIDEOS(SMPLTRKMATES.out.greenFlag.collect())
		if (params.unpoolTrackingQcParams.videoFormat.avi)
			MKAVIVIDEOS(
				Channel.from(1..params.unpoolTrackingQcParams.videoFormat.nBatchesOfAvi),
				SMPLTRKMATES.out.greenFlag.collect())
	}

	if (params.trackingParams.localRefParams.validatedSetOk) { 
		TRACKSTG3(cycleCh
				.filter{ it != 1 }
				.combine(wellCh),
			TRACKSTG2.out.greenFlag.collect(), 
			SMPLTRKMATES.out.greenFlag.collect())
		TRACKINGQC1(cycleCh
				.filter{ it != 1 },
			TRACKSTG3.out.greenFlag.collect())
		SMPLTRKMATES1(TRACKSTG3.out.wellCycTpl
				.groupTuple()
				.map{ it -> it[0]})
		if (params.unpoolTrackingQcParams.videoFormat.mp4) 
			MKMP4VIDEOS1(SMPLTRKMATES1.out.greenFlag.collect())
		if (params.unpoolTrackingQcParams.videoFormat.avi)
			MKAVIVIDEOS1(
				Channel.from(1..params.unpoolTrackingQcParams.videoFormat.nBatchesOfAvi),
				SMPLTRKMATES1.out.greenFlag.collect())
	}
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




