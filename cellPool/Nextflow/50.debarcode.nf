#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*****************************************************************************************
/ Channels feeding the wells
/****************************************************************************************/

wellsFile = file(Path.of("${params.analysisDir}", "${params.listOfWells}"))
wells = wellsFile.readLines()
wellCh = Channel.fromList(wells)

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
	mkdir -p !{params.analysisDir}/!{params.barcodeAnalysisDir}
	mkdir -p !{params.analysisDir}/!{params.qualityControlDir}/barcodes
	'''
}

//--------------------------
process PROFILEBCD {
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.csv'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '0.5h'
	cpus 1
	scratch true

	input:
	each well
	val greenFlag

	output:
	path '*.csv'
	val "Done with profiling the barcodes for another well ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail
	
	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/60.profileBarcodes.py \
		!{well} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process THRESHPOP {
	publishDir "${params.analysisDir}/${params.qualityControlDir}/barcodes", pattern: "${params.debarcodingParams.saveNotes}"
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.txt'
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.txt.gz'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/barcodes", pattern: '*.png'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/barcodes", pattern: '*.pdf'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	val greenFlag

	output:
	path '*.txt.gz'
	path '*.txt'
	path '*.png'
	path '*.pdf'
	val "Finished population level thresholding ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/62.organize_barcode_data.R \
		!{params.analysisDir}/!{params.cellPoolConfig}

	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/63.backbone_prot_gating.py \
		!{params.analysisDir}/!{params.cellPoolConfig}

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/65.refmt_barcode_metadata.py \
		!{params.analysisDir}/!{params.cellPoolConfig}

	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/64.population_level_epitope_thresholding.R \
		!{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process PREPKNOWNS {
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '6h'
	cpus 1
	scratch true

	input:
	each cut
	val greenFlag

	output:
	path '*.txt'
	val "${cut}", emit: cut

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/66.prepKnowns.py \
		!{cut} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process FITUNKNOWNS {
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '2h'
	cpus 1
	scratch true

	input:
	each cutInitTpl

	output:
	path '*.txt'
	val "Model fitting finished for this cut and initialization ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_py3
	set -euo pipefail

	python !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/67.fitUnknowns.py \
		!{cutInitTpl[0]} !{cutInitTpl[1]} !{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process MODELFITQC {
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.txt'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/barcodes", pattern: '*.pdf'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/barcodes", pattern: '*.png'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	val greenFlag

	output:
	path '*.txt'
	path '*.pdf'
	path '*.png'
	val "Saved iterations and initializations that yielded highest likelihood ...", emit: greenFlag

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/68.smrzFitPfmnc.R \
		!{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

//--------------------------
process CLUSTERBCD {
	publishDir "${params.analysisDir}/${params.barcodeAnalysisDir}", pattern: '*.txt.gz'
	publishDir "${params.analysisDir}/${params.qualityControlDir}/barcodes", pattern: '*.txt'
	containerOptions "-B ${params.analysisDir}"
	memory '16 GB'
	time '1h'
	cpus 1
	scratch true

	input:
	val greenFlag

	output:
	path '*.txt.gz'
	path '*.txt'

	shell:
	'''
	source /opt/conda/etc/profile.d/conda.sh
	set +euo pipefail #Switch off strict mode to activate conda
	conda activate cellPool_R
	set -euo pipefail

	Rscript !{params.reqdCellpoolContainerPath_LET_ME_BE}/cellpool_procedures/69.model_based_clustering.R \
		!{params.analysisDir}/!{params.cellPoolConfig}
	'''
}

/******************************************************************************************
/ Workflow
/*****************************************************************************************/

workflow {
	MKOUTPUTDIR()
	PROFILEBCD(wellCh, MKOUTPUTDIR.out.greenFlag)
	THRESHPOP(PROFILEBCD.out.greenFlag.collect())
	PREPKNOWNS(Channel.from(1..params.debarcodingParams.nCuts),
		THRESHPOP.out.greenFlag.collect())
	FITUNKNOWNS(PREPKNOWNS.out.cut
			.combine(Channel.from(1..params.debarcodingParams.nInit)))
	MODELFITQC(FITUNKNOWNS.out.greenFlag.collect())
	CLUSTERBCD(MODELFITQC.out.greenFlag.collect())
}

/******************************************************************************************
/ Completion handler
******************************************************************************************/
workflow.onComplete {
	log.info ( workflow.success ? '\nDone!' : '\nOops ... something went wrong ...' )
}




