This repository contains the code for the manuscript entitled "Scalable imaging-based profiling of CRISPR perturbations with protein barcodes" by Krishna Choudhary and Michael McManus. All code in this repository was developed by Krishna Choudhary as part of his postdoctoral research at the Michael McManus lab.

To reproduce figures in the manuscript, download the cellPool Docker container from `docker://krriisshhna/cellpool:latest` or build using the Docker file in the cellPool folder here. Then, edit the paths for raw data and output directories in the `00.config/config.json` folder. Finally, run the .nf files in the numerical order following the Nextflow commands, e.g.,

```
docker pull docker://krriisshhna/cellpool:latest
nextflow run 10.maxproj_and_vigcorr.nf -with-singularity cellpool.sif -params-file 00.config/config.json
```

See [Nextflow documentation](https://www.nextflow.io/docs/latest/index.html) for more parameter options. 
