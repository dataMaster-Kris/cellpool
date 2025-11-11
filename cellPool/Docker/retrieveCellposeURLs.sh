#!/bin/bash --login

source /opt/conda/etc/profile.d/conda.sh
set +euo pipefail
conda activate cellPool_cellpose
set -euo pipefail
 
python retrieveCellposeURLs.py
