import cv2
import os
import sys
import json
from datetime import datetime

now = datetime.now()
date_time_str = now.strftime("%Y-%m-%d_%H%M%S")

batch = int(sys.argv[1])
smpl = sys.argv[2]
paramsFile = sys.argv[3]

with open(paramsFile, 'r') as inFile:
	params = json.load(inFile)

nPerBatch = params['unpoolTrackingQcParams']['videoFormat']['nWellsPerAvi']
inDir = os.path.join(params['analysisDir'], params['qualityControlDir'], \
	'tracking', params['unpoolTrackingQcParams']['showUnpoolDir'])
images = [img for img in os.listdir(inDir) if img.endswith('.' + smpl + '.png') & \
		 (int(img[:2]) <= batch*nPerBatch) & (int(img[:2]) > nPerBatch * (batch - 1))]
images = sorted(images)
frame = cv2.imread(os.path.join(inDir, images[0]))
height, width, layers = frame.shape

video = cv2.VideoWriter('showing' + smpl + 'Mates.' + str(batch) + date_time_str + '.avi', \
				0, 1, (width, height))

for image in images:
    video.write(cv2.imread(os.path.join(inDir, image)))

video.release()

