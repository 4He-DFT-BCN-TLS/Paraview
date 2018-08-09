#!/bin/bash

INDIR=${1}
FIRST_FRAME_NR=0
FPS=15
ID_WIDTH=4

ffmpeg -start_number ${FIRST_FRAME_NR} -r ${FPS} -f image2 -i ${INDIR}/denxyz.%0${ID_WIDTH}d.png \
	-vcodec libx264 -preset ultrafast -crf 18 -pix_fmt yuv420p Movie-3D.mp4
