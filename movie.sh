#!/bin/bash

INDIR=${1}
FPS=${2}

ffmpeg -start_number 120 -r ${FPS} -f image2 -i ${INDIR}/denxyz.%03d.png -vcodec libx264 -preset ultrafast -crf 18 -pix_fmt yuv420p movie-3D.mp4
