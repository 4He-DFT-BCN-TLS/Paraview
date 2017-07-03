#!/bin/bash

INDIR=${1}
FPS=${2}

ffmpeg -r ${FPS} -f image2 -i ${INDIR}/denxyz.%04d.png -vcodec libx264 -preset ultrafast -crf 18 -pix_fmt yuv420p movie-3D.mp4
