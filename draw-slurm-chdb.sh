#!/bin/bash
#source /opt/intel/bin/compilervars.sh intel64

# Call: ./draw-slurm-chdb.sh ../path/to/density.xxxx.dat output_directory

dt=0.5
ID=${1##*density.}   # density.0001.dat => 0001.dat
ID=${ID%.dat}      # 0001.dat         => 0001
INDIR=${2}
OUTDIR=${3}

(cd ${OUTDIR}; mkdir denxyz_py impurity readwf_dat readwf_err vtk xyz3Dplot 2> /dev/null )

printf -v t "%.1f" $(echo "${ID}*${dt}"|bc)
cp readwf.dat readwf.${ID}.dat
sed -i "/File10=/c\ File10='${INDIR}/density.${ID}.dat'" readwf.${ID}.dat
sed -i "/File31=/c\ File31='position.${ID}.dat'" readwf.${ID}.dat
sed -i "/File41=/c\ File41='denxyz.${ID}.vtk'" readwf.${ID}.dat
[ ! -f $OUTDIR/vtk/denxyz.${ID}.vtk ] && ./readwf-slurm < readwf.${ID}.dat 2> Err.${ID}.out
[ -f position.${ID}.dat ] && mv position.${ID}.dat $OUTDIR/impurity/
[ -f denxyz.${ID}.vtk ] && mv denxyz.${ID}.vtk $OUTDIR/vtk/
[ -f Err.${ID}.out ] && mv Err.${ID}.out $OUTDIR/readwf_err/
[ -f readwf.${ID}.dat ] && mv readwf.${ID}.dat $OUTDIR/readwf_dat/
r=($(cat $OUTDIR/impurity/position.${ID}.dat))
rr=($(awk '{ print +$1, +$2, +$3 }' <<< $( echo ${r[0]}  ${r[1]}  ${r[2]} )))
ximp=${rr[0]}
yimp=${rr[1]}
zimp=${rr[2]}
cp denxyz.py denxyz.${ID}.py
sed -i "/^denxyz/c\denxyzvtk = LegacyVTKReader(FileNames=['$OUTDIR/vtk/denxyz.${ID}.vtk'])" denxyz.${ID}.py
sed -i "/text1.Text =/c\text1.Text = ' t = ${t} ps'" denxyz.${ID}.py
sed -i "/sphere1.Center = /c\sphere1.Center = [${ximp}, ${yimp}, ${zimp}]" denxyz.${ID}.py
sed -i "/SaveScreenshot/c\SaveScreenshot('denxyz.${ID}.png', magnification=1, quality=100, view=renderView1)" denxyz.${ID}.py
[ ! -f $OUTDIR/xyz3Dplot/denxyz.${ID}.png ] && pvbatch denxyz.${ID}.py 2> /dev/null
[ -f denxyz.${ID}.png ] && mv denxyz.${ID}.png $OUTDIR/xyz3Dplot/
[ -f denxyz.${ID}.py ] && mv denxyz.${ID}.py $OUTDIR/denxyz_py/
