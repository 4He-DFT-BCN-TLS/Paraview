#!/bin/bash
source /opt/intel/bin/compilervars.sh intel64

# Call: ./draw-slurm-chdb.sh density.xxxx.dat path/to/he-wfs/ output-directory


declare -a X
declare -a Y
declare -a Z

dt=0.5
ID=${1##*density.}   # density.0001.dat => 0001.dat
ID=${ID%.dat}      # 0001.dat         => 0001
printf -v t "%.1f" $(echo "${ID}*${dt}"|bc) # store current time in $t
INDIR=${2}
OUTDIR=${3}

(cd ${OUTDIR}; mkdir denxyz_py impurity readwf_dat readwf_err vtk xyz3Dplot 2> /dev/null )

cp readwf.dat readwf.${ID}.dat
gsed -i "/File10=/c\ File10='${INDIR}/density.${ID}.dat'" readwf.${ID}.dat
gsed -i "/File31=/c\ File31='position.${ID}.dat'" readwf.${ID}.dat
gsed -i "/File41=/c\ File41='denxyz.${ID}.vtk'" readwf.${ID}.dat
[ ! -f $OUTDIR/vtk/denxyz.${ID}.vtk ] && ./readwf-slurm < readwf.${ID}.dat 2> Err.${ID}.out
[ -f position.${ID}.dat ] && mv position.${ID}.dat $OUTDIR/impurity/
[ -f denxyz.${ID}.vtk ] && mv denxyz.${ID}.vtk $OUTDIR/vtk/
[ -f Err.${ID}.out ] && mv Err.${ID}.out $OUTDIR/readwf_err/
[ -f readwf.${ID}.dat ] && mv readwf.${ID}.dat $OUTDIR/readwf_dat/

coords=($(cat $OUTDIR/impurity/position.${ID}.dat))
num_imp=$(cat nimp.dat)
for ((i=0;i<num_imp;i++))
do
	Xi=$((i*3+0))
	Yi=$((i*3+1))
	Zi=$((i*3+2))
	X[$i]=${coords[$Xi]}
	Y[$i]=${coords[$Yi]}
	Z[$i]=${coords[$Zi]}
done

cp denxyz-MI-top.py top.${ID}.py
cp denxyz-MI-bottom.py bottom.${ID}.py
gsed -i "/^denxyz/c\denxyzvtk = LegacyVTKReader(FileNames=['$OUTDIR/vtk/denxyz.${ID}.vtk'])" top.${ID}.py
gsed -i "/text1.Text =/c\text1.Text = ' t = ${t} ps'" top.${ID}.py
### HERE STARTS THE PART FOR THE IMPURITY SPHERES
for ((i=0;i<num_imp;i++))
do
	CPO=$((i+1))
	gsed "s/sphere1/sphere${CPO}/g" denxyz-MI-middle.py >> middle.${ID}.py
	gsed -i "/sphere${CPO}.Center/c\sphere${CPO}.Center = [${X[$i]}, ${Y[$i]}, ${Z[$i]}]" middle.${ID}.py
done
gsed -i "/SaveScreenshot/c\SaveScreenshot('denxyz.${ID}.png', renderView1, ImageResolution=[2833, 1800])" bottom.${ID}.py
cat top.${ID}.py middle.${ID}.py bottom.${ID}.py > denxyz.${ID}.py
rm top.${ID}.py middle.${ID}.py bottom.${ID}.py nimp.dat

[ ! -f $OUTDIR/xyz3Dplot/denxyz.${ID}.png ] && pvbatch denxyz.${ID}.py 2> /dev/null
[ -f denxyz.${ID}.png ] && mv denxyz.${ID}.png $OUTDIR/xyz3Dplot/
[ -f denxyz.${ID}.py ] && mv denxyz.${ID}.py $OUTDIR/denxyz_py/
