#!/bin/bash
#
### Script assumes 'output-directory' has been created before execution
### CALL:	./draw-slurm-chdb.sh density.xxxx.dat path/to/he-wfs/ output-directory
### EXMPL:	./draw-slurm-chdb.sh density.005.dat ../he-wfs test

#source /opt/intel/bin/compilervars.sh intel64 ## Only uncomment on a Mac
SED=sed			## 'sed' on Linux, 'gsed' on a Mac with Homebrew

USAGE="Usage:		./draw-slurm-chdb.sh <density.NNN.dat> <path/to/he-wfs/> <output-directory/>"
EXAMPLE="Example:	./draw-slurm-chdb.sh density.005.dat ../he-wfs/ test/"
[ "$1" = "" ] && echo $USAGE && echo  $EXAMPLE && exit 1
[ "$2" = "" ] && echo $USAGE && exit 1
[ "$3" = "" ] && echo $USAGE && exit 1



DT=0.5
ID=${1##*density.}	# density.0001.dat => 0001.dat
ID=${ID%.dat}		# 0001.dat         => 0001
INDIR=${2}
OUTDIR=${3}
HEIGHT=1080			# Rescale the output image to this height
printf -v t "%.1f" $(echo "${ID}*${DT}"|bc) # store current time in $t
declare -a X
declare -a Y
declare -a Z

(cd ${OUTDIR}; mkdir denxyz_py impurity nimp readwf_dat readwf_err vtk xyz3Dplot 2> /dev/null)

cp readwf.dat readwf.${ID}.dat
$SED -i "/File10=/c\ File10='${INDIR}/density.${ID}.dat'" readwf.${ID}.dat
$SED -i "/File31=/c\ File31='position.${ID}.dat'" readwf.${ID}.dat
$SED -i "/File41=/c\ File41='denxyz.${ID}.vtk'" readwf.${ID}.dat
$SED -i "/nimpFile=/c\ nimpFile='nimp.${ID}.dat'" readwf.${ID}.dat
[ ! -f $OUTDIR/vtk/denxyz.${ID}.vtk ] && ./readwf-slurm < readwf.${ID}.dat 2> Err.${ID}.out
[ -f position.${ID}.dat ] && mv position.${ID}.dat $OUTDIR/impurity/
[ -f nimp.${ID}.dat ] && mv nimp.${ID}.dat $OUTDIR/nimp/
[ -f denxyz.${ID}.vtk ] && mv denxyz.${ID}.vtk $OUTDIR/vtk/
[ -f Err.${ID}.out ] && mv Err.${ID}.out $OUTDIR/readwf_err/
[ -f readwf.${ID}.dat ] && mv readwf.${ID}.dat $OUTDIR/readwf_dat/

coords=($(cat $OUTDIR/impurity/position.${ID}.dat))
num_imp=$(cat $OUTDIR/nimp/nimp.${ID}.dat)
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
$SED -i "/^denxyzvtk =/c\denxyzvtk = LegacyVTKReader(FileNames=['$OUTDIR/vtk/denxyz.${ID}.vtk'])" top.${ID}.py
$SED -i "/text1.Text =/c\text1.Text = ' t = ${t} ps'" top.${ID}.py

### HERE STARTS THE PART FOR THE IMPURITY SPHERES
for ((i=0;i<num_imp;i++))
do
	CPO=$((i+1))
	$SED "s/sphere1/sphere${CPO}/g" denxyz-MI-middle.py >> middle.${ID}.py
	$SED -i "/sphere${CPO}.Center/c\sphere${CPO}.Center = [${X[$i]}, ${Y[$i]}, ${Z[$i]}]" middle.${ID}.py
done
$SED -i "/SaveScreenshot/c\SaveScreenshot('denxyz.${ID}.png', renderView1, ImageResolution=[2833, 1800])" bottom.${ID}.py
cat top.${ID}.py middle.${ID}.py bottom.${ID}.py > denxyz.${ID}.py
rm top.${ID}.py middle.${ID}.py bottom.${ID}.py
### HERE ENDS THE PART FOR THE IMPURITY SPHERES

[ ! -f $OUTDIR/xyz3Dplot/denxyz.${ID}.png ] && pvbatch --no-mpi denxyz.${ID}.py 2> /dev/null
convert denxyz.${ID}.png -resize x${HEIGHT} denxyz.${ID}.png
WIDTH=$(convert denxyz.${ID}.png -print "%w" /dev/null)
REM=$((WIDTH%2))
[ $REM -ne 0 ] && convert denxyz.${ID}.png -chop 1x0 denxyz.${ID}.png
[ -f denxyz.${ID}.png ] && mv denxyz.${ID}.png $OUTDIR/xyz3Dplot/
[ -f denxyz.${ID}.py ] && mv denxyz.${ID}.py $OUTDIR/denxyz_py/
