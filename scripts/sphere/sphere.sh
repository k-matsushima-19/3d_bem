#!/bin/bash
# example: sh sphere.sh id radius filename
# id: 正20面体の各三角形の各辺の分割数
# radius: 半径
# filename: 出力ファイル名 (.off形式)

workdir=`dirname $0`
tmpdir=`mktemp -d`


gfortran -o $tmpdir/kyu $workdir/kyu.f
gfortran -o $tmpdir/merge $workdir/merge.f
gfortran -o $tmpdir/convert $workdir/convert.f90
#id=20
#rad=0.5

id=$1
rad=$2
filename=$3

echo "1 $id $rad 0 0 0" | $tmpdir/kyu > $tmpdir/tmp
#mv kyu_${id}_${rad}.dat tmp
#echo "tmp" | ./merge > kyu_${id}_${rad}.dat
echo "$tmpdir/tmp" | $tmpdir/merge > $tmpdir/tmp2
# rm -f tmp

echo "$tmpdir/tmp2" | $tmpdir/convert

mv $tmpdir/tmp2 $filename
