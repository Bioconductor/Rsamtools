#! /bin/sh

repos0="git://github.com/samtools/samtools.git"
repos1="git://github.com/samtools/tabix.git"
wd=`pwd`

dest="../inst/extdata"
if test ! -d $dest; then
	echo "directory does not exist: '$dest'"
	exit 1
fi

dest0="../src/samtools/"
dest1="../src/samtools/bcftools/"
if test ! -d $dest0; then
	echo "directory does not exist: '$dest0'"
	exit 1
fi
if test ! -d $dest1; then
	echo "directory does not exist: '$dest1'"
	exit 1
fi
if test ! -d $dest2; then
	echo "directory does not exist: '$dest2'"
	exit 1
fi
unpack=`mktemp -d -u` ## unsafe: get but don't create temp dir name
git clone $repos0 $unpack

for f in `ls $dest0`; do
	echo "updating file $f"
	cp $unpack/$f $dest0/$f
done
for f in `ls $dest1`; do
	echo "updating file $f"
	cp $unpack/bcftools/$f $dest1/$f
done

dest2="../src/tabix/"
unpack1=`mktemp -d -u` ## unsafe: get but don't create temp dir name
git clone $repos1 $unpack1
for f in `ls $dest2`; do
	echo "updating file $f"
	cp $unpack1/$f $dest2/$f
done
rm -rf $unpack


