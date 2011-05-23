#! /bin/sh

repos="https://samtools.svn.sourceforge.net/svnroot/samtools/trunk/samtools/"
wd=`pwd`

dest="../inst/extdata"
if test ! -d $dest; then
	echo "directory does not exist: '$dest'"
	exit 1
fi

svn info $repos > $dest/samtools-svninfo.txt
rev=`sed -n 's/Revision: //p' < $dest/samtools-svninfo.txt`

dest0="../src/samtools/"
dest1="../src/samtools/bcftools/"
dest2="../src/tabix/"
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
svn export -r $rev $repos $unpack
for f in `ls $dest0`; do
	echo "updating file $f"
	cp $unpack/$f $dest0/$f
done
for f in `ls $dest1`; do
	echo "updating file $f"
	cp $unpack/bcftools/$f $dest1/$f
done
for f in `ls $dest2`; do
	echo "updating file $f"
	cp $unpack/tabix/$f $dest2/$f
done
rm -rf $unpack


