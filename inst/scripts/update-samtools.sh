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

dest="../src/samtools/"
if test ! -d $dest; then
	echo "directory does not exist: '$dest'"
	exit 1
fi
unpack=`mktemp -d -u` ## unsafe: get but don't create temp dir name
svn export -r $rev $repos $unpack
for f in `ls $dest`; do
	echo "updating file $f"
	cp $unpack/$f $dest/$f
done
rm -rf $unpack


