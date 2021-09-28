#! /bin/sh

for i in {1..25}
do
    bsub -q l -o  simplerejection/x1Bkg/logfile_overlay_"$i".txt  -e simplerejection/x1Bkg/errorfile_overlay_"$i".txt "root -l -b -q create_overlay.cc"
done