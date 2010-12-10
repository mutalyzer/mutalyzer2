#!/bin/sh

skel="../Gen2Phen_Presentation_Skeleton"

if ! [ -e $skel ]; then
  echo Make sure you are in an empty directory and $skel exists.
  exit
fi

if [ -f presentation.tex ]; then
  echo Make sure you are in an empty directory and $skel exists.
  exit
fi

ln -s $skel/bg.eps
ln -s $skel/Makefile
cp $skel/presentation.tex .
