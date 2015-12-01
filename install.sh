#!/usr/bin/env bash
#
#
# Install 
#

INSTALL_DIR=$HOME/bin

cp -f python_scripts/* $INSTALL_DIR
cp -f bash_scripts/* $INSTALL_DIR
ln -s $PWD/ancient_dna_funcs.sh  $INSTALL_DIR
ln -s $PWD/ancient_pipeline.sh $INSTALL_DIR

