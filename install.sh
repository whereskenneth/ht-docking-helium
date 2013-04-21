#!/bin/bash

cp ~/.bashrc ~/.bashrc.backup
echo "Backing up .bashrc"

cp env/.bashrc ~
cp env/.chemenv.bashrc ~


cd src/dlgtosdf

echo "Compiling dlgtosdf"
make
ln -s dlgtosdf ../../bin

cd ../sdftopdbqt

echo "Compiling sdftopdbqt"
make
ln -s sdftopdbqt ../../bin

cd ../../

echo "Finished"




