export CHEMAXON_HOME=~/shared/opt/ChemAxon
export JCHEM_HOME=$CHEMAXON_HOME/JChem
export JCHEM_PATH=$JCHEM_HOME/bin

###############################
#Amber Variables

export AMBERHOME=~/shared/opt/amber11
export AMBER_PATH=$AMBERHOME/bin

################################
#Openeye variables

export OE_LICENSE=~/shared/opt/openeye/oe_license.txt
export OE_PATH=~/shared/opt/openeye/bin
export OE_ARCH=redhat-RHEL5-x64:centos-4.2-x86_64:centos-4.2-g++3.4-x86_64:redhat-RHEL6-x64


#################################
#Openbabel variables

export LD_LIBRARY_PATH=~/shared/usr/lib:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=~/shared/usr/include:$C_INCLUDE_PATH
export CPP_INCLUDE_PATH=~/shared/usr/include:$CPP_INCLUDE_PATH


#################################
#Dovis variables

export DOVIS=~/shared/usr/dovis
export MGL_ROOT=$DOVIS/MGLTools
export MGL_ARCHOSV=$MGL_ROOT/i86Linux2
#export DOVIS_PATH=$DOVIS/bin:$DOVIS/scripts:$MGL_ROOT/i86Linux2/bin:$MGL_ROOT/MGLToolsPckgs/AutoDockTools/Utilities24:$PATH
#export PYTHONHOME=$MGL_ROOT/share:$MGL_ARCHOSV
#export PYTHONPATH=$MGL_ROOT/MGLToolsPckgs:~/usr/local/lib/python2.6

#################################
#Dock6 variables

export DOCK6_PATH=~/shared/opt/dock6/bin

#################################
#Final load/export

export PATH=$DOCK6_PATH:$OE_PATH:$AMBER_PATH:$JCHEM_PATH:$PATH
