# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions


module load python_2.6.4 
module load intel

export PATH=~/shared/usr/bin:~/shared/bin:$PATH
export PS1="[bash-\u@\h \W]\\$ "

source ~/.chemenv.bashrc

export PYTHONPATH=$PYTHONPATH:~/shared/usr/openbabel/lib/python2.6/site-packages:/opt/Python-2.6.4/lib/python2.6/site-packages:~/shared/usr/lib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel/mkl/10.2.5.035/lib/em64t
alias qstat='qstat -q AS'
