setup_local_env
setup_root_afs_gcc46 

basedir=$HOME/local

export BATINSTALLDIR=$basedir

export PATH=$PATH:${basedir}/bin:$BATINSTALLDIR
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${basedir}/lib

