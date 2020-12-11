# Loads python environment on CITA/mussel

source venv-3.6.4/bin/activate
module load python/3.6.4

alias python="python3"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/python/3.6.4/lib
#export PYTHONPATH="${PYTHONPATH}:/mnt/raid-cita/etyhurst/.local/lib/python3.6/site-packages/:/mnt/raid-cita/etyhurst/venv-3.6.4/lib/python3.6/site-packages/"

