#!/bin/bash
#

### Job Name (will be used as prefix later on!)
#PBS -N icar_batch_run
### Project code
#PBS -A P48500028
#PBS -l walltime=00:15:00
#PBS -q regular
### Merge output and error files
#PBS -o job_output/log.out
### job error file (PBS will not allow use of ${PBS_JOBID} here? )
#PBS -e job_output/log.err
### Select X nodes with 36 CPUs each for a total of 72 MPI processes
#PBS -l select=1:ncpus=36:mpiprocs=36:ompthreads=1

### PBS options for automation: https://gif.biotech.iastate.edu/submitting-dependency-jobs-using-pbs-torque

# otherwise xarray is not available:
module load conda/latest
source /glade/u/apps/opt/conda/bin/activate

# echo ${PBS_JOBID::7}

# Set OpenMP variables
export OMP_NUM_THREADS=1
# export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

# the easy way
# icar icar_options.nml

# the complex way (allows a continuous sequence of jobs)
PREFIX=$PBS_JOBNAME

# it is useful to keep all other filenames relative to $PREFIX
# note that this is not required anywhere though
OUTDIR=$PREFIX
OPTFILE=${PREFIX}_options.nml
BATCHFILE=${PREFIX}_batch_submit.sh
TEMPLATE=${PREFIX}_template.nml

# specify the location of the icar executable to use:
EXE=${HOME}/bin/icar

# various useful helper scripts (SETUP_RUN is critical)
SETUP_RUN=${HOME}/icar/helpers/setup_next_run.py
MAKE_TEMPLATE=${HOME}/icar/helpers/make_template.py
MKOUTDIR=mkdir #<user_defined_path>/mkscratch.py # mkscratch creates the directory on scratch and links to it


# --------------------------------------------------
# SHOULD NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
# --------------------------------------------------

# if the template file doesn't exist yet, make it
if [[ ! -e $TEMPLATE ]]; then
    $MAKE_TEMPLATE $OPTFILE $TEMPLATE > job_output/py_mktemp.out
fi

# if the output directory doesn't exist, create it
if [[ ! -e $OUTDIR ]]; then
    $MKOUTDIR $OUTDIR
fi

# if we didn't finish yet we have to continue -BK: but we print this in line 87, so 2 jobs max?
if [[ ! -e ${PREFIX}_finished ]]; then
    # first submit the next job dependant on this one
    qsub -W depend=afterany:${PBS_JOBID} ${BATCHFILE}

    # if we have run before, setup the appropriate restart options
    if [[ -e ${PREFIX}_running ]]; then
        # echo "setting up next run (setup_next_run.py)"
        $SETUP_RUN $OPTFILE $TEMPLATE > job_output/py_setup.out
    fi

    # declare that we have run before so the next job will know
    touch ${PREFIX}_running

    # run the actual executable (e.g. icar options.nml)
    cafrun -n 36 $EXE $OPTFILE >> job_output/icar${PBS_JOBID::7}.out
    # typically the job will get killed while icar is running
    # but for some reason bkilling the job still lets it touch _finished...
    # maybe this will give it a chance to really kill it first?
    sleep 20

    # if icar completes, we are done, tell the next job that we finished
    touch ${PREFIX}_finished
else
    # if the last job ran to completion, delete the inter-job communication files and exit
    rm ${PREFIX}_running
    rm ${PREFIX}_finished
fi
