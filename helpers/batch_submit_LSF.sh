#!/bin/bash
#
# LSF batch script to run ICAR
#
#BSUB -P P48500028           # project code
#BSUB -W 12:00               # wall-clock time (hrs:mins)
#BSUB -n 1                   # number of tasks in job
#BSUB -R "span[ptile=16]"    # run 16 MPI tasks per node
#BSUB -J            run_name            # job name
#BSUB -o job_output/run_name.%J.out     # job output file (%J is replaced by the job ID)
#BSUB -e job_output/run_name.%J.err     # job error file (%J is replaced by the job ID)
#BSUB -q regular             # queue

# Set OpenMP variables
export OMP_NUM_THREADS=16
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

# the easy way
# icar icar_options.nml

# the complex way (allows a continuous sequence of jobs)
PREFIX=run_name

# it is useful to keep all other filenames relative to $PREFIX
# note that this is not required anywhere though
OUTDIR=$PREFIX
OPTFILE=${PREFIX}_options.nml
BATCHFILE=${PREFIX}_batch_submit.sh
TEMPLATE=${PREFIX}_template.nml

EXE=/glade/u/home/gutmann/bin/icar

# various useful helper scripts (SETUP_RUN is critical)
SETUP_RUN=/glade/u/home/gutmann/src/icar/helpers/setup_next_run.py
MAKE_TEMPLATE=/glade/u/home/gutmann/src/icar/helpers/make_template.py
MKOUTDIR=/glade/u/home/gutmann/bin/python/mkscratch.py # mkscratch creates the directory on scratch and links to it

# --------------------------------------------------
# SHOULD NOT NEED TO MODIFY ANYTHING BELOW THIS LINE
# --------------------------------------------------

# if the template file doesn't exist yet, make it
if [[ ! -e $TEMPLATE ]]; then
    $MAKE_TEMPLATE $OPTFILE $TEMPLATE
fi

# if the output directory doesn't exist, create it
if [[ ! -e $OUTDIR ]]; then
    $MKOUTDIR $OUTDIR
fi

# if we didn't finish yet we have to continue
if [[ ! -e ${PREFIX}_finished ]]; then
    # first submit the next job dependant on this one
    bsub -w "ended(${LSB_JOBID})" < $BATCHFILE

    # if we have run before, setup the appropriate restart options
    if [[ -e ${PREFIX}_running ]]; then
        $SETUP_RUN $OPTFILE $TEMPLATE
    fi

    # declare that we have run before so the next job will know
    touch ${PREFIX}_running

    # run the actual executable (e.g. icar options.nml)
    $EXE $OPTFILE
    # typically the job will get killed while icar is running
    # but for some reason bkilling the job still lets it touch _finished...
    # maybe this will give it a chance to really kill it first?
    sleep 10

    # if icar completes, we are done, tell the next job that we finished
    touch ${PREFIX}_finished
else
    # if the last job ran to completion, delete the inter-job communication files and exit
    rm ${PREFIX}_running
    rm ${PREFIX}_finished
fi
