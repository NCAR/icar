#!/bin/bash
### Job Name (will be used as prefix later on!)
#SBATCH --job-name="ICAR_tst"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=00:05:00
#SBATCH --constraint=haswell
#SBATCH --qos=debug
### Project code
#SBATCH --account=m4062
### error and output files in separate folder, name with jobid (%x) an job name (%j)
### N.B: create the job_output folder before submitting this job!
#SBATCH --output=job_output/log-%x.%j.out
#SBATCH --error=job_output/log-%x.%j.err

# Make sure a python environment with xarray is available:
module load python
conda activate myenv

# Set OpenMP variables
export OMP_NUM_THREADS=1
# export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS

# the easy way
# icar icar_options.nml

# the complex way (allows a continuous sequence of jobs)
PREFIX=tst  ##$SBATCH_JOB_NAME

# it is useful to keep all other filenames relative to $PREFIX
# note that this is not required anywhere though
OUTDIR=$PREFIX
OPTFILE=options.nml  #${PREFIX}_options.nml
BATCHFILE=batch_submit_SLURM.sh #${PREFIX}_batch_submit.sh
TEMPLATE=${PREFIX}_template.nml

# the ICAR executable to use
EXE=$HOME/bin/icar_dbs

# load any environmental settings to run icar properly (system dependent):
. /global/cfs/cdirs/m4062/env_scripts/UO-GNU-env.sh


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

# # if the output directory doesn't exist, create it
# if [[ ! -e $OUTDIR ]]; then
#     $MKOUTDIR $OUTDIR
# fi

# if we didn't finish yet we have to continue -BK: but we print this in line 87, so 2 jobs max?
if [[ ! -e ${PREFIX}_finished ]]; then
    # first submit the next job dependant on this one
    # sub -w "ended(${PBS_JOBID})" < $BATCHFILE
    # qsub -W depend=afterany:${PBS_JOBID} ${BATCHFILE}  ## PBS version
    sbatch --dependency=afternotok:$SLURM_JOB_ID ${BATCHFILE}

    # if we have run before, setup the appropriate restart options
    if [[ -e ${PREFIX}_running ]]; then
        # echo "setting up next run (setup_next_run.py)"
        $SETUP_RUN $OPTFILE $TEMPLATE > job_output/py_setup.out
    fi

    # declare that we have run before so the next job will know
    touch ${PREFIX}_running

    # run the actual executable (e.g. icar options.nml)
    # cafrun -n 36 $EXE $OPTFILE > job_output/icar_$SLURM_JOB_ID.out
    cafrun -n 36 $EXE $OPTFILE >> job_output/icar.out ### if you prefer one log file for the icar output

    # typically the job will get killed while icar is running
    # but for some reason bkilling the job still lets it touch _finished...
    # maybe this will give it a chance to really kill it first?
    sleep 10

    # if icar completes, we are done, tell the next job that we finished
    # BK dont understand this: wont it prevent the next (or after-next job from starting (ln 63))
    touch ${PREFIX}_finished
else
    # if the last job ran to completion, delete the inter-job communication files and exit
    rm ${PREFIX}_running
    rm ${PREFIX}_finished
fi
