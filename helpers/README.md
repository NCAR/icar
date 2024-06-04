# Install Python Dependencies
The following instructions and dependecy files work for the core ICAR scripts.
Tools in `make_domain.py` and ccsm, cesm, cmip, erai, and wrf directories will require the mygis packages as well.
The Python script `ideal_linear.py` will require Nio to be installed with `pip install nio`.

## Setup Environment
### Install With Conda
```bash
$ conda env create -f environment.yml --prefix /path/to/install/icar_env
$ conda activate icar_env

set PYTHONPATH, this will be saved for future use by the environment
$ conda env config vars set PYTHONPATH=$(pwd)/lib:$PYTHONPATH

reactivate environment
$ conda activate icar_env
```

### Install With Pip
```bash
$ pip install -r requirements.txt
```
Make sure the `lib` directory is in the `PYTHONPATH`, add to `.bashrc` or other startup files for repeat use.
```bash
$ export PYTHONPATH=$(pwd)/lib:$PYTHONPATH
```


## Deprecated Scripts
The [Nio package](https://www.pyngl.ucar.edu/Nio.shtml) used in `create_geo_testfiles.py` is not installed as a dependency.
