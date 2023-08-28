## Developing
Users are encouraged to modify ICAR for their own specific purposes, and to make those changes available for others by submitting them back to the main repository.  The ICAR code base have been set up to make some additions (e.g. a new physics package) very easy to add, and the developers will do their best to work with anyone who wishes to add more sophisticated changes.

It is worth reading a discussion of the [git workflow](howto/icar_git_workflow.md) used with ICAR.

Then there is a description of [how to work with git and ICAR](howto/icar_and_git_howto.md).

For an outline of the basic code structure see the [ICAR code overview](icar_code_overview.md)

For a more complete documentation of ICAR, see the [detailed ICAR code description](http://ncar.github.io/icar/).  This full description is being updated as the doxygen markup is added to the code.  It also has detailed interactive diagrams of the inter-relationship between all functions.  For the overview, start by looking at the [documentation for the main program](http://ncar.github.io/icar/driver_8f90.html).

ICAR main is now a significant change to the internal workings of ICAR, and has made some components vastly more dynamic and easier to implement/modify.  However, to support that, some elements became more complicated.  For example, adding a new variable to the output file is incredibly easy if that variable already exists.  Just add it to a call (almost anywhere) to `output%add_variables()` however, if that variables has not been added before there are several steps that need to occur.  Below is a quick outline of the process:

1. The variable has to be a a “variable_t” type, not just an array (and preferably part of the domain object though not technically required, but other steps below would be different if it is not.)
2. That variable needs to have an entry in the kVARS structure in the `constants/icar_constants.f90` file.
3. That variable needs to have an entry in the default metadata array `io/default_output_metadata.f90`.
4. That variable needs to be initialized somewhere appropriate (for example in `objects/domain_obj.f90` `create_variables`.)
5. Somewhere you need to tell it to write that variable too.  The easiest place to do this is in the domain_obj routine `var_request`, add it to the list of variables to be allocated, **and** to the list of “restart” variables (i.e. output variables for now).
6. The one change that needs to be made in `io/output_obj.f90` is to add a single line for that variable to the `add_variables` routine.

I know that ended up being more complicated than it could be.  Once a variable is in that format, it is really easy to add and remove variables (and for that matter create multiple output files with different variable lists in each, it is almost magical), but the initial setup is non-trivial to get there.

A few of the steps above should be simplified drastically eventually… it is just complicated by the use of co-arrays, because of that the compiler doesn’t permit a more dynamic specification of these variables (at least the exchangeable variables).  I have a few ways around that, so it is on my list of things to do… it is just much lower priority than, say, making the model work better at the moment. In particular, the add_variables and create_variables steps should be more dynamic and not require this step. The intention of that whole complicated process was to make it easier… but it hasn’t quite had that effect yet.
