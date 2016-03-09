-------------------------------------------------------------------------
Readme document on changes made to the model when implementing YSU-scheme
-------------------------------------------------------------------------

The YSU-scheme needs additional input variables and extra variables for initialization 
that did not exist in this model before. Normal procedure is to get relevant input from the 
surface layer scheme, however a surface layer scheme is not yet implemented such that needed 
variables are computed and or initialized in time_step.f90 and pbl_driver.f90.
The YSU scheme is a non-local closure scheme based on Richardson-Bulk # as stability measure 
to calculate the respective integrated stability functions from the similarity theory.

The stability classes as it is now are:
1. Rib >= 0.2 (stable/nighttime conditions)
2. 0 < Rib < 0.2 (damped, mechanical turbulence regime)
3. Rib = 0 (forced convection)
4. Rib < 0 free convection

Based on theses classes the integrated stability functions psim and psih are computed. The 
computations of Rib, psim and psih are done in time_step.f90. Rib, psim and psih are pure 
ingoing variables and are not being changed by YSU-scheme.

An important variable is the pbl-height which has to be initialized and then grows to a height 
dependant on the physical state of the pbl (stability). As I understand now the pbl-height is 
being altered in the right manner by the YSU-scheme but has to be initialized thus this variable 
is of type inout.

Other variables that need to be calculated are hol, ust, hpbl, znt, wspd, zol. These are meant 
to be inout. However, I do not see in the code changes made to ust which would be given back 
out of the function. The same is true for zol, znt, wspd (although this is used to calculate 
the wspd at 10m and there is also a one-dimensional change to wspd1 which is only a column. 
As far as I understand this does not feedback on the full 2D wspd field). I consider therefore 
ust, zol, znt, wspd variables of type in rather than inout.

hol is altered in the YSU-code also because of the change of pbl-height. The
pbl-height is initialized in pbl_driver.f90 as a constant and not yet using
the equation (1) recommended in Hong et al. 2006.

In order to obtain zol and hol not only the similarity ust is needed but also the temperature 
scale thstar. thstar again depends on the integrated similarity stability function psih described 
above. the same is true for ust which of course needs to be altered by psim. To my understanding 
these updates need to be done outside of YSU-scheme (e.g. in the surface layer scheme) and then 
fed back to the YSU-scheme.

Further the virtual potential temperature is needed "thv" as well as the potential
temperature near the ground "thvg". Both of which are calculated in time_step.f90.
After Hing et al. 2006 thvg is the sum of thv and thT (= virtual temperature
excess near the ground) and initialized without thT as thv. For the next step
thT should be taken into account.
In Jiminez et al. 2012 they do not elaborate on this but they describe thvg as the virtual 
potential temperature near the ground which I calculate in time_step.f90 using the 
exner-function to transform t2m into th2m and then add the humidity contribution resulting 
from use of a moist gas.

A primary goal of the surface layer parameterization is the calculation of the dimenensionless 
bulk transfer coefficients for momentum, heat and moisture. There is a parameter called exch_h 
which is going into the YSU-scheme which I assume to be the transfer coefficient for heat. 
For the YSU-scheme the transfer coefficient for heat is needed as a variable of type inout. 
It seems also to be altered in the YSU-code. However, since this is usually a parameter calculated 
by the surface scheme I am not yet sure whether exch_h should only be initialized and then left 
alone or whether I should update it along with other variables. Currently exch_h is initialized 
as a 2D field of constants exch_h = 0.1. 

u10 and v10 are variables of type in and should be updated using the stability fcts rather than 
the log-law. This is done in the surface layer scheme. Currently the calculations are still 
done using the log law.

---------------- !!! Changes to be done next !!! --------------------------------
- Separate initialization and next time steps as suggested in Hing et al. 2006
  and Jiminez et al. 2012
- Include the update and calculation of the dimenensionless bulk transfer coefficients for heat
  (momentum and moisture secondary since they do not seem to be needed by YSU). !!! DONE !!!
- Include update and calculation of ustar from similarity function (already
  existing in time_step.f90 but not yet used) !!! DONE !!!
- Include update of thvg using equation from Hong et al. 2006 (a bit difficult
  since I do not have the virtual heat flux from the surface), might not be
  crucial if I use a different approach e.g. using th near the ground.
