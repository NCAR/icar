# FAQ
* [Are you planning on add the Cumulus, PBL or/and land surface schemes to the ICAR?](#Has-ICAR-include-Cumulus,-PBL-and-land-surface-schemes)
* [Can I say that ICAR is an Hydrostatic model or a quasi-geostrophic approximated model?](#Is-ICAR-an-Hydrostatic-model-or-a-quasi-geostrophic-approximated-model)
* [What are the limitations on using the 1st order approximation to solve the advection term?](#What-are-the-limitations-on-using-the-1st-order-approximation-to-solve-the-advection-term?)
* [Is the LUT (Look Up Table) for the dry N^2 or both dry and moist N^2?](#Is-the-LUT-(Look-Up-Table)-for-the-dry-N^2-or-both-dry-and-moist-N^2)
* [How do you include the soundings as an input?](#How-do-you-include-the-soundings-as-an-input?)

## Has ICAR include Cumulus, PBL and land surface schemes?

There is already a cumulus (Tiedtke [LINK]), PBL (greatly simplified form of YSU [LINK]), and LSM (Noah [LINK]) in ICAR. 

## Is ICAR an Hydrostatic model or a quasi-geostrophic approximated model?

Not really. It isn’t a hydrostatic model, or a non-hydrostatic model.  
That assumption doesn’t really come into the modeling framework. (ditto quasi-geostrophic).  
To some extend it inherits those features from whatever model you use to force it though.  

## What are the limitations on using the 1st order approximation to solve the advection term?

Mostly that it ends up creating far too much diffusion, so fields such as water vapor get blurred out a little. 
If you look up the “upwind” scheme, or read one of the MP-DATA papers they will describe this effect.  
Actually most papers that deal with advection will probably talk about this effect.  

## Is the LUT (Look Up Table) for the dry N^2 or both dry and moist N^2?
Both (and neither).  The atmosphere behaves (sort of) the same for a given N^2 whether it is dry or moist. So the LUT itself knows nothing about dry or moist; however, when ICAR runs it calculates the dry or moist N^2 depending on the atmospheric conditions, and uses that to get the correct perturbation from the LUT. 

## How do you include the soundings as an input?
You just have to set up a 3D model grid that is initialized with that sounding profile.  It is not something that ICAR handles directly, it has to be done as a pre-processing step.  However, ICAR has an option to run an “ideal” simulation, in which it assumes that you want the boundary conditions to stay constant over time (e.g. to run a sounding to the steady state solution). 
