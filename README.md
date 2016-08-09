# SLCTM: SimpleLightCurveTransitModel
Currently, the SLCTM.LightCurve (Python 2.7) class encapsulates the orbital and transit parameters for a single planet. The class also provides storage for transit times, errors, fluxes and flux timestamps. In addition to the transit data, transit start/stop indices are saved. By importing the SLCTM module, you also gain access to a few functions that operate on a SLCTM.LightCurve object. This Readme gives some basic information on the class, with usage tips. Provided along with this documentation file are a few examples of how to use the class and the module's functions. You can find the project's github here.

The batman module (Bad-Ass Transit Model Calculation), developed by Laura Kreidberg is used this within the SLCTM module to create an ensemble of different models. The SLCTM module provides some custom wrapper functions for batman, but the LightCurve class is abstract enough to be extended to take data - simulated or otherwise - from any source.

The emcee module developed by Dan Foreman-Mackey is used to regress a function for the light curve



File List:

SLCTM.py -> Contains SLCTM Class and several useful functions that apply to it.

ex_plotdata.py -> This example showcases the production of a light curve inside the SLCTM.LightCurve class.

ex_chisqs.py -> This example showcases the production of a simulated light curve and a large number of model light curves
varying one parameter, which then have their ChiSquareds evaluated and plotted.

ex_fit.py -> This example shows the usage of emcee with SLCTM.


#Acknowledgments

The SLCTM Module was started by Joseph Sheldon Jr as part of a summer project supervised by Prof. Eric B. Ford. Prof. Ford provided the initial idea and lots of advice/mentoring.

Credits to Laura Kreidberg for batman and to Dan Foreman-Mackey for emcee
