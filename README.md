# literature

To help validate the pyrolysis models, comparisons were made with results from literature. Each folder (listed below) in this repository contains code examples which utilize formulas or theories proposed in literature. These examples were compared to results or figures in the paper and later applied to other models. References are also listed in the comments of each code file.

## Blasi2003
Correlation for devolatilization time as tv where tv = time for 95% conversion. Based on Equation 2 and Figure 12 in the paper.
- Di Blasi, Colomba, and Carmen Branca. "Temperatures of wood particles in a hot sand bed fluidized by 
nitrogen." Energy & fuels 17.1 (2003): 247-254.

## Diego2003
Correlation for devolatilization time that considers moisture content within a wood particle. Refer to Equations 2-7 in the paper for more details.
- De Diego, L. F., et al. "Effect of moisture content on devolatilization times of pine wood particles in a 
fluidized bed." Energy & fuels 17.2 (2003): 285-290.

## Churchill1974
Correlations from section "Functions Which Cross One Limiting Solution". Plots should be same as in Figures 9 and 10 on page 41.
- Churchill, Stuart W., and R. Usagi. "A standardized procedure for the production of correlations in the form of a common empirical equation." Industrial & Engineering Chemistry Fundamentals 13.1 (1974): 39-44.

## Ganser1993
Drag coefficient, Cd, for non-spherical particles. Calculate terminal velocity, ut, from drag coefficient, Cd, and particle sphericity. See page 150, Equation 18 and Table 7.
- Ganser, Gary H. "A rational approach to drag prediction of spherical and nonspherical particles." Powder Technology 77.2 (1993): 143-152.

## Keshavarz2006
An improved lumped capacitance method for Biot <= 20. Traditional method was typically restricted to Bi < 0.1.
- Keshavarz, P., and M. Taheri. "An improved lumped analysis for transient heat conduction by using the polynomial approximation method." Heat and mass transfer 43.11 (2007): 1151-1156.

## Kunii1991
Calculate terminal velocity, ut, from ut* and dp* for different particle sphericities. See pages 80-83, Equations 28-29, and Equations 31-33 in the book. Applicable for sphericity from 0.5 to 1.0.
- Kunii, D., and O. Levenspiel. "Fluidization Engineering, 2nd Ed.", 1991.

## Ranzi2008
Plots reaction rate constant, K, and chemical species as function of temperature, T from cellulose, hemicellulose, and lignin components. Note that Table 3 lists 7 reactions but 8 reaction rates are given.
- Ranzi, Eliseo, et al. "Chemical kinetics of biomass pyrolysis." Energy & Fuels 22.6 (2008): 4292-4300.

## Ranzi2013
Plots chemical species as a function of time for cellulose, hemicellulose, and lignin biomass components. Lignin is divided into three groups as oxygen (lig-O), carbon (lig-C), and hydrogen (lig-H) rich components. Reaction rate constants are also plotted.
- Ranzi, Eliseo, et al. "Kinetic modeling of the thermal degradation and combustion of biomass." Chemical Engineering Science 110 (2014): 2-12.

## Sadhukhan2009
1-D transient heat conduction model compared to experimental values for a cylindrical and spherical wood shapes. Experimental data taken from Figures 1 (cylinder) and 2 (sphere).
- Sadhukhan, Anup Kumar, Parthapratim Gupta, and Ranajit Kumar Saha. "Modelling of pyrolysis of large wood particles." Bioresource technology 100.12 (2009): 3134-3139.

## Santos2010
Horio and Nonaka bubble diameter correlation using Equation 14.30 on page 322. Bed expansion facter, fbexp, and expanded bed height, zexp, for a fluidized bed using Equations 14.7 and 14.18 from pages 318-320. Calculate terminal velocity, ut, for a near spherical particle using Equations 4.12-4.15 on page 86.
- Marcio L. de Souza-Santos, "Solid Fuels Combustion and Gasification: Modeling, Simulation, and Equipment Operations", 2nd Ed., 2010.

## License
Code in this repository is available under the MIT license. See the LICENSE file for more info.
