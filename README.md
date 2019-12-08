# RANS_pipe
RANS code to solve turbulent flows over a flat plate or through a pipe or channel with variable properties

## References 
### Turbulent viscosity models
* [Spalart-Almaras](https://arc.aiaa.org/doi/pdf/10.2514/6.1992-439)
* [V2-f]()
* [k-epsilon](https://www.sciencedirect.com/science/article/pii/0017931072900762)
* [k-omega](https://arc.aiaa.org/doi/abs/10.2514/3.12149)

### Turbulent heat flux models
* [Tang et al. (2016)](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205) 
* [Kays](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205) 
* [Kays & Crawford](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205)
* [Irrenfried](https://www.sciencedirect.com/science/article/pii/S0142727X17304083?via%3Dihub)
* [Bae](https://www.infona.pl/resource/bwmeta1.element.elsevier-e6af3d8b-9871-32d3-a9f6-4972a82f5f76)
* [Deng](https://www.sciencedirect.com/science/article/abs/pii/S0017931000001319)
* [Nagano & Kim](https://asmedigitalcollection.asme.org/heattransfer/article-abstract/110/3/583/382763/A-Two-Equation-Model-for-Heat-Transport-in-Wall?redirectedFrom=fulltext)
## To Do
* Implement two equation models for the turbulent heat flux
* Write out residuals
* Fix the non-equidistant grid for the turbulence models


## Results
<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/clean/notebooks/bl.png"
     style="float: center; margin-right: 10px;" />
     
<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/clean/notebooks/channel.png"
     style="float: center; margin-right: 10px;" />
     
<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/clean/notebooks/pipe.png"
     style="float: center; margin-right: 10px;" />

## References
<a id="1">[1]</a> 
Tang, G. et al (2016).
A variable turbulent Prandtl number model for simulating supercritical pressure CO2 heat transfer.
International Journal of Heat and Mass Transfer, 102, 1082 - 1092

<a id="2">[2]</a> 
Deng, B. et al (2001).
A near-wall two-equation heat transfer model for wall turbulent flows.
International Journal of Heat and Mass Transfer, 44(4), 691 - 698

<a id="3">[3]</a> 
Irrenfried, C., Steiner, H. (2017).
DNS based analytical P-function model for RANS with heat transfer at high Prandtl numbers.
International Journal of Heat and Fluid Flow, 66, 217-225

<a id="4">[4]</a> 
Jones, W.P., Launder, B.E. (1972).
The prediction of laminarization with a two-equation model of turbulence.
International Journal of Heat and Mass Transfer, 15(2), 301-314

<a id="5">[5]</a> 
Menter, F.R. (1994) 
Two-equation eddy-viscosity turbulence models for engineering applications. 
AIAA journal, 32(8), 1598-1605.

<a id="5">[5]</a> 
Spalart, P., Allmaras, S. (1992).
A one-equation turbulence model for aerodynamic flows.
30th aerospace sciences meeting and exhibit. 

<a id="6">[6]</a>
Nagano, Y., Kim, C. (1988).
A two-equation model for heat transport in wall turbulent shear flows. 
Journal Heat Transfer, 110(3), 583-589.

<a id="7">[7]</a>
Kays, W.M., Crawford, M.E., Weigand, B. (2005). 
Convective heat and mass transfer. 
Boston: McGraw-Hill Higher Education, 76,81-89. 

