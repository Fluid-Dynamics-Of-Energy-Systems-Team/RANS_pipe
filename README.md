# Canonical RANS Solver
RANS code to solve turbulent flows over a flat plate or through a pipe and channel with variable properties as a function of temperature only ([low mach number approximation](https://ccse.lbl.gov/Research/LowMach/lowMach.html)). Periodic boundary conditions for the inlet and outlet will produce fully developed profiles that can be used as an inlet boundary condition for a developing flow with heating (using either an isoflux or isothermal boundary condition).
## Mesh
The code supports a non-equidistant mesh in wall-normal direction as well as streamwise direction with refinement where the boundary conditions change.

<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/blstart/notebooks/non-equidistant.svg?sanitize=true"
     style="float: center; margin-right: 10px;" />
## Models
### Turbulent viscosity models
To model the reynolds stress, the following turbulent viscosity models are included:
* Spalart & Allmaras (1992) [[6]](#6)
* Laurence et al. (2004) [[9]](#9)
* Jones & Launder (1972) [[4]](#4) 
* Myong & Kasagi (1990)  [[5]](#5)

### Turbulent heat flux models
To model the turbulent heat flux, the following turbulent diffusivity models are included:
* Tang et al. (2016) [[1]](#1)
* Kays (2005) [[11]](#11)
* Kays & Crawford (2005) [[8]](#8)
* Irrenfried & Steiner (2017) [[3]](#3)
* Bae (2016) [[10]](#10)
* Deng et al. (2001) [[2]](#2) 
* Nagano & Kim (1988) [[7]](#5)

## To Do
* Find the bug in the k-epsilon turbulence model
* Convert the v2f model to a non-equidistant grid
* Write out residuals
* Write out the turbulent prandtl number/the two variable associated with the turbulent diffusivity

## Results
<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/clean/notebooks/bl.png"
     style="float: center; margin-right: 10px;" />
     
<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/clean/notebooks/channel.png"
     style="float: center; margin-right: 10px;" />
     
<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/clean/notebooks/pipe.png"
     style="float: center; margin-right: 10px;" />

## References
<a id="1">[1]</a> 
[Tang, G., Shi, H., Wu, Y., Lu, J., Li, Z., Liu, Q., Zhang, H. (2016).
A variable turbulent Prandtl number model for simulating supercritical pressure CO2 heat transfer.
International Journal of Heat and Mass Transfer, 102, 1082 - 1092](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205)

<a id="2">[2]</a> 
[Deng, B. et al (2001).
A near-wall two-equation heat transfer model for wall turbulent flows.
International Journal of Heat and Mass Transfer, 44(4), 691 - 698](https://www.sciencedirect.com/science/article/abs/pii/S0017931000001319)

<a id="3">[3]</a> 
[Irrenfried, C., Steiner, H. (2017).
DNS based analytical P-function model for RANS with heat transfer at high Prandtl numbers.
International Journal of Heat and Fluid Flow, 66, 217-225](https://www.sciencedirect.com/science/article/pii/S0142727X17304083?via%3Dihub)

<a id="4">[4]</a> 
[Jones, W.P., Launder, B.E. (1972).
The prediction of laminarization with a two-equation model of turbulence.
International Journal of Heat and Mass Transfer, 15(2), 301-314](https://www.sciencedirect.com/science/article/pii/0017931072900762)

<a id="5">[5]</a> 
[Myong, H.K.,Kasagi ,N. (1990)
A new approach to the improvement of k-epsilon turbulence model for wall-bounded shear flows. 
JSME International Journal, 33, 63–72](https://arc.aiaa.org/doi/abs/10.2514/3.12149)

<a id="6">[6]</a> 
[Spalart, P., Allmaras, S. (1992).
A one-equation turbulence model for aerodynamic flows.
30th aerospace sciences meeting and exhibit.](https://arc.aiaa.org/doi/pdf/10.2514/6.1992-439)

<a id="7">[7]</a>
[Nagano, Y., Kim, C. (1988).
A two-equation model for heat transport in wall turbulent shear flows. 
Journal Heat Transfer, 110(3), 583-589.](https://asmedigitalcollection.asme.org/heattransfer/article-abstract/110/3/583/382763/A-Two-Equation-Model-for-Heat-Transport-in-Wall?redirectedFrom=fulltext)

<a id="8">[8]</a>
[Kays, W.M., Crawford, M.E., Weigand, B. (2005). 
Convective heat and mass transfer. 
Boston: McGraw-Hill Higher Education, 76, 81-89.](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205)

<a id="9">[9]</a>
[Laurence, D.R., Uribe J.C., Utyuzhnikov, S.V. (2004). 
A Robust Formulation of the v2-f Model. 
Flow Turbulence and Combustion, 73, 169-185](https://link.springer.com/article/10.1007/s10494-005-1974-8)

<a id="10">[10]</a>
[Bae, Y.Y. (2016).
A new formulation of variable turbulent Prandtl number for heat transfer to supercritical fluids. 
International Journal of Heat and Mass Transfer, 92, 792-806.](https://www.infona.pl/resource/bwmeta1.element.elsevier-e6af3d8b-9871-32d3-a9f6-4972a82f5f76)

<a id="11">[11]</a>
[Kays, W.M. (1994).
Turbulent Prandtl Number—Where Are We?
Journal of Heat Transfer, 116(2), 284-295.](https://asmedigitalcollection.asme.org/heattransfer/article/116/2/284/383190/Turbulent-Prandtl-Number-Where-Are-We)
