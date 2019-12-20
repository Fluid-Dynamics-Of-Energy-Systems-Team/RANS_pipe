# Canonical RANS Solver
RANS code to solve turbulent flows over a flat plate or through a pipe and channel with variable properties as a function of temperature only ([low mach number approximation](https://ccse.lbl.gov/Research/LowMach/lowMach.html)). Periodic boundary conditions for the inlet and outlet will produce fully developed profiles that can be used as an inlet boundary condition for a developing flow with heating (using either an isoflux or isothermal boundary condition).
## Mesh
The code supports a non-equidistant mesh in wall-normal direction as well as streamwise direction with refinement where the boundary conditions change.

<img src="https://github.com/Fluid-Dynamics-Of-Energy-Systems-Team/RANS_pipe/raw/blstart/notebooks/non-equidistant.svg?sanitize=true"
     style="float: center; margin-right: 10px;" />
## Models
### Turbulent viscosity models
The reynolds shear stress is modeled using:
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\overline{\rho{u_i''}{u_i''}}=-{\mu_{t}}\left(2\frac{\partial{\tilde{u_i}}}{{\partial}{x_i}}-\frac{2}{3}{\nabla}{\cdot}\tilde{u}\right)+\frac{2}{3}\overline{\rho}\tilde{k}" />
</p>

To model the turbulent viscosity <img src="https://latex.codecogs.com/svg.latex?\Large&space;\mu_t" />, the following turbulent viscosity models can be used:
* Spalart & Allmaras (1992) [[6]](#6)
* Durbin (1995) [[9]](#9)
* Jones & Launder (1972) [[4]](#4) 
* Myong & Kasagi (1990)  [[5]](#5)
* Abe (1994)             [[14]](#14)

To account for the varying properties, corrections for these models can be used as presented by Catris & Aupoix (2000) [[12]](#12) and Otero et al. (2018) [[13]](#13).

### Turbulent heat flux models
The turbulent heat flux is modeled according to:

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\overline{\rho{u_j''}{h''}}=-{\rho}{\alpha_t}\frac{\partial\tilde{h}}{\partial{x_j}},"/>
</p>

where the turbulent diffusivity <img src="https://latex.codecogs.com/svg.latex?\Large&space;\alpha_t" /> can be calculated using the zero equations models, given by:
<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;\alpha_t=\frac{\mu_t}{Pr_t}," />
</p>

where <img src="https://latex.codecogs.com/svg.latex?\Large&space;{\mu}_t" /> is the turbulent viscosity and <img src="https://latex.codecogs.com/svg.latex?\Large&space;{Pr}_{t}" /> the turbulent Prandtl number. The turbulent Prandtl number can be modeled using the following models:

* Tang et al. (2016) [[1]](#1)
* Kays (2005) [[11]](#11)
* Kays & Crawford (2005) [[8]](#8)
* Irrenfried & Steiner (2017) [[3]](#3)
* Bae (2016) [[10]](#10)

To calculate the thermal diffusivity using 2 extra transport equations, the following models are included:
* Deng et al. (2001) [[2]](#2) 
* Nagano & Kim (1988) [[7]](#5)

## To Do
* fix the dissipation bc for the Deng model and the NK model
* Convert the v2f model to a non-equidistant grid
* Write out residuals
* Restart with the turbulent diffusivity variables interpolated
* Add the convective outlet bc for the scalars

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
[Durbin, P.A. (1995)
Separated flow computations with the k-epsilon-v-squared model.
AIAA journal, 33(4), 659-664.](https://link.springer.com/article/10.1007/s10494-005-1974-8)

<a id="10">[10]</a>
[Bae, Y.Y. (2016).
A new formulation of variable turbulent Prandtl number for heat transfer to supercritical fluids. 
International Journal of Heat and Mass Transfer, 92, 792-806.](https://www.infona.pl/resource/bwmeta1.element.elsevier-e6af3d8b-9871-32d3-a9f6-4972a82f5f76)

<a id="11">[11]</a>
[Kays, W.M. (1994).
Turbulent Prandtl Number—Where Are We?
Journal of Heat Transfer, 116(2), 284-295.](https://asmedigitalcollection.asme.org/heattransfer/article/116/2/284/383190/Turbulent-Prandtl-Number-Where-Are-We)

<a id="12">[12]</a>
[Catris, S., Aupoix, B. (2000).
Density corrections for turbulence models.
Aerospace Science and Technology, 4.1,1-11.](https://www.sciencedirect.com/science/article/pii/S1270963800001127)


<a id="13">[13]</a>
[Otero R., G.J., Patel, A., Diez S., R., Pecnik, R. (2018).
Turbulence modelling for flows with strong variations in thermo-physical properties.
International Journal of Heat and Fluid Flow, 73, 114-123.](https://www.sciencedirect.com/science/article/pii/S0142727X18301978)

<a id="14">[14]</a>
[Abe, K., Kondoh, T., Nagano, Y. (1994).
A new turbulence model for predicting fluid flow and heat transfer in separating and reattaching flows—I. Flow field calculations.
International journal of heat and mass transfer, 37(1), 139-151.](https://www.sciencedirect.com/science/article/pii/0017931094901686)
