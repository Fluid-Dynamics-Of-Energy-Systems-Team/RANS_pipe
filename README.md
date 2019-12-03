# RANS_pipe
RANS code to solve turbulent flows over a flat plate or through a pipe or channel with variable properties

## References 
### Turbulent viscosity models
* [Spalart-Almaras](https://arc.aiaa.org/doi/pdf/10.2514/6.1992-439)
* [V2-f]()
* [k-epsilon](https://www.sciencedirect.com/science/article/pii/0017931072900762)
* [k-omega](https://arc.aiaa.org/doi/abs/10.2514/3.12149)

### Turbulent heat flux models
* [Tang](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205) 
* [Kays](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205) 
* [Kays & Crawford](https://www.sciencedirect.com/science/article/pii/S0017931016300734#b0205)
* [Irrenfried](https://www.sciencedirect.com/science/article/pii/S0142727X17304083?via%3Dihub)
* [2 equations](https://www.infona.pl/resource/bwmeta1.element.baztech-article-BAT3-0010-0041)

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
