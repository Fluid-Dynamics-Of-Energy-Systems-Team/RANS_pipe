module mod_common2
  use mod_turbmodels
  use mod_eosmodels

  class(EOSModel),  allocatable :: eos_model
  class(TurbModel), allocatable :: turb_model

end module
