module mod_common2
  use mod_tm
  use mod_eosmodels

  class(EOSModel),  allocatable :: eos_model
  class(TurbModel), allocatable :: turb_model

end module
