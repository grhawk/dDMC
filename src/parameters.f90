MODULE parameter
  USE precision
  IMPLICIT NONE
  
  ! Conversion Coefficients form Atomic Units to Common Value
  real(kr),parameter :: BohrAngst     = 0.529177249      ! 1 bohr    --> <value> Angstrom
  real(kr),parameter :: HartKcalMol   = 627.509          ! 1 hartree --> <value> Kcal/Mol

  ! Pure Number
  real(kr),parameter :: avogadro = 6.02214199E-23   ! Avogadro Number

END MODULE parameter
