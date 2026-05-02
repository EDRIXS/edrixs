!> Global numeric constants and I/O unit numbers used throughout EDRIXS.
!!
!! All floating-point arithmetic uses double precision (dp). Complex
!! constants are constructed with dcmplx so they carry the correct kind.
module m_constants
    implicit none

    integer, public, parameter :: dp = kind(1.0d0)   !< Double-precision kind parameter
    integer, public, parameter :: mystd = 6           !< Standard output unit (stdout)
    integer, public, parameter :: mytmp = 100         !< Scratch file unit number

    real(dp), public, parameter :: zero = 0.0_dp                    !< Real zero
    real(dp), public, parameter :: one  = 1.0_dp                    !< Real one
    real(dp), public, parameter :: pi   = 3.1415926535897932_dp     !< Pi

    complex(dp), public, parameter :: czi   = dcmplx(0.0_dp, 1.0_dp)  !< Imaginary unit i
    complex(dp), public, parameter :: cone  = dcmplx(1.0_dp, 0.0_dp)  !< Complex one
    complex(dp), public, parameter :: czero = dcmplx(0.0_dp, 0.0_dp)  !< Complex zero

    real(dp), public, parameter :: ev2k = 11604.505008098_dp  !< Conversion factor: eV to Kelvin

end module m_constants
