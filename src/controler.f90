!!!note that the rank 0 do not do any work.. so use n+1 cores, where n is the number of jobs (if that many cores are available..)

!!!! make
!! mpirun -np 4 ./exe >& exe_out
!!global array xhi, TK, TS, TB


SUBROUTINE wrapping_fn(npa, nzl, nkb, mbin, zcut, arr, z_arr, GLB_arr, k_arr, PSKZ_arr, BSDZRarr, BSDZ_arr)
use param
use subr_main
implicit none
integer, intent(in) :: npa, nzl, nkb, mbin
real(8), intent(in) :: zcut
real(8), dimension(npa), intent(in) ::arr
real(8), dimension(nzl), intent(out) :: z_arr
real(8), dimension(nzl, 4), intent(out) :: GLB_arr
real(8), dimension(nkb), intent(out) :: k_arr
real(8), dimension(nzl, nkb), intent(out) :: PSKZ_arr
real(8), dimension(mbin), intent(out) :: BSDZRarr
real(8), dimension(nzl, mbin), intent(out) :: BSDZ_arr

call reionization_21cm(arr, zcut, z_arr, GLB_arr, k_arr, PSKZ_arr, BSDZRarr, BSDZ_arr)

end SUBROUTINE wrapping_fn


