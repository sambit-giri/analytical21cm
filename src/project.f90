!!!note that the rank 0 do not do any work.. so use n+1 cores, where n is the number of jobs (if that many cores are available..)

!!!! make
!! mpirun -np 4 ./exe >& exe_out

program analytical_21cm
use param
use subr_main
	implicit none 
	real(dp), dimension(nparam) ::arr
	real(dp), dimension(nzloop) :: z_arr
	real(dp), dimension(nzloop, 4) :: GLB_arr	!!global array xhi, TK, TS, TB
	real(dp), dimension(nkbin) :: k_arr
	real(dp), dimension(nzloop, nkbin) :: PSKZ_arr
	real(dp), dimension(massdim_comm) :: BSD_R_out
	real(dp), dimension(nzloop, massdim_comm) :: BSDZ_arr

	arr=1.d0
	arr(2) = 0.5d0

	call reionization_21cm(arr, z_arr, GLB_arr, k_arr, PSKZ_arr, BSD_R_out, BSDZ_arr)

write(*,*) BSD_R_out
end program analytical_21cm


