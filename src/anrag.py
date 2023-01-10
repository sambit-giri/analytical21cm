import numpy
import ANRAGpylib

nzl = 120 # this is nzloop used in the fortran code
nkb = 10 # this is nkbin used in fortran code
mbin = 100 # this is massdim_comm
arr = [1.0, 1.0, 1.0, 1.0]	# this is the parameter array 

## Parameter array should be in the frm of fX, Alpha, N_ion/4000, T_vit/10000.
zcut = 15.0 # zcut is the redshift where the sim will stop
z_arr, GLB_arr, k_arr, PSKZ_arr, BSD_R_ARR, BSDZ_arr = ANRAGpylib.wrapping_fn(nzl, nkb, mbin, zcut, arr,npa=len(arr))
# z_arr = execute.wrapping_fn(npa, nzl, nkb, mbin, arr)
# z_arr = execute.wrapping_fn(nzl,nkb,mbin,arr,npa=len(arr))
print(z_arr)     # This output array is z_arr(nzl) carries the redshfits bin
print(GLB_arr)   # This output array is GLB_arr(nzl, 4) which carries the global info at different z: xhii, Tk, Ts, Tb 
print(k_arr)     # This output array is k_arr(nkb) which represents the k-bins
print(PSKZ_arr)  # PSKZ_arr(nzl, nkb) carries power spectrum at different scales at different redshift
print(BSD_R_ARR) # BSD_R_ARR(mbin) is the R-bin at which BSD is estimated
print(BSDZ_arr)  # BSDZ_arr(nzl, mbin) carries BSDs at different mass scales at different redshift







