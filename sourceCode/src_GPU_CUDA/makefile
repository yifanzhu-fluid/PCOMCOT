all:
	rm -f pcomcot-gpu	
	nvfortran -g -Mpreprocess -cudalib=cublas,cusparse VariableDefination.cuf okada.cuf pcomcotLIB_CUDA.cuf \
	solveSWEs.cuf dispersion.cuf breaker.cuf \
	pcomcotNetCDFlibEmpty.cuf pcomcotLIB.cuf pcomcot.cuf -o pcomcot-gpu
	rm -f *.mod
double:
	rm -f pcomcot-gpu	
	nvfortran -g -Mpreprocess -DDOUBLE -cudalib=cublas,cusparse VariableDefination.cuf okada.cuf pcomcotLIB_CUDA.cuf \
	solveSWEs.cuf dispersion.cuf breaker.cuf \
	pcomcotNetCDFlibEmpty.cuf pcomcotLIB.cuf pcomcot.cuf -o pcomcot-gpu
	rm -f *.mod
