all:
	rm -f pcomcot-cpu
	mpifort -g -I/usr/include VariableDefination.f90 okada.f90 pcomcotNetCDFlib.f90 pcomcotLIB.f90 MPICommunicationLIB.f90 \
	BiCGStabLIB.f90 dispersion.f90 breaker.f90 solveSWEs.f90 pcomcot.f90 -L/usr/lib/x86_64-linux-gnu -lnetcdff -o pcomcot-cpu
	rm -f VariableDefination.mod variabledefination.mod
nocdf:
	rm -f pcomcot-cpu
	mpifort -g VariableDefination.f90 okada.f90 pcomcotNetCDFlibEmpty.f90 pcomcotLIB.f90 MPICommunicationLIB.f90 \
	BiCGStabLIB.f90 dispersion.f90 breaker.f90 solveSWEs.f90 pcomcot.f90 -o pcomcot-cpu
	rm -f VariableDefination.mod variabledefination.mod
