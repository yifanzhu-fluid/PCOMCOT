all:
	rm -f pcomcot
	mpifort -g -I/usr/include VariableDefination.f90 okada.f pcomcotNetCDFlib.f90 pcomcotLIB.f90 MPICommunicationLIB.f90 BiCGStabLIB.f90 nonlinearFlux.f90 dispersion.f90 breaker.f90 solvEQNs.f90 pcomcot.f90 -L/usr/lib/x86_64-linux-gnu -lnetcdff -o pcomcot
	rm -f VariableDefination.mod variabledefination.mod
nocdf:
	rm -f pcomcot
	mpifort -g VariableDefination.f90 okada.f pcomcotNetCDFlibEmpty.f90 pcomcotLIB.f90 MPICommunicationLIB.f90 BiCGStabLIB.f90 nonlinearFlux.f90 dispersion.f90 breaker.f90 solvEQNs.f90 pcomcot.f90 -o pcomcot
	rm -f VariableDefination.mod variabledefination.mod
checkbounds:
	rm -f pcomcot
	mpifort -g -fcheck=bounds -I/usr/include VariableDefination.f90 okada.f pcomcotNetCDFlib.f90 pcomcotLIB.f90 MPICommunicationLIB.f90 BiCGStabLIB.f90 nonlinearFlux.f90 dispersion.f90 breaker.f90 solvEQNs.f90 pcomcot.f90 -L/usr/lib/x86_64-linux-gnu -lnetcdff -o pcomcot
	rm -f VariableDefination.mod variabledefination.mod
