# f2py -c jenkinslib.f90 -m jenkinslib
gfortran -c lib_array.f90 -o lib_array.o -fPIC
gfortran -c jenkinslibtools.f90 -o jenkinslibtools.o -fPIC
# f2py -c jenkinslib.f90 lib_array.o -m jenkinslib
#f2py -c jenkinslib.f90 jenkinslibtools.o lib_array.o -m jenkinslib
f2py -m jenkinslib -c jenkinslib.f90 jenkinslibtools.o lib_array.o -L/usr/lib -llapack
