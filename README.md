
this is the fjordmelt tool after the plume model (jenkins,2011) and Cowton
see publication: beckmann, 2018, thecryosphere
and a tool to create a Temperature profile as in the
glacier model. Wiritten in Fortran but wrapped for python,
it  can be used in python for direct investigations.
when dowloaded in to (e.g. "/home/user/gitrepos/fjord/"")
compile with ./build.sh
then can be used interactivley in ipython:
sys.path.insert(0, "/home/user/gitrepos/fjord/")
import jenkinslib as jenk
m, y = jjj.jenkinslib.melting_process(X,Z,ta,sa,E0,sina,y0,len(X),'cone')
with:
X (1d array) - path under glacier tongue, X=-Z for veritcal cliff
Z (1d array)  - depth, Z is negaitve
Ta (1d array) - Ta (Z) temepertur profile from bottom to water surface
Sa (1d array) - Sa (Z) salinity profile from bottom to water surface
E0 - Enrtainement coefficeint (0.025-0.16)
sina (1d array) - sina (Z) angele under tongue, sina = 1 for vertical cliff
y0  - array, length 4 initial starting values, for velocity, thickness, salinity
temperture ,times flux  see code
'cone' or '1d' for type of plume
output:
m (1d array)- melt rate along the glacier tongue
y (4, 1d array) - plume proerties along the glacier tongue
( for velocity, thickness, salinity
temperture ,times flux  see code)
