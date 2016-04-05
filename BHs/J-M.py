from scipy import genfromtxt,linspace
from pylab import *

V=0.3
def plotRH(rh,c):
    inputFile = 'gnuplot-V0='+str(V)+'-rh='+rh+'.dat'
    dat = genfromtxt(inputFile)
    J = dat[:,4]
    M = dat[:,2]
    plot(J,M,c+'.',ms=2.5)


plotRH('0.0200','b')
plotRH('0.0250','b')
plotRH('0.0500','b')
plotRH('0.0750','b')
plotRH('0.1000','b')
plotRH('0.1250','b')
plotRH('0.1500','b')
plotRH('0.1750','b')
plotRH('0.2000','b')
plotRH('0.2150','b')

datExtremal = genfromtxt('gnuplot-extremal-V0=0.3.dat')
JExtremal = datExtremal[:,4]
MExtremal = datExtremal[:,2]
plot(JExtremal,MExtremal,'g-',ms=2.5)

datSolitons = genfromtxt('BS.txt')
JSolitons = datSolitons[:,2]
MSolitons = datSolitons[:,1]
plot(JSolitons,MSolitons,'r-',ms=2.5)

datZM = genfromtxt('zero-mode-V0=0.3.dat')
rh = datZM[:,-3]
a = datZM[:,-2]            
w = a/(a*a+rh*rh)
Qe = V * (a*a+rh*rh)/rh
M = (a*a+rh*rh+Qe*Qe)/(2*rh)
J=a*M
plot(J,M,'b-',ms=2.5)

# xlim(0.64,1)
# ylim(0,2)

xlabel(r'$J$')
ylabel(r'$M$')
savefig('J-M.pdf')
