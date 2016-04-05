from scipy import genfromtxt,linspace
from pylab import *

datZM = genfromtxt('zero-mode-V0=0.3.dat')
V = 0.3
rh = datZM[:,-3]
a = datZM[:,-2]
w = a/(a*a+rh*rh)
Qe = V * (a*a+rh*rh)/rh
M = (a*a+rh*rh+Qe*Qe)/(2*rh)
plot(w,M,'b--',ms=2.5)
extremalLine = lambda w: sqrt((1.0 - 4.0*V*V + sqrt(1.0+8.0*V*V))/w)/(2.0*sqrt(2.0*w))
ws = linspace(0.64,1,100)
plot(ws,extremalLine(ws),'k-',ms=2.5)


def plotRH(rh,c):
    inputFile = 'gnuplot-V0='+str(V)+'-rh='+rh+'.dat'
    dat = genfromtxt(inputFile)
    w = dat[:,0]
    M = dat[:,2]
    plot(w,M,c+'.',ms=2.5)


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

datBS = genfromtxt('BS.txt')
wBS = datBS[:,0]
MBS = datBS[:,1]
plot(wBS,MBS,'r-',ms=2.5)

datExtremal = genfromtxt('gnuplot-extremal-V0=0.3.dat')
wExtremal = datExtremal[:,0]
MExtremal = datExtremal[:,2]
plot(wExtremal,MExtremal,'g-',ms=2.5)

xlim(0.64,1)
ylim(0,1.4)

xlabel(r'$w$')
ylabel(r'$M$')
savefig('w-M.pdf')
