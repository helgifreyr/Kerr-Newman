from scipy import genfromtxt,linspace
from pylab import *


V=0.3
def plotRH(rh,c):
    inputFile = 'gnuplot-V0='+str(V)+'-rh='+rh+'.dat'
    dat = genfromtxt(inputFile)
    w = dat[:,0]
    g = dat[:,-1]
    plot(w,g,c+'.',ms=2.5)


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
wExtremal = datExtremal[:,0]
gExtremal = datExtremal[:,-1]
plot(wExtremal,gExtremal,'g-',ms=2.5)

xlim(0.64,1)
ylim(0,2)

xlabel(r'$w$',fontsize=18)
ylabel(r'$g$',fontsize=18)
savefig('w-g.pdf')
