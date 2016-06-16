from scipy import genfromtxt,linspace
from pylab import *

def extremalPlot(V,c):
    datZM = genfromtxt('zero-mode-V0='+str(V)+'.dat')
    rh = datZM[:,-3]
    a = datZM[:,-2]
    w = a/(a*a+rh*rh)
    Qe = V * (a*a+rh*rh)/rh
    M = (a*a+rh*rh+Qe*Qe)/(2*rh)
    plot(w,M,c+'-',ms=2.5)
    extremalLine = lambda w: sqrt((1.0 - 4.0*V*V + sqrt(1.0+8.0*V*V))/w)/(2.0*sqrt(2.0*w))
    ws = linspace(0.75,1,100)
    plot(ws,extremalLine(ws),c+'-',ms=2.5)

    datExtremal = genfromtxt('gnuplot-extremal-V0='+str(V)+'.dat')
    wExtremal = datExtremal[:,0]
    MExtremal = datExtremal[:,2]
    plot(wExtremal,MExtremal,c+'--',ms=2.5)

extremalPlot(0.3,'b')
extremalPlot(0.4,'g')
extremalPlot(0.5,'m')
extremalPlot(0.6,'y')

V = 0
ws = linspace(0.75,1,100)
extremalLine = lambda w: sqrt((1.0 - 4.0*V*V + sqrt(1.0+8.0*V*V))/w)/(2.0*sqrt(2.0*w))
plot(ws,extremalLine(ws),'r-',ms=2.5)
datExtKerr = genfromtxt('eHBH-kerr.dat')
wExtKerr = datExtKerr[:,0]
MExtKerr = datExtKerr[:,1]
plot(wExtKerr,MExtKerr,'r--',ms=2.5)

datZMKerr = genfromtxt('ZM-kerr.dat')
wZMKerr = datZMKerr[:,5]
MZMKerr = datZMKerr[:,1]
plot(wZMKerr,MZMKerr,'r-',ms=2.5)

xlim(0.75,1)
ylim(0,1.3)

xlabel(r'$w$')
ylabel(r'$M$')
savefig('comparison.pdf')
