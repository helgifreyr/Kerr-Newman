from scipy import genfromtxt,linspace
from pylab import *
from glob import glob
import sys

V = float(sys.argv[1])
datZM = genfromtxt('zero-mode-V0='+str(V)+'.dat')
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

# if V==0.3:
#     rhs = [0.0200,0.0250,0.0500,0.0750,0.1000,0.1250,0.1500,0.1750,0.2000,0.2150]
#     for rh in rhs:
#         plotRH('%5.4f'%rh,'b')
# if V==0.6:
#     rhs = [0.0200,0.0300,0.0400,0.0500,0.0600,0.0700]
#     for rh in rhs:
#         plotRH('%5.4f'%rh,'b')

datBS = genfromtxt('BS.txt')
wBS = datBS[:,0]
MBS = datBS[:,1]
plot(wBS,MBS,'r-',ms=2.5)

datExtremal = genfromtxt('gnuplot-extremal-V0='+str(V)+'.dat')
wExtremal = datExtremal[:,0]
MExtremal = datExtremal[:,2]
plot(wExtremal,MExtremal,'g-',ms=2.5)

for file in sorted(glob('Qe=*.dat')):
    datQe = genfromtxt(file)
    wQe = datQe[:,0]
    MQe = datQe[:,1]
    plot(wQe,MQe,'m--',ms=1.5)

xlim(0.64,1)
ylim(0,1.4)

xlabel(r'$w$')
ylabel(r'$M$')
savefig('w-M-V0='+str(V)+'-const-Qe.pdf')
