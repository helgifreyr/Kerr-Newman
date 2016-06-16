from scipy import genfromtxt
from glob import glob
from scipy.interpolate import interp1d
import sys

QeVal = float(sys.argv[1])

files = sorted(glob('gnuplot-w=*dat'))

for file in files:
    dat = genfromtxt(file)
    w = dat[0,0]
    V = dat[0,1]
    M = dat[:,2]
    Qe = -dat[:,14]

    if QeVal<=max(Qe) and QeVal>=min(Qe):
        iM = interp1d(Qe,M)
        print(w,iM(QeVal))
