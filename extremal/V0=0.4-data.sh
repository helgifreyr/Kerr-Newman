for file in m=1.0/V0=0.4000/rh=0.1000/*/w=0.*/tmp.txt; do cat $file; done > V0=0.4-rh=0.1.dat
cat V0=0.4-rh=0.1.dat > extremal-V0=0.4.dat
cp extremal-V0=0.4.dat ../BHs/
