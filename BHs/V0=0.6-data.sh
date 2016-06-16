rhs=('0.0200' '0.0300' '0.0400' '0.0500' '0.0600' '0.0700')

for rh in $rhs
do
 rm V0=0.6-rh=$rh.dat
 if [[ -d m=1.0/V0=0.600/rh=$rh/1st/ ]]
 then
   echo 'Extracting the 1st branch for', $rh
   for file in m=1.0/V0=0.600/rh=$rh/1st/w=0.*/tmp.txt(On)
   do 
     cat $file >> V0=0.6-rh=$rh.dat
   done
 fi
 if [[ -d m=1.0/V0=0.600/rh=$rh/2nd/ ]]
 then
   echo 'Extracting the 2nd branch for', $rh
   for file in m=1.0/V0=0.600/rh=$rh/2nd/w=0.*/tmp.txt
   do 
     cat $file >> V0=0.6-rh=$rh.dat
   done
 fi
 math -script convert.m "V0=0.6-rh=$rh.dat"
done

rm V0=0.6.dat
for file in V0=0.6-rh=*.dat
do
  cat $file >> V0=0.6.dat
done
math -script convert.m V0=0.6.dat

python w-m.py 0.6
python w-g.py
