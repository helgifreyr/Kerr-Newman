rhs=('0.0175' '0.0200' '0.0250' '0.0500' '0.0750' '0.1000' '0.1250' '0.1500' '0.1750' '0.2000' '0.2150')

for rh in $rhs
do
 rm V0=0.3-rh=$rh.dat
 if [[ -d m=1.0/V0=0.300/rh=$rh/1st/ ]]
 then
   echo 'Extracting the 1st branch for', $rh
   for file in m=1.0/V0=0.300/rh=$rh/1st/w=0.*/tmp.txt(On)
   do 
     cat $file >> V0=0.3-rh=$rh.dat
   done
 fi
 if [[ -d m=1.0/V0=0.300/rh=$rh/2nd/ ]]
 then
   echo 'Extracting the 2nd branch for', $rh
   for file in m=1.0/V0=0.300/rh=$rh/2nd/w=0.*/tmp.txt
   do 
     cat $file >> V0=0.3-rh=$rh.dat
   done
 fi
 math -script convert.m "V0=0.3-rh=$rh.dat"
done

rm V0=0.3.dat
for file in V0=0.3-rh=*.dat
do
  cat $file >> V0=0.3.dat
done
math -script convert.m V0=0.3.dat

python w-m.py 0.3
python w-g.py
