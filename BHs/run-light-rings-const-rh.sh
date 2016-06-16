for dir in m=1.0/V0=0.300/rh=0.{0175,0200,0250,0500,0750,1000,1250,1500,1750,2000,2150}/1st/w=0.*(On)
do
  echo $dir
  cp light-rings-extract.m $dir
  cd $dir
  math -script light-rings-extract.m > light-ring-output.txt
  rm light-rings-extract.m
  cd ../../../../../
done
for dir in m=1.0/V0=0.300/rh=0.{0175,0200,0250,0500,0750,1000,1250,1500,1750,2000,2150}/2nd/w=0.*
do
  echo $dir
  cp light-rings-extract.m $dir
  cd $dir
  math -script light-extractnst-rh.m > light-ring-output.txt
  rm light-rings-extract.m
  cd ../../../../../
done
