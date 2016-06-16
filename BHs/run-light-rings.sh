for dir in m=1.0/V0=0.300/rh=*/1st/w=0.*
do
  echo $dir
  cp light-rings.m $dir
  cd $dir
  math -script light-rings.m > light-ring-output.txt
  rm light-rings.m
  cd ../../../../../
done
for dir in m=1.0/V0=0.300/rh=*/2nd/w=0.*(On)
do
  echo $dir
  cp light-rings.m $dir
  cd $dir
  math -script light-rings.m > light-ring-output.txt
  rm light-rings.m
  cd ../../../../../
done
