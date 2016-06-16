for dir in "$@"; do
  if [ ! -f $dir/tmp.txt ]; # comment out this line to run on already run folders
  then # and this
    echo $dir
    cp extract-global-data-ungauged-BHs.m $dir/
    cd $dir
    rm tmp.txt
    math -script extract-global-data-ungauged-BHs.m > math-output.txt
    rm extract-global-data-ungauged-BHs.m
    cd ../../../../../
  fi # and this
done
