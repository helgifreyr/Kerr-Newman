rsync -avzsh blafis:kerr-newman/ungauged/extremal/m=1.0 .
./re-run.sh m=1.0/V0=0.*/*/*/*
grep Int m=1.0/V0=0.*/*/*/*/tmp.txt
