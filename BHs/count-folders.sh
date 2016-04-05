touch temp.txt
for dir in "$@"
do
  echo $dir >> temp.txt
done

wc temp.txt
rm temp.txt

# 100.0: 2597
# 350.0: 2078
# 750.0: 1859
