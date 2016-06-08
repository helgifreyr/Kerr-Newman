# sed -i 's/\\text//g' $1
# sed -i 's/F1/F_1/g' $1
# sed -i 's/F2/F_2/g' $1
# sed -i 's/F0/F_0/g' $1
# sed -i 's/At/A_t/g' $1
# sed -i 's/Aphi/A_\\phi/g' $1
# sed -i 's/phi/\\phi/g' $1
sed -i 's/\^{(1,0)}/_{,r}/g' $1
sed -i 's/\^{(0,1)}/_{,\\theta}/g' $1
sed -i 's/\^{(2,0)}/_{,rr}/g' $1
sed -i 's/\^{(0,2)}/_{,\\theta\\theta}/g' $1
sed -i 's/\^{(1,1)}/_{,r\\theta}/g' $1

