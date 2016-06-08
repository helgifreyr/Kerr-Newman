(* ::Package:: *)

Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];
Off[General::spell1]
<<VariationalMethods`



(* 1  2   3   4  5  6   7*)
(* nr,w,alpha,c1,c2,c3,rh*)
conf=ReadList["res.txt",{Number,Number ,Number ,Number ,Number,Number,Number }];

 
nr=conf[[1]][[1]];
rh= conf[[1]][[4]];
V= conf[[1]][[2]];
w= conf[[1]][[3]];
gr=ReadList["gridx.dat",{Number}];


lgr=Length[gr];
nx=lgr;

listar=Table[gr[[k]][[1]],{k,1,lgr}] ;
listalogr=Table[Log[10,gr[[k]][[1]]],{k,1,lgr}];



dat=ReadList["f-tPi2.txt",Number,RecordLists->True];
lung1=Length[dat] ;

 nF1=Table[dat[[i]][[2]],{i,1,lung1}];
 nF2=Table[dat[[i]][[3]],{i,1,lung1}];
nF0=Table[dat[[i]][[4]],{i,1,lung1}];
nZ=Table[dat[[i]][[5]],{i,1,lung1}];
nW=Table[dat[[i]][[6]],{i,1,lung1}];

order=3;
Subscript[F, 1]=ListInterpolation[nF1,{listar},InterpolationOrder->order];
Subscript[F, 2]=ListInterpolation[nF2,{listar},InterpolationOrder->order];
Subscript[F, 0]=ListInterpolation[nF0,{listar},InterpolationOrder->order];
Z=ListInterpolation[nZ,{listar},InterpolationOrder->order];
W=ListInterpolation[nW,{listar},InterpolationOrder->order];

(* Here I compute first derivatives *)
 
dat=ReadList["fx-tPi2.txt",Number,RecordLists->True];
lung1=Length[dat] ;

 nF1x=Table[dat[[i]][[2]],{i,1,lung1}];
 nF2x=Table[dat[[i]][[3]],{i,1,lung1}];
nF0x=Table[dat[[i]][[4]],{i,1,lung1}];
nZx=Table[dat[[i]][[5]],{i,1,lung1}];
nWx=Table[dat[[i]][[6]],{i,1,lung1}];

order=3;
Subscript[F, 1]'=ListInterpolation[1/(1+listar)^2 nF1x,{listar},InterpolationOrder->order];
Subscript[F, 2]'=ListInterpolation[1/(1+listar)^2 nF2x,{listar},InterpolationOrder->order];
Subscript[F, 0]'=ListInterpolation[1/(1+listar)^2 nF0x,{listar},InterpolationOrder->order];
Z'=ListInterpolation[1/(1+listar)^2 nZx,{listar},InterpolationOrder->order];
W'=ListInterpolation[1/(1+listar)^2 nWx,{listar},InterpolationOrder->order];

derVplus[r_] = (-2 W[r] Derivative[1][g][r] + 2 g[r] Derivative[1][W][r] + (E^(Subscript[F, 0][r] -Subscript[F, 2][r]) (g[r] Derivative[1][H][r] -2 H[r] (Derivative[1][g][r] +g[r] (-Derivative[1][Subscript[F, 0]][r] +Derivative[1][Subscript[F, 2]][r]))))/Sqrt[H[r]])/(2 g[r]^2)
derVminus[r_] = (-2 W[r] Derivative[1][g][r] + 2 g[r] Derivative[1][W][r] + (E^(Subscript[F, 0][r] -Subscript[F, 2][r]) (-g[r] Derivative[1][H][r] +2 H[r] (Derivative[1][g][r] +g[r] (-Derivative[1][Subscript[F, 0]][r] +Derivative[1][Subscript[F, 2]][r]))))/Sqrt[H[r]])/(2 g[r]^2)


g[r_] := r^2+rh^2;
H[r_] := 1 - rh/r;
gtt[r_]:=-E^(2 Subscript[F,0][r]) r^2/g[r] H[r]+E^(2 Subscript[F,2][r]) r^2 g[r](W[r]/g[r])^2
rh


findRoots[f_,xMax_]:=Block[{zeros,soln,y,x},
  zeros=Reap[soln=y[x]/.First[NDSolve[{y'[x]==Evaluate[D[f[x],x]],y[xMax]==(f[xMax])},y[x],{x,xMax,10^-20},Method->{"EventLocator","Event"->y[x],"EventAction":>Sow[{x,y[x]}]}]]][[2,1]];
{Plot[f[x],{x,rh,xMax},Epilog->{PointSize[Medium],Red,Point[zeros]}],zeros[[All,1]],soln}]


{dVplusPlot,dVplusRoots,sol1}=findRoots[derVplus,10]
{dVminusPlot,dVminusRoots,sol2}=findRoots[derVminus,10]
{ergoPlot,ergoRoots,sol3}=findRoots[gtt,20]

dat=ReadList["fx-inf.txt",Number,RecordLists->True];
lung1=Length[dat] ;
unghi0=ReadList["gridy.dat",{Number}];
ny=Length[unghi0];


infF0=Table[dat[[i]][[4]],{i,1,lung1}]

constINF=Sum[infF0[[i]],{i,1,ny}]/ny

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  computation Mass from asymptotics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

Mc=  constINF;

MSch=rh/2 ;
Mass=MSch+Mc;

zMax = r/.NMinimize[{Z[r],r>0},r][[2]]

output={Flatten[{{w,rh,PadRight[ergoRoots,3],PadRight[dVplusRoots,3],PadRight[dVminusRoots,3],Mass,zMax}},2]};
Print[output]
file = OpenAppend["/home/h/usb/work/kerr-newman/ungauged/BHs/LR-ER-data/rh="<>ToString[rh]<>".dat"]

Export[file, output, "Table"];
WriteString[file, "\n"];

Close[file]
