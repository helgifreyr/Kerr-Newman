(* ::Package:: *)

Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];
Off[General::spell1]


data=ReadList["functf.dat",Number,RecordLists->True];
nOfFs=7;
radialCoord=ReadList["gridx.dat",Number];
angularCoord=ReadList["gridy.dat",Number];
nx=Length[radialCoord];
ny=Length[angularCoord];


dat=ArrayReshape[data,{nx ny,nOfFs}];


(*1 2 3 4 5 6 7*)(*nr,w,alpha,c1,c2,c3,rh*)
conf=ReadList["res.txt",{Number,Number,Number,Number,Number,Number,Number}];

nr=conf[[1]][[1]];
V0=conf[[1]][[2]];
w=conf[[1]][[3]];
rh=conf[[1]][[4]];

Xtor[X_]:=If[X!=1,Sqrt[(X/(1-X))^2],1000]
f1=Table[{radialCoord[[j]],angularCoord[[i]],dat[[j+(i-1)*nx,1]]},{i,1,ny},{j,1,nx}];
f2=Table[{radialCoord[[j]],angularCoord[[i]],dat[[j+(i-1)*nx,2]]},{i,1,ny},{j,1,nx}];
f0=Table[{radialCoord[[j]],angularCoord[[i]],dat[[j+(i-1)*nx,3]]},{i,1,ny},{j,1,nx}];
Z=Table[{radialCoord[[j]],angularCoord[[i]],dat[[j+(i-1)*nx,4]]},{i,1,ny},{j,1,nx}];
W=Table[{radialCoord[[j]],angularCoord[[i]],dat[[j+(i-1)*nx,5]]},{i,1,ny},{j,1,nx}];
if1=Interpolation[Flatten[f1,1],InterpolationOrder->4];
if2=Interpolation[Flatten[f2,1],InterpolationOrder->4];
if0=Interpolation[Flatten[f0,1],InterpolationOrder->4];
iZ=Interpolation[Flatten[Z,1],InterpolationOrder->4];
iW=Interpolation[Flatten[W,1],InterpolationOrder->4];
f1=if1;
f2=if2;
f0=if0;
Z=iZ;
W=iW;


derVplus[r_,t_]=-((2 W[r, t]g'[r] + (E^(f0[r, t] - f2[r, t])Csc[t] (-2 r H[r] g'[r]+g[r] (r H'[r]+2 H[r](1+r D[f0[r, t],r] - r D[f2[r, t],r]))))/Sqrt[H[r]] - 2 g[r] D[W[r, t],r])/(2 g[r]^2));
derVminus[r_,t_]=(-2 W[r, t] g'[r]+( E^(f0[r, t] - f2[r, t])Csc[t] (-2 r H[r]g'[r]+g[r] (r H'[r]+2 H[r] (1 + r D[f0[r, t],r] - r D[f2[r, t],r]))))/Sqrt[H[r]] + 2 g[r] D[W[r, t],r])/(2 g[r]^2);


g[r_] := r^2+rh^2;
H[r_] := 1 - rh/r;


findRoots[f_,xMax_]:=Block[{zeros,soln,y,x},
zeros=Reap[soln=y[x]/.First[NDSolve[{y'[x]==Evaluate[D[f[x],x]],y[xMax]==(f[xMax])},y[x],{x,xMax,10^-16},Method->{"EventLocator","Event"->y[x],"EventAction":>Sow[{x,y[x]}]}]]][[2,1]];
{Plot[f[x],{x,rh,xMax},Epilog->{PointSize[Medium],Red,Point[zeros]}],zeros[[All,1]],soln}]


derVplusEq[r_]=derVplus[r,Pi/2];
derVminusEq[r_]=derVminus[r,Pi/2];


{dVplusPlot,dVplusRoots,sol1}=findRoots[derVplusEq,10]
{dVminusPlot,dVminusRoots,sol2}=findRoots[derVminusEq,10]



dat=ReadList["fx-inf.txt",Number,RecordLists->True];
lung1=Length[dat] ;


infF0=Table[dat[[i]][[4]],{i,1,lung1}]

constINF=Sum[infF0[[i]],{i,1,ny}]/ny

(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  computation Mass from asymptotics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

Mc=  constINF;
MSch=rh/2 ;
Mass=MSch+Mc;

zMax = r/.NMinimize[{Z[r,Pi/2],r>0},r][[2]]

output={Flatten[{{w,rh,PadRight[dVplusRoots,3],PadRight[dVminusRoots,3],Mass,zMax}},2]};
Print[output]
file = OpenAppend["/home/h/usb/work/kerr-newman/ungauged/BHs/LR-data/rh="<>ToString[rh]<>".dat"]

Export[file, output, "Table"];
WriteString[file, "\n"];

Close[file]

