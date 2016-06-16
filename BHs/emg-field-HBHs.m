(* ::Package:: *)

Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];
Off[General::spell1];


coords={t,r,\[Theta],\[Phi]};


(* ::Text:: *)
(*Just as with the analytical case, we define the metric.*)


 gdd={{g00,0,0,g03},{0,g11,0,0},{0,0,g22,0},{g03,0,0,g33}};
MatrixForm[gdd]
g00=-E^(2f0[r,\[Theta]]) r^2/g[r] H[r]+E^(2f2[r,\[Theta]]) Sin[\[Theta]]^2 g[r](W[r,\[Theta]]/g[r])^2;
g03=-E^(2f2[r,\[Theta]]) Sin[\[Theta]]^2 W[r,\[Theta]];
g11=E^(2f1[r,\[Theta]])/H[r];
g22=E^(2f1[r,\[Theta]]) g[r];
g33=E^(2f2[r,\[Theta]]) g[r]Sin[\[Theta]]^2;MatrixForm[gdd]


(* ::Text:: *)
(*Its inverse and determinant*)


guu=Inverse[gdd];
dim=Length[gdd];


Do[g[i,j]=Simplify[gdd[[i,j]]];,{i,1,dim},{j,1,dim}];
Do[gu[i,j]=Simplify[guu[[i,j]]];,{i,1,dim},{j,1,dim}];


DET=Simplify[Det[gdd]];


(* ::Text:: *)
(*The tetrad is the same in terms of the metric elements*)


euu={{e00,0,0,e03},{0,e11,0,0},{0,0,e22,0},{e30,0,0,e33}};


(*we=0;*)(*Static observer*)
we=-g03/g33;(*Zero angular momentum observer*)
(*we=a/(r^2+a^2);*)(*Carter observer*)
(*we=Sqrt[M]/(r^(3/2)+Sqrt[M]a)*)(*Freely orbiting stationary observer*)
ut=1/Sqrt[-g00-2 g03 we-g33 we^2];
e00=ut;
e03=ut we;
e11=Sqrt[gu[2,2]];
e22=Sqrt[gu[3,3]];
e30=-Sqrt[ut^2+gu[1,1]];
e33=Sqrt[ut^2 we^2+gu[4,4]];
edd=Simplify[Inverse[euu],Csc[\[Theta]]>0&&r>0&&f2[r,\[Theta]]>0&&f1[r,\[Theta]]>0&&f0[r,\[Theta]]>0&&H[r]>0];
dim=Length[gdd];


Do[e[i,j]=edd[[i,j]];,{i,1,dim},{j,1,dim}];


Do[eu[i,j]=euu[[i,j]];,{i,1,dim},{j,1,dim}];


MatrixForm[euu]


(* ::Text:: *)
(*And we find the electromagnetic tensor w.r.t the tetrad in the same way*)


A = {At[r,\[Theta]]-Aphi[r,\[Theta]] Sin[\[Theta]] W[r,\[Theta]]/g[r],0,0,Aphi[r,\[Theta]]Sin[\[Theta]]};
MatrixForm[A]


(* F_\[Mu]\[Nu] = d_\[Mu] A_\[Nu] - d_\[Nu] A_\[Mu] *)
F = Table[D[A[[\[Nu]]],coords[[\[Mu]]]] - D[A[[\[Mu]]],coords[[\[Nu]]]],{\[Mu],1,dim},{\[Nu],1,dim}];
MatrixForm[F];


(* F_ab = F_\[Mu]\[Nu] e_a^\[Mu] e_b^\[Nu] *)
tetradF = Simplify[Table[Sum[F[[\[Mu],\[Nu]]]euu[[b,\[Mu]]]euu[[c,\[Nu]]],{\[Mu],1,dim},{\[Nu],1,dim}],{b,1,dim},{c,1,dim}],Sin[\[Theta]]>0];
MatrixForm[tetradF];


\[CurlyEpsilon]=LeviCivitaTensor[4];


b[r_,\[Theta]_]=Table[Sum[\[CurlyEpsilon][[1,i,j,k]] tetradF[[j,k]],{j,1,dim},{k,1,dim}],{i,1,dim-1}];


e[r_,\[Theta]_]=Table[tetradF[[1,i+1]],{i,1,dim-1}];


H[r_]:=1-rh/r;
g[r_] := r^2+rh^2;


(* ::Text:: *)
(*Here I read the solution and interpolate it onto a (r,\[Theta]) grid.*)


dat=ReadList["funct.dat",Number,RecordLists->True];
nx=251;ny=30;
xmax=dat[[251,1]];xmin=dat[[1,1]];
ymax=dat[[nx*ny,2]];ymin=dat[[1,2]];


conf=ReadList["res.txt",{Number,Number,Number,Number,Number,Number,Number}];
nr=conf[[1]][[1]];
rh=conf[[1]][[4]];
V=conf[[1]][[2]];
w=conf[[1]][[3]];


Xtor[X_]:=If[X!=1,Sqrt[(X/(1-X))^2+rh^2],1000]


f1=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,3]]},{i,1,ny},{j,1,nx}];
f2=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,4]]},{i,1,ny},{j,1,nx}];
f0=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,5]]},{i,1,ny},{j,1,nx}];
Z=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,6]]},{i,1,ny},{j,1,nx}];
W=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,7]]},{i,1,ny},{j,1,nx}];
Aphi=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,8]]},{i,1,ny},{j,1,nx}];
At=Table[{Xtor[dat[[j+(i-1)*nx,1]]],dat[[j+(i-1)*nx,2]],dat[[j+(i-1)*nx,9]]},{i,1,ny},{j,1,nx}];


if1=Interpolation[Flatten[f1,1],InterpolationOrder->4];
if2=Interpolation[Flatten[f2,1],InterpolationOrder->4];
if0=Interpolation[Flatten[f0,1],InterpolationOrder->4];
iZ=Interpolation[Flatten[Z,1],InterpolationOrder->4];
iW=Interpolation[Flatten[W,1],InterpolationOrder->4];
iAphi=Interpolation[Flatten[Aphi,1],InterpolationOrder->4];
iAt=Interpolation[Flatten[At,1],InterpolationOrder->4];


f1=if1;
f2=if2;
f0=if0;
Z=iZ;
W=iW;
Aphi=iAphi;
At=iAt;


(* ::Text:: *)
(*Just as with the analytical case, we define some transformations from the tetrad coordinates to Cartesian coordinates*)


r=.;\[Theta]=.;\[Phi]=.;x=.;y=.;z=.;
x[r_,\[Theta]_,\[Phi]_]:=r Sin[\[Theta]] Cos[\[Phi]]
y[r_,\[Theta]_,\[Phi]_]:=r Sin[\[Theta]] Sin[\[Phi]]
z[r_,\[Theta]_,\[Phi]_]:=r Cos[\[Theta]]


ex[r_,\[Theta]_,\[Phi]_]:=e[r,\[Theta]][[1]] Sin[\[Theta]] Cos[\[Phi]]+e[r,\[Theta]][[2]] Cos[\[Theta]] Cos[\[Phi]]-e[r,\[Theta]][[3]] Sin[\[Phi]]
ey[r_,\[Theta]_,\[Phi]_]:=e[r,\[Theta]][[1]] Sin[\[Theta]] Sin[\[Phi]]+e[r,\[Theta]][[2]] Cos[\[Theta]] Sin[\[Phi]]+e[r,\[Theta]][[3]] Cos[\[Phi]]
ez[r_,\[Theta]_,\[Phi]_]:=e[r,\[Theta]][[1]] Cos[\[Theta]]-e[r,\[Theta]][[2]] Sin[\[Theta]]


bx[r_,\[Theta]_,\[Phi]_]:=b[r,\[Theta]][[1]] Sin[\[Theta]] Cos[\[Phi]]+b[r,\[Theta]][[2]] Cos[\[Theta]] Cos[\[Phi]]-b[r,\[Theta]][[3]] Sin[\[Phi]]
by[r_,\[Theta]_,\[Phi]_]:=b[r,\[Theta]][[1]] Sin[\[Theta]] Sin[\[Phi]]+b[r,\[Theta]][[2]] Cos[\[Theta]] Sin[\[Phi]]+b[r,\[Theta]][[3]] Cos[\[Phi]]
bz[r_,\[Theta]_,\[Phi]_]:=b[r,\[Theta]][[1]] Cos[\[Theta]]-b[r,\[Theta]][[2]] Sin[\[Theta]]


Piecewise[{
{a,r>rh&&z>0},
{b,r>rh&&z<0}},0]


PlotYZEfield[n_]:=Block[{r,\[Theta],\[Phi]=\[Pi]/2,x=0,y,z},
r[x_,y_,z_]:=Sqrt[x^2+y^2+z^2];\[Theta][x_,y_,z_]:=ArcCos[z/r[x,y,z]];
rMax=n rh;
YZey=Interpolation[Flatten[Table[{y,z,
Piecewise[{
{ey[r[x,y,z],\[Theta][x,y,z],\[Phi]],r[x,y,z]^2>(rh+0.01)^2&&z>0},
{ey[r[x,y,-z],\[Theta][x,y,-z],\[Phi]],r[x,y,-z]^2>(rh+0.01)^2&&z<0}},
0]
},{y,-rMax,rMax,rMax/100},{z,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
YZez=Interpolation[Flatten[Table[{y,z,
Piecewise[{
{ez[r[x,y,z],\[Theta][x,y,z],\[Phi]],r[x,y,z]^2>(rh+0.01)^2&&z>0},
{ez[r[x,y,-z],\[Theta][x,y,-z],\[Phi]],r[x,y,-z]^2>(rh+0.01)^2&&z<0}},
0]
},{y,-rMax,rMax,rMax/100},{z,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
YZe[y_,z_]=Sqrt[YZey[y,z]^2+YZez[y,z]^2];
DYZe[y_,z_]=D[YZe[y,z],{{y,z}}];
cP=ContourPlot[Piecewise[{{
YZe[Y,Z],Y^2+Z^2>(rh+0.01)^2}},0]
,{Y,-rMax,rMax},{Z,-rMax,rMax},Contours->8,ContourStyle->Red,ContourShading->Automatic,PlotLegends->Placed[BarLegend[Automatic,LegendMarkerSize->150,LegendMargins->{{-5,0},{-10,0}}],Right],ColorFunction->"Pastel"];
sP=StreamPlot[Piecewise[{{
DYZe[Y,Z],Y^2+Z^2>(rh+0.01)^2}},{0,0}]
,{Y,-rMax,rMax},{Z,-rMax,rMax},StreamPoints->Coarse,StreamStyle->Blue];
Show[cP,sP,FrameLabel->{y,z}]]


PlotYZBfield[n_]:=Block[{r,\[Theta],\[Phi]=\[Pi]/2,x=0,y,z},
r[x_,y_,z_]:=Sqrt[x^2+y^2+z^2];\[Theta][x_,y_,z_]:=ArcCos[z/r[x,y,z]];
rMax=n rh;
YZby=Interpolation[Flatten[Table[{y,z,
Piecewise[{
{by[r[x,y,z],\[Theta][x,y,z],\[Phi]],r[x,y,z]^2>(rh+0.01)^2&&z>0},
{by[r[x,y,-z],\[Theta][x,y,-z],\[Phi]],r[x,y,-z]^2>(rh+0.01)^2&&z<0}},
0]
},{y,-rMax,rMax,rMax/100},{z,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
YZbz=Interpolation[Flatten[Table[{y,z,
Piecewise[{
{bz[r[x,y,z],\[Theta][x,y,z],\[Phi]],r[x,y,z]^2>(rh+0.01)^2&&z>0},
{bz[r[x,y,-z],\[Theta][x,y,-z],\[Phi]],r[x,y,-z]^2>(rh+0.01)^2&&z<0}},
0]
},{y,-rMax,rMax,rMax/100},{z,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
YZb[y_,z_]=Sqrt[YZby[y,z]^2+YZbz[y,z]^2];
DYZb[y_,z_]=D[YZb[y,z],{{y,z}}];
cP=ContourPlot[Piecewise[{{
YZb[Y,Z],Y^2+Z^2>(rh+0.01)^2}},0]
,{Y,-rMax,rMax},{Z,-rMax,rMax},Contours->8,ContourStyle->Red,ContourShading->Automatic,PlotLegends->Placed[BarLegend[Automatic,LegendMarkerSize->150,LegendMargins->{{-5,0},{-10,0}}],Right],ColorFunction->"Pastel"];
sP=StreamPlot[Piecewise[{{
DYZb[Y,Z],Y^2+Z^2>(rh+0.01)^2}},{0,0}]
,{Y,-rMax,rMax},{Z,-rMax,rMax},StreamPoints->Coarse,StreamStyle->Blue];
Show[cP,sP,FrameLabel->{y,z}]]


PlotXYEfield[n_]:=Block[{r,\[Theta]=\[Pi]/2,\[Phi],x,y,z=0,rMax},
r[x_,y_,z_]:=Sqrt[x^2+y^2+z^2];\[Theta]=\[Pi]/2;\[Phi][x_,y_,z_]:=ArcTan[y/x];
rMax=n rh;
XYex=Interpolation[Flatten[Table[{x,y,
Piecewise[{{ex[r[x,y,z],\[Theta],\[Phi][x,y,z]],r[x,y,z]^2>(rh+0.01)^2&&x!=0}},0]
},{x,-rMax,rMax,rMax/100},{y,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
XYey=Interpolation[Flatten[Table[{x,y,
Piecewise[{{ey[r[x,y,z],\[Theta],\[Phi][x,y,z]],r[x,y,z]^2>(rh+0.01)^2&&x!=0}},0]
},{x,-rMax,rMax,rMax/100},{y,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
XYe[x_,y_]=Sqrt[XYex[x,y]^2+XYey[x,y]^2];
DXYe[x_,y_]=D[XYe[x,y],{{x,y}}];
cP=ContourPlot[Piecewise[{{XYe[X,Y],X^2+Y^2>(rh+0.01)^2}},0],{X,-rMax,rMax},{Y,-rMax,rMax},Contours->8,ContourStyle->Red,ContourShading->Automatic,PlotLegends->Placed[BarLegend[Automatic,LegendMarkerSize->150,LegendMargins->{{-5,0},{-10,0}}],Right],ColorFunction->"Pastel"];
sP=StreamPlot[Piecewise[{{DXYe[X,Y],X^2+Y^2>(rh+0.01)^2}},{0,0}],{X,-rMax,rMax},{Y,-rMax,rMax},StreamPoints->Coarse,StreamStyle->Blue];
Show[cP,sP,FrameLabel->{x,y}]]


PlotXYBfield[n_]:=Block[{r,\[Theta]=\[Pi]/2,\[Phi],x,y,z=0,rMax},
r[x_,y_,z_]:=Sqrt[x^2+y^2+z^2];\[Theta]=\[Pi]/2;\[Phi][x_,y_,z_]:=ArcTan[y/x];
rMax=n rh;
XYbx=Interpolation[Flatten[Table[{x,y,
Piecewise[{{bx[r[x,y,z],\[Theta],\[Phi][x,y,z]],r[x,y,z]^2>(rh+0.01)^2&&x!=0}},0]
},{x,-rMax,rMax,rMax/100},{y,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
XYby=Interpolation[Flatten[Table[{x,y,
Piecewise[{{by[r[x,y,z],\[Theta],\[Phi][x,y,z]],r[x,y,z]^2>(rh+0.01)^2&&x!=0}},0]
},{x,-rMax,rMax,rMax/100},{y,-rMax,rMax,rMax/100}],1],InterpolationOrder->1];
XYb[x_,y_]=Sqrt[XYbx[x,y]^2+XYby[x,y]^2];
DXYb[x_,y_]=D[XYb[x,y],{{x,y}}];
cP=ContourPlot[Piecewise[{{XYb[X,Y],X^2+Y^2>(rh+0.01)^2}},0],{X,-rMax,rMax},{Y,-rMax,rMax},Contours->8,ContourStyle->Red,ContourShading->Automatic,PlotLegends->Placed[BarLegend[Automatic,LegendMarkerSize->150,LegendMargins->{{-5,0},{-10,0}}],Right],ColorFunction->"Pastel"];
sP=StreamPlot[Piecewise[{{DXYb[X,Y],X^2+Y^2>(rh+0.01)^2}},{0,0}],{X,-rMax,rMax},{Y,-rMax,rMax},StreamPoints->Coarse,StreamStyle->Blue];
Show[cP,sP,FrameLabel->{x,y}]]


n=5;


plotYZEn5=PlotYZEfield[n];
plotYZBn5=PlotYZBfield[n];
plotXYEn5=PlotXYEfield[n];
plotXYBn5=PlotXYBfield[n];


plotsn5=GraphicsGrid[{{plotYZEn5,plotYZBn5},{plotXYEn5,plotXYBn5}},ImageSize->600,Spacings->{-50,-100},PlotLabel->"V="<>ToString[V]<>", w="<>ToString[w]<>", rh="<>ToString[rh]];


Export["HBH-EM-fields-n=5.pdf",plotsn5];


n=50;


plotYZEn50=PlotYZEfield[n];
plotYZBn50=PlotYZBfield[n];
plotXYEn50=PlotXYEfield[n];
plotXYBn50=PlotXYBfield[n];


plotsn50=GraphicsGrid[{{plotYZEn50,plotYZBn50},{plotXYEn50,plotXYBn50}},ImageSize->600,Spacings->{-50,-100},PlotLabel->"V="<>ToString[V]<>", w="<>ToString[w]<>", rh="<>ToString[rh]];


Export["HBH-EM-fields-n=50.pdf",plotsn50]


n=150;


plotYZEn150=PlotYZEfield[n];
plotYZBn150=PlotYZBfield[n];
plotXYEn150=PlotXYEfield[n];
plotXYBn150=PlotXYBfield[n];


plotsn150=GraphicsGrid[{{plotYZEn150,plotYZBn150},{plotXYEn150,plotXYBn150}},ImageSize->600,Spacings->{-50,-100},PlotLabel->"V="<>ToString[V]<>", w="<>ToString[w]<>", rh="<>ToString[rh]];


Export["HBH-EM-fields-n=150.pdf",plotsn150]
