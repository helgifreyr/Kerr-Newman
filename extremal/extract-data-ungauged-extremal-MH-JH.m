(* ::Package:: *)

(* 8 Martie 14 *)


 (*<<NumericalMath`ListIntegrate` *)


 
Directory[]


Off[General::spell1]
Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];

(* 1  2   3   4 *)
(* nr,V0,w,rh*)
conf=ReadList["res.txt",{Number,Number ,Number ,Number   }]

 
nr=conf[[1]][[1]];
V0= conf[[1]][[2]]; 
w= conf[[1]][[3]]
 rh= conf[[1]][[4]]; 

Print["winding number n = ",nr];
Print["rh   = ",rh];
Print["V0 = ",V0];
Print["w = ",w]; 

gr=ReadList["gridx.dat",{Number}];


lgr=Length[gr];
nx=lgr;
Print["nx = ",nx];

listar=Table[gr[[k]][[1]],{k,1,lgr}] ;
listalogr=Table[Log[10,gr[[k]][[1]]],{k,1,lgr}];

 unghi0=ReadList["gridy.dat",{Number}];
ny=Length[unghi0];
Print["ny = ",ny];
unghi=Table[unghi0[[k]][[1]],{k,1,ny}]; 

(*unghi=Table[(k-1)*Pi/2/(ny-1),{k,1,ny}];*)

ntot=nx*ny;

a=ReadList["functf.dat",{Number,Number,Number,Number ,Number,Number,Number   }];


lung1=Length[a];


(*Datele sunt salvate direct cu indexarea globala*)
F1=Table[a[[k]][[1]],{k,1,lung1}];
F2=Table[a[[k]][[2]],{k,1,lung1}];
F0=Table[a[[k]][[3]],{k,1,lung1}]; 
Z=Table[a[[k]][[4]],{k,1,lung1}]; 
W=Table[a[[k]][[5]],{k,1,lung1}]; 
Afi=Table[a[[k]][[6]],{k,1,lung1}]; 
At=Table[a[[k]][[7]],{k,1,lung1}]; 
  
(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[

F1u[k]=Table[F1[[i]],{i,(k-1)*nx+1,k*nx}];
F2u[k]=Table[F2[[i]],{i,(k-1)*nx+1,k*nx}];
F0u[k]=Table[F0[[i]],{i,(k-1)*nx+1,k*nx}]; 
Wu[k]=Table[W[[i]],{i,(k-1)*nx+1,k*nx}];
Zu[k]=Table[Z[[i]],{i,(k-1)*nx+1,k*nx}]; 
Afiu[k]=Table[Afi[[i]],{i,(k-1)*nx+1,k*nx}]; 
Atu[k]=Table[At[[i]],{i,(k-1)*nx+1,k*nx}]; 
,{k,1,ny}];
 
 
as1=2;
as2=IntegerPart[ny/2];
as3=ny-1;

sa1=3;
sa2=IntegerPart[nx/2];
sa3=nx-1;


Print["rmax = ",gr[[nx]][[1]]];


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Here I extract the relevant coefficients from the asymptotics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


ct=1;
cut=1;

ini=5; 

Do[

data=Table[{listar[[i]] ,F0u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{ 1/x ,1/x^2     } ,x];
cf01[k ]=Coefficient[u,1/x ];

data=Table[{listar[[i]] ,F1u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2       } ,x]; 
cf11[k ]=Coefficient[u,1/x ];


data=Table[{listar[[i]] ,F2u[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2     } ,x];
cf21[k ]=Coefficient[u,1/x];

 data=Table[{listar[[i]] ,Wu[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{ 1/x ,1/x^2     } ,x];
cW[k ]=Coefficient[u,1/x];


data=Table[{listar[[i]] ,Atu[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2     } ,x];
cV[k ]=Coefficient[u,1/x];

data=Table[{listar[[i]] ,Afiu[k][[i]]   },{i,nx-ini,nx-cut}];
u=Fit[data,{1/x ,1/x^2     } ,x];
cAfi[k ]=Coefficient[u,1/x];



,{k,1,ny }]

f01=Table[cf01[k],{k,1,ny }] ;
f11=Table[cf11[k],{k,1,ny }] ;
f21=Table[cf21[k],{k,1,ny }] ;
W3=Table[cW[k],{k,1,ny }] ;
cQe=Table[cV[k],{k,1,ny }] ;
mu=Table[cAfi[k]/Sin[unghi[[k]]],{k,2,ny }] ;


 constJINF=Sum[W3[[i]],{i,1,ny}]/ny


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 crucial numerical test 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

test1=1+(2f01+f11 )/(f21);

err1=Max[Abs[test1]]


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA infinity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


datx=ReadList["fx-inf.txt",Number,RecordLists->True];
lung1=Length[datx] ;


r=Table[datx[[i]][[1]],{i,1,lung1}];
infF1x=Table[datx[[i]][[2]],{i,1,lung1}];
infF2x=Table[datx[[i]][[3]],{i,1,lung1}];
infF0x=Table[datx[[i]][[4]],{i,1,lung1}];
infWx=Table[datx[[i]][[6]],{i,1,lung1}];
infAfix=Table[datx[[i]][[7]],{i,1,lung1}];
infAtx=Table[datx[[i]][[8]],{i,1,lung1}];


ctINF=Sum[infF0x[[i]],{i,1,ny}]/ny;
Print["Mass computed from the data at infinity = ",ctINF]
 

JINF=Sum[infWx[[i]],{i,1,ny}]/ny;
Print["J computed from the data at infinity = ",JINF];


QeINF=Sum[infAtx[[i]],{i,1,ny}]/ny;
Print["Qe computed from the data at infinity   = ",QeINF];

muINF=Sum[infAfix[[i]]/Sin[unghi[[i]]],{i,2,ny}]/(ny-1);
Print["Qm computed from the data at infinity   = ",muINF];



constINF=ctINF


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
computation Mass from asymptotics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(* th result:*) 
(*
Series[-g[4,4],{r,Infinity,1}]
1+(2 const-rh)/r+O[1/r]^2: non-extremal

Series[gtt,{r,Infinity,1}]
-1+(-2 ct+2 rh)/r+O[1/r]^2: extremal
*)
(* atentie!! la extremal este diferit! *)


const=Sum[f01[[i]],{i,1,ny}]/ny;

Mc=  constINF;

(* non-extremal:
Mc=  constINF;
MSch=rh/2 ;
Mass=MSch+Mc;
*)


MSch=rh ;
Mass=MSch+Mc;

Print["Mass Schw     = ",MSch];
Print["Mass correction= ",Mc];
Print["total Mass     = ",Mass];

(*Print[ MSAdS/Mc ]*)


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE  DATA t-0 -- no conical singularities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


dat=ReadList["f-t0.txt",Number,RecordLists->True];
lung1=Length[dat] ;

  
r=Table[dat[[i]][[1]],{i,1,lung1}];
t0F1=Table[dat[[i]][[2]],{i,1,lung1}];
t0F2=Table[dat[[i]][[3]],{i,1,lung1}];
t0F0=Table[dat[[i]][[4]],{i,1,lung1}];
t0W=Table[dat[[i]][[6]],{i,1,lung1}];
t0Z=Table[dat[[i]][[5]],{i,1,lung1}];
t0Afi=Table[dat[[i]][[7]],{i,1,lung1}];
t0At=Table[dat[[i]][[8]],{i,1,lung1}];

ratio= t0F1 -t0F2;
Max[Abs[ratio]]





(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PRELUCRARE HORIZON DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


dat=ReadList["f-0.txt",Number,RecordLists->True];
lung1=Length[dat] ;

unghi=Table[dat[[i]][[1]],{i,1,lung1}];
hF1=Table[dat[[i]][[2]],{i,1,lung1}]; 
hF2=Table[dat[[i]][[3]],{i,1,lung1}]; 
hF0=Table[dat[[i]][[4]],{i,1,lung1}]; 
hZ=Table[dat[[i]][[5]],{i,1,lung1}]; 
hW=Table[dat[[i]][[6]],{i,1,lung1}]; 
hAfi=Table[dat[[i]][[7]],{i,1,lung1}];
hAt=Table[dat[[i]][[8]],{i,1,lung1}];



(*%%%%%%% event horizon area %%%%%%%%*)
AH0=4 Pi rh^2;

iAHc=Table[{unghi[[k]],1/2  Sin[unghi[[k]]] E^((hF1[[k]] +hF2[[k]] )) },{k,1,lung1 }];

(* 2 because I integrate between 0, Pi/2 *)
 AHc= 2Integrate[Interpolation[iAHc,InterpolationOrder->1][x],{x,0,Pi/2}];

AH=AH0 AHc;

Print["Schw area  AH0= ",AH0];
Print["correction AHc= ",AHc];
Print["Event horizon area =",AH];







(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
computation Le, Lp -- see MATH code for derivation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(*Le=2 \[ExponentialE]^F2[0,\[Pi]/2] \[Pi] rH*)
Le=2 E^hF2[[ny]] \[Pi] rh;
Print["Le = ",Le ];

(*Lp=2 \!\(
\*SubsuperscriptBox[\(\[Integral]\), \(0\), \(\[Pi]\)]\(\(
\*SuperscriptBox[\(\[ExponentialE]\), \(F1[0, t]\)]\ rH\)\[DifferentialD]t\)\)*)
Lp1= Table[{unghi[[k]],2 rh    E^hF1[[k]]   },{k,1,lung1 }];
(*Lp=  2ListIntegrate[Lp1,2]//N;(*factor 2: because I integrate [0,pi/2]*) *)
Lp=2 Integrate[Interpolation[Lp1,InterpolationOrder->3][x],{x,Min[Lp1[[All,1]]],Max[Lp1[[All,1]]]}];
Print["Lp = ",Lp ];
Print["  " ];


(* computation MASS from the energy momentum tensor - Smarr relation *)

asa1=2;
asa2=IntegerPart[ny/2];
asa3=ny-1;


 (* ordinea este: { T34,T44,Ttot }  *)
q=ReadList["T44.dat",{Number,Number,Number }];
lungq=Length[q];

diference=lungq-ntot;
Print["It must be zero! ",diference];
 

T34=Table[q[[k]][[1]],{k,1,lungq}]; 
ro=Table[q[[k]][[2]],{k,1,lungq}]; 
 Ttot=Table[q[[k]][[3]],{k,1,lungq}]; 

(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[ 
 
 T34u[k]=Table[T34[[i]],{i,(k-1)*nx+1,k*nx}];  
rou[k]=Table[ro[[i]],{i,(k-1)*nx+1,k*nx}]; 
Ttotu[k]=Table[Ttot[[i]],{i,(k-1)*nx+1,k*nx}]; 
 
,{k,1,ny}] 



ni=2;


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I compute the scalar field contribution to the total mass
& Smarr law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(* sqrt = \[ExponentialE]^(F0[r,t]+2 F1[r,t]+F2[r,t]) r Sqrt[g[r]] Sin[t]*)

(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)

Do[

 Mio2[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2]  Ttotu[k][[i]]},{i,ni,nx-1}];

 Mio3[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2]  T34u[k][[i]]},{i,ni,nx-1}];

,{k,2,ny-1}];




(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)
Do[
(*Ma2[k]=ListIntegrate[Mio2[k],2]//N;
Ma3[k]=ListIntegrate[Mio3[k],2]//N;*)
Ma2[k]=Integrate[Interpolation[Mio2[k] 1,InterpolationOrder->3][x],{x,Min[Mio2[k][[All,1]]],Max[Mio2[k][[All,1]]]}]//N;
Ma3[k]=Integrate[Interpolation[Mio3[k] 1,InterpolationOrder->3][x],{x,Min[Mio3[k][[All,1]]],Max[Mio3[k][[All,1]]]}]//N;
 ,{k,2,ny-1}];

 
Ma2[1]=Ma2[2];
Ma2[ny]=Ma2[ny-1];
 
Ma3[1]=Ma3[2];
Ma3[ny]=Ma3[ny-1];

 Minn2=Table[{unghi[[k]], Sin[unghi[[k]]] Ma2[k]},{k,1,ny}];
  Minn3=Table[{unghi[[k]], Sin[unghi[[k]]] Ma3[k]},{k,1,ny}];

Mint=Integrate[Interpolation[Minn2,InterpolationOrder->3][x],{x,Min[Minn2[[All,1]]],Max[Minn2[[All,1]]]}]//N;
Jint=Integrate[Interpolation[Minn3,InterpolationOrder->3][x],{x,Min[Minn3[[All,1]]],Max[Minn3[[All,1]]]}]//N;
 (*Mint=ListIntegrate[Minn2,2]//N; 
Jint=ListIntegrate[Minn3,2]//N; *)
 Print["Mintegral = ",Mint];  
 Print["Jintegral = ",Jint];  



(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
computation MH, JH -- see MATH code for derivation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
dat=ReadList["fxx-0.txt",Number,RecordLists->True];
lung1=Length[dat] ;

unghi=Table[dat[[i]][[1]],{i,1,lung1}];
hF1xx=Table[dat[[i]][[2]],{i,1,lung1}]; 
hF2xx=Table[dat[[i]][[3]],{i,1,lung1}]; 
hF0xx=Table[dat[[i]][[4]],{i,1,lung1}]; 
hZxx=Table[dat[[i]][[5]],{i,1,lung1}]; 
hWxx=Table[dat[[i]][[6]],{i,1,lung1}]; 

hW2= hWxx;

(* tr44EH=tr44EH=(\[ExponentialE]^(-f00[t]+3 f20[t]) w0 Sin[t]^3 (-w0+rh^2 w2[t]))/rh;*)
 itr44EH=Table[{unghi[[k]], 1/rh  E^(-hF0[[k]] +3hF2[[k]] )  hW[[k]] Sin[unghi[[k]]]^3(-hW[[k]] +(1/2)rh^2 hW2[[k]]) },{k,1,lung1 }];
 

(*tr34EH=-\[ExponentialE]^(-f00[t]+3 f20[t]) rh Sin[t]^3 (-w0+rh^2 w2[t]);*)
 itr34EH=Table[{unghi[[k]], -E^(-hF0[[k]] +3hF2[[k]] )  rh Sin[unghi[[k]]]^3(-hW[[k]] +(1/2)rh^2 hW2[[k]])  },{k,1,lung1 }];






(* 2 because I integrate between 0, Pi/2 *)
(*
MH= 2 2 Pi Integrate[Interpolation[itr44EH,InterpolationOrder\[Rule]1][x],{x,0,Pi/2}];
 JH= 2 2 Pi   Integrate[Interpolation[itr34EH,InterpolationOrder\[Rule]1][x],{x,0,Pi/2}];

(* I consider without factor of 4 Pi because I take 4 Pi G=2 2 Pi 
*)
*)

MH=- Integrate[Interpolation[itr44EH,InterpolationOrder->1][x],{x,0,Pi/2}];JH=  1/2 Integrate[Interpolation[itr34EH,InterpolationOrder->1][x],{x,0,Pi/2}];
 


(* MH=2 TH1/4 AH-2  w   (1/2constJINF -Jint) *)
alfa=1;

J=1/2 constJINF;
MInt=-Mint 1/2 2alfa^2;

errM=1-(MH+MInt)/Mass;
errJ=1-(JH+Jint)/J;
Print["M   = ",Mass];
Print["MInt= ",MInt];
Print["MH  = ",MH];
Print["Merror = ",errM];

Print[""];

Print["J   = ",J];
Print["Jint= ",Jint];
Print["JH  = ",JH];
Print["Jerror= ",errJ];










alfa=1;

(* Smarr law Schw BH;
MSch-2 TH01/4 AH0 
 *)

Print["total mass     = ",Mass];
Print["total J     = ",1/2 constJINF ];
Print["mass integral  = ",-Mint 1/2 2alfa^2];
Print["J integral  = ",Jint];
Print["AH    = ",  AH];
Print["factor TH S    = ",2 TH 1/4 AH];
Print["factor OmegaHJ = ",2  w   (1/2 constJINF -Jint)];

Print[" "];

(* smarr relation *)
errSmarr=1-(Mass -2  w   (1/2 constJINF -Jint)  )/(-Mint 1/2 2alfa^2);
Print[" errSmarr= ",errSmarr]; 


TH=0;
g=-2 Mass (muINF /(QeINF J))  


asa= Table[{ w ,V0,Mass ,MH,J ,JH,Mint, Jint, TH,AH,QeINF,muINF,rh,g}] 






stmp=OpenAppend["tmp.txt"];
Write[stmp,asa];
Close[stmp] ;
 


(* gyromagnetic ratio*)


