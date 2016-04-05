(* ::Package:: *)

(* 17 Jan 16 *)



Off[SetDelayed::"write"]


 (*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	          extract global data gauged BSs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)



 (*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	       Here I read the parameters + the functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


Directory[]


Off[General::spell1]
Remove["Global`*"];
Unprotect[In,Out];
Clear[In,Out];

(* 1  2   3 4*)
(* nr,V0,w,rh*)
conf=ReadList["res.txt",{Number,Number ,Number,Number   }];

 
nr=conf[[1]][[1]];
V0= conf[[1]][[2]]; 
wf= conf[[1]][[3]];
rh= conf[[1]][[4]]; 

Print["winding number n = ",nr];
Print["frequency w = ",wf];
Print["V0 = ",V0]; 
Print["rh = ",rh]; 


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

a=ReadList["functf.dat",{Number,Number,Number,Number,Number ,Number,Number }];


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


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
plot quantities which show the quality of numerics
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

 
 (* ordinea este: {Eq11,Eq12,Ricci,Kr}  *)
q=ReadList["eq1.dat",{Number,Number,Number,Number  } ];
lungq=Length[q];

diference=lungq-ntot;
Print["It must be zero! ",diference];
 

Eq11=Table[q[[k]][[1]],{k,1,lungq}]; 
Eq12=Table[q[[k]][[2]],{k,1,lungq}]; 
 Ricci=Table[q[[k]][[3]],{k,1,lungq}]; 
Kr=Table[q[[k]][[4]],{k,1,lungq}];   

(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[ 
Eq11u[k]=Table[Eq11[[i]],{i,(k-1)*nx+1,k*nx}]; 
Eq12u[k]=Table[Eq12[[i]],{i,(k-1)*nx+1,k*nx}]; 
 Ricciu[k]=Table[Ricci[[i]],{i,(k-1)*nx+1,k*nx}]; 
 Kru[k]=Table[Kr[[i]],{i,(k-1)*nx+1,k*nx}];  
 
,{k,1,ny}]


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 Compute mass&co from asymptotics 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


ct=1;
cut=1;

ini= 5; 

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
u=Fit[data,{1/x,1/x^2       } ,x];
cw2[k ]=Coefficient[u,1/x];

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
w2=Table[cw2[k],{k,1,ny }] ;
cQe=Table[cV[k],{k,1,ny }] ;
mu=Table[cAfi[k]/Sin[unghi[[k]]],{k,2,ny }] ;
 


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 crucial numerical test 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

test1=1+(2f01+f11 )/(f21);

err1=Max[Abs[test1]]


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   DATA infinity -first derivatives
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
Print["Mass computed from the data at infinity = ",ctINF];

JINF=Sum[infWx[[i]],{i,1,ny}]/ny;
Print["J computed from the data at infinity = ",JINF];


QeINF=Sum[infAtx[[i]],{i,1,ny}]/ny;
Print["Qe computed from the data at infinity   = ",QeINF];

muINF=Sum[infAfix[[i]]/Sin[unghi[[i]]],{i,2,ny}]/(ny-1);
Print["Qm computed from the data at infinity   = ",muINF];


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
computation Mass from asymptotics
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)

(* th result:*) 
(*
Series[-g[4,4],{r,Infinity,1}]
1+(2 const-rh)/r+O[1/r]^2
*)


Mc=  ctINF;

MSch=rh/2 ;
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
 HORIZON DATA
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


(*%%%%%%% Hawking temperature %%%%%%%%*)
TH0=1/(4 Pi rh) ;

THc=  1/lung1 Sum[ E^( (hF0[[i]]-hF1[[i]])),{i,1,lung1}];(* correction due to the scalar field*)
TH=TH0 THc;
Print["Schw temp. TH0= ",TH0]
Print["correction THc= ",THc]
Print["TH= ",TH]


errTH=1-Abs[Min[hF0-hF1] /Max[hF0-hF1]];
Print["error TH= ",errTH//N]




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







(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 computation MASS from the energy momentum tensor - Smarr relation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)


asa1=2;
asa2=IntegerPart[ny/2];
asa3=ny-1;


 (* ordinea este: {T34,T34M,T34s,T44,Ttot,TtotM,Ttots}  *)
q=ReadList["T44.dat",{Number,Number,Number,Number,Number,Number,Number }];
lungq=Length[q];

diference=lungq-ntot;
Print["It must be zero! ",diference];
 

T34=Table[q[[k]][[1]],{k,1,lungq}]; 
T34M=Table[q[[k]][[2]],{k,1,lungq}]; 
T34s=Table[q[[k]][[3]],{k,1,lungq}]; 
ro=Table[q[[k]][[4]],{k,1,lungq}]; 
 Ttot=Table[q[[k]][[5]],{k,1,lungq}]; 
 TtotM=Table[q[[k]][[6]],{k,1,lungq}]; 
 Ttots=Table[q[[k]][[7]],{k,1,lungq}]; 

(*Se construiesc ny liste pt. marimi de interes la unghiuri fixate *)
(*foarte util in reprezentari grafice *)

Do[ 
T34u[k]=Table[T34[[i]],{i,(k-1)*nx+1,k*nx}]; 
T34Mu[k]=Table[T34M[[i]],{i,(k-1)*nx+1,k*nx}]; 
T34su[k]=Table[T34s[[i]],{i,(k-1)*nx+1,k*nx}]; 
rou[k]=Table[ro[[i]],{i,(k-1)*nx+1,k*nx}]; 
Ttotu[k]=Table[Ttot[[i]],{i,(k-1)*nx+1,k*nx}]; 
TtotMu[k]=Table[TtotM[[i]],{i,(k-1)*nx+1,k*nx}]; 
Ttotsu[k]=Table[Ttots[[i]],{i,(k-1)*nx+1,k*nx}]; 
(*	Print[T44u[k]];*)
,{k,1,ny}]

ni=2;


(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 now I compute the scalar field contribution to the total mass
& Smarr law
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
(* sqrt = \[ExponentialE]^(F0[r,t]+2 F1[r,t]+F2[r,t]) r Sqrt[g[r]] Sin[t]*)

(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)

Do[

 Mio2[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2 ]  Ttotu[k][[i]]},{i,ni,nx-1}];

 Mio2M[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2 ]  TtotMu[k][[i]]},{i,ni,nx-1}];

 Mio2s[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2 ]  Ttotsu[k][[i]]},{i,ni,nx-1}];


Jio2[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2 ]  T34u[k][[i]]},{i,ni,nx-1}];

Jio2M[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2 ]  T34Mu[k][[i]]},{i,ni,nx-1}];

Jio2s[k]=Table[{listar[[i]], E^(F0u[k][[i]]+2 F1u[k][[i]]+F2u[k][[i]])  listar[[i]]  Sqrt[listar[[i]]^2+rh^2 ]  T34su[k][[i]]},{i,ni,nx-1}];




,{k,2,ny-1}];



(*Se construiesc ny liste pt. integralele marimilor de interes la unghiuri fixate *)
Do[
(*Ma2[k]=ListIntegrate[Mio2[k],2]//N;
Ma2M[k]=ListIntegrate[Mio2M[k],2]//N;
Ma2s[k]=ListIntegrate[Mio2s[k],2]//N;*)
Ma2[k]=Integrate[Interpolation[Mio2[k]1,InterpolationOrder->3][x],{x,Min[Mio2[k][[All,1]]],Max[Mio2[k][[All,1]]]}]//N;
Ma2M[k]=Integrate[Interpolation[Mio2M[k]1,InterpolationOrder->3][x],{x,Min[Mio2M[k][[All,1]]],Max[Mio2M[k][[All,1]]]}]//N;
Ma2s[k]=Integrate[Interpolation[Mio2s[k]1,InterpolationOrder->3][x],{x,Min[Mio2s[k][[All,1]]],Max[Mio2s[k][[All,1]]]}]//N;
(*Ja2[k]=ListIntegrate[Jio2[k],2]//N;
Ja2M[k]=ListIntegrate[Jio2M[k],2]//N;
Ja2s[k]=ListIntegrate[Jio2s[k],2]//N;*)
Ja2[k]=Integrate[Interpolation[Jio2[k]1,InterpolationOrder->3][x],{x,Min[Jio2[k][[All,1]]],Max[Jio2[k][[All,1]]]}]//N;
Ja2M[k]=Integrate[Interpolation[Jio2M[k]1,InterpolationOrder->3][x],{x,Min[Jio2M[k][[All,1]]],Max[Jio2M[k][[All,1]]]}]//N;
Ja2s[k]=Integrate[Interpolation[Jio2s[k]1,InterpolationOrder->3][x],{x,Min[Jio2s[k][[All,1]]],Max[Jio2s[k][[All,1]]]}]//N;
 ,{k,2,ny-1}];

 
Ma2[1]=Ma2[2];
Ma2[ny]=Ma2[ny-1];
Ma2M[1]=Ma2M[2];
Ma2M[ny]=Ma2M[ny-1];
Ma2s[1]=Ma2s[2];
Ma2s[ny]=Ma2s[ny-1];

 
Ja2[1]=Ja2[2];
Ja2[ny]=Ja2[ny-1];
Ja2M[1]=Ja2M[2];
Ja2M[ny]=Ja2M[ny-1];
Ja2s[1]=Ja2s[2];
Ja2s[ny]=Ja2s[ny-1];



 Minn2=Table[{unghi[[k]], Sin[unghi[[k]]] Ma2[k]},{k,1,ny}];
 Minn2M=Table[{unghi[[k]], Sin[unghi[[k]]] Ma2M[k]},{k,1,ny}];
 Minn2s=Table[{unghi[[k]], Sin[unghi[[k]]] Ma2s[k]},{k,1,ny}];
Jinn2=Table[{unghi[[k]], Sin[unghi[[k]]] Ja2[k]},{k,1,ny}];
Jinn2M=Table[{unghi[[k]], Sin[unghi[[k]]] Ja2M[k]},{k,1,ny}];
Jinn2s=Table[{unghi[[k]], Sin[unghi[[k]]] Ja2s[k]},{k,1,ny}];

 
 (*Mint=ListIntegrate[Minn2,2]//N; 
 MintM=ListIntegrate[Minn2M,2]//N; 
 Mints=ListIntegrate[Minn2s,2]//N; *)
Mint=Integrate[Interpolation[Minn2,InterpolationOrder->3][x],{x,Min[Minn2[[All,1]]],Max[Minn2[[All,1]]]}]//N;
MintM=Integrate[Interpolation[Minn2M,InterpolationOrder->3][x],{x,Min[Minn2M[[All,1]]],Max[Minn2M[[All,1]]]}]//N;
Mints=Integrate[Interpolation[Minn2s,InterpolationOrder->3][x],{x,Min[Minn2s[[All,1]]],Max[Minn2s[[All,1]]]}]//N;
(*Jint=ListIntegrate[Jinn2,2]//N; 
JintM=ListIntegrate[Jinn2M,2]//N; 
Jints=ListIntegrate[Jinn2s,2]//N; *)
Jint=Integrate[Interpolation[Jinn2,InterpolationOrder->3][x],{x,Min[Jinn2[[All,1]]],Max[Jinn2[[All,1]]]}]//N;
JintM=Integrate[Interpolation[Jinn2M,InterpolationOrder->3][x],{x,Min[Jinn2M[[All,1]]],Max[Jinn2M[[All,1]]]}]//N;
Jints=Integrate[Interpolation[Jinn2s,InterpolationOrder->3][x],{x,Min[Jinn2s[[All,1]]],Max[Jinn2s[[All,1]]]}]//N;

 Print["Mintegral = ",Mint];  
 Print["Mintegral (Maxwell) = ",MintM];  
 Print["Mintegral (scalar) = ",Mints];  
 Print["Jintegral = ",Jint]; 
 Print["Jintegral (Maxwell) = ",JintM]; 
 Print["Jintegral (scalar) = ",Jints]; 



(*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
computation MH, JH -- see MATH code for derivation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*)
dat=ReadList["fxx-0.txt",{Number,Number ,Number ,Number,Number ,Number,Number,Number  }];
lung1=Length[dat] ;

unghi=Table[dat[[i]][[1]],{i,1,lung1}];
hF1xx=Table[dat[[i]][[2]],{i,1,lung1}]; 
hF2xx=Table[dat[[i]][[3]],{i,1,lung1}]; 
hF0xx=Table[dat[[i]][[4]],{i,1,lung1}]; 
hZxx=Table[dat[[i]][[5]],{i,1,lung1}]; 
hWxx=Table[dat[[i]][[6]],{i,1,lung1}]; 

hW2= hWxx;

(* tr44EH=-(1/2) \[ExponentialE]^(f00[t]+f20[t]) rh Sin[t]+(\[ExponentialE]^(-f00[t]+3 f20[t]) w0 Sin[t]^3 (-w0+rh^2 w2[t]))/rh;;*)
 itr44EH=Table[{unghi[[k]],-(1/2) E^(hF0[[k]] +hF2[[k]] ) rh Sin[unghi[[k]]]+  1/rh  E^(-hF0[[k]] +3hF2[[k]] )  hW[[k]] Sin[unghi[[k]]]^3(-hW[[k]] +(1/2)rh^2 hW2[[k]]) },{k,1,lung1 }];
 it=Table[{unghi[[k]],-(1/2) E^(hF0[[k]] +hF2[[k]] ) rh Sin[unghi[[k]]]},{k,1,lung1 }];



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
qt=   -Integrate[Interpolation[it ,InterpolationOrder->1][x],{x,0,Pi/2}];






(* MH=2 TH1/4 AH-2  w   (1/2constJINF -Jint) *)


J= -(1/2)JINF;
MInt=-Mint 1/2 2 ;

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






(* Smarr law Schw BH;
MSch-2 TH01/4 AH0 
 *)

Print["total mass     = ",Mass];
Print["mass integral  = ",-Mint 1/2 2 ];
Print["factor TH S    = ",2 TH 1/4 AH];
Print["factor OmegaHJ = ",2  wf   (J -Jint)];

Print[" "];

(* smarr relation *)
errSmarr=1-(Mass-2 TH 1/4 AH-2  wf   (J -Jint) )/(-Mint 1/2 2 );
Print[" errSmarr= ",errSmarr]; 


(* gyromagnetic ratio *)
g=-2 Mass (muINF /(QeINF J))  


 asa=Table[{ wf,V0,Mass ,MH,J ,JH,Mint,MintM,Mints,Jint,JintM,Jints,TH,AH,QeINF,muINF,rh,g}] 


stmp=OpenAppend["tmp.txt"];
Write[stmp,asa];
Close[stmp];

