$ontext
=================================================================================================
Modelo de destilaci�n extractiva para la obtenci�n de etanol anhidro con glicerol como entrainer
Solucionado mediante variables externas y el algoritmo D-SDA
=================================================================================================

$offtext
*-------------------------------------------------------------------------------
*                                Seccion 1
*           Conjuntos para definir operacion en estado estable
*-------------------------------------------------------------------------------
*Se usa solo el elemento inicial con el punto de colocacion inicial: estado estable
sets j "1 punto de colocacion (0)" /1/
     N "1 elemento finito" /1/;
*-------------------------------------------------------------------------------
*                                Seccion 2
*       Conjuntos, variables, par�metros y ecuaciones principales del sistema
*-------------------------------------------------------------------------------
*N�mero m�ximo de etapas
*Conjuntos
set comp "lista de componentes que intervienen en el sistema" /wat, eth, gly/;
set Net "Todas las etapas de la columna reales y sobrantes incluyendo con y rev" /1*30/;
alias(Net,Net1);

*Variables principales
positive variables
L(N,j,Net)  "Flujo de l�quido [mol/hr]"
V(N,j,Net)  "Flujo de vapor [mol/hr]"
x(N,j,comp,Net)     "Porcentaje molar en el l�quido [%]"
y(N,j,comp,Net)     "Porcentaje molar en el vapor [%]"
Temp(N,j,Net)       "Temperatura de operaci�n [K]"
P(N,j,Net)  "Presi�n por etapa [atm]"
RR(N,j)      "Relaci�n molar de reflujo [-]"
Qc(N,j)      "Carga t�rmica del condensador [kJ/hr]"
Qr(N,j)      "Carga t�rmica del rehervidor [kJ/hr]"
BR(N,j)       "Boil up [-]"
;


*Par?metros hidr?ulicos
parameter
d_hole      "Di�metro de los agujeros [m]"  /0.0127/
tray_t      "Espesor del plato [m]" /0.002/
hw      "Weir height [m]" /0.0254/
Lw      "Weir length [m]" /0.578/
HS      "Altura de cada plato [m]" /0.61/
Sfactor "Factor de seguridad altura de la columna [-]" /0.15/
K0    "Coeficiente de orificio [-]"
;
K0=(880.6-(67.7*d_hole/tray_t)+(7.32*((d_hole/tray_t)**2))-(0.338*((d_hole/tray_t)**3)))*1e-3;

*Variables hidraulicas
positive variable D "Di�metro de la columna [m]";
positive variable Htotal "Altura total de la columna [m]";
positive variable At "Area activa [m2]";
positive variable Ad  "Area de derramadero [m2]";
positive variable A0 "Area agujerada [m2]";
positive variable poro  "Porosidad del plato [-]";
Positive variable pitch  "Distancia entre agujeros [m]";
Positive variable A_col "(m2)";

*Ecuaciones hidraulicas
equation EqAt;
EqAt.. At=e=sqr(D/2)*(pi-(1.854590-0.96));
At.l= 0.299;
equation EqAd;
EqAd.. Ad=e=sqr(D/2)*(0.5*(1.854590-0.96));
Ad.l=0.1*0.37;
Equation EqPitch;
EqPitch.. pitch=e=sqr(D/2)*pi*0.12;
pitch.l=0.12*0.6;
equation EqA0;
EqA0.. A0=e=At*poro;
A0.l=0.3*0.04;
Equation EqPoro;
EqPoro.. poro=e=0.907*sqr(d_hole/pitch);
poro.l=0.907*sqr(d_hole/(0.12*0.5));
Equation EqDiam;
EqDiam..A_col=e=sqr(D/2)*pi;
D.l = 1.23;

*Alimentacion  azeotropica
parameter FAz "Flujo de alimentaci�n de aze�tropo [kmol/h]" /100/;
parameter zAz(N,j,comp) "Porcentaje molar en la alimentaci�n de butenos";
zAz(N,j,'eth')=85;
zAz(N,j,'wat')=100-zAz(N,j,'eth');
zAz(N,j,'gly')=0;

*Alimentacion glicerol
*variable FG  "Flujo de alimentaci�n de glicerol [kmol/h]";
parameters
FG  "Flujo de alimentaci�n de glicerol [kmol/h]" /52/
zg(comp) "Porcentaje molar en la alimentaci?n de glicerol"
/
wat 0
eth 0
gly 100
/
;
*Parametros de operaci�n
parameter
Pop     "Presi�n de operaci�n condensador [atm]"    /1/
Tali1   "Temperatura de alimentaci�n Feed 1 (Glicerol)[K]"  /305/
Tali2   "Temperatura de alimentaci�n Feed 2 (Azeotr�pica) [K]"  /293.15/
xReth  "Composici�n molar de etanol en fondos deseada"   /99.5/
R "Constante de los gases [kJ/K*kmol]" /8.31447215/
;

*-------------------------------------------------------------------------------
*                                Secci?n 4
*                          Restricciones de pureza
*-------------------------------------------------------------------------------

equations pureza0(Net);
pureza0(Net)$(ord(Net) eq 1).. x('1','1','eth',Net)=g=xReth;

*-------------------------------------------------------------------------------
*                                Secci?n 5
*    C?lculo de presiones de saturaci?n por medio de la ecuaci?n de Antoine
*-------------------------------------------------------------------------------
*Constantes de la ecuaci?n de Antoine expandida
parameters
C1a(comp)
/
wat      62.12291155
eth      61.77791155
gly      88.45991155
/
C2a(comp)
/
wat      -7258.200000
eth      -7122.300000
gly      -13808.00000
/
C3a(comp)
/
wat      0
eth      0
gly      0
/
C4a(comp)
/
wat      0
eth      0
gly      0
/
C5a(comp)
/
wat      -7.303700000
eth      -7.142400000
gly      -10.08800000
/
C6a(comp)
/
wat      4.16530000e-6
eth      2.88530000e-6
gly      3.5712000e-19
/
C7a(comp)
/
wat      2.0
eth      2.0
gly      6.0
/
;
positive variables Psat(N,j,comp,Net) presion de saturacion (atm);
equations EqPsat(N,j,comp,Net);
EqPsat(N,j,comp,Net).. Psat(N,j,comp,Net)=e=1/1.01325*exp( C1a(comp) + (C2a(comp)/(Temp(N,j,Net)+C3a(comp))) + (C4a(comp)*Temp(N,j,Net)) + (C5a(comp)*log(Temp(N,j,Net)) + (C6a(comp)*(Temp(N,j,Net)**C7a(comp)))) );

*-------------------------------------------------------------------------------
*                                Secci?n 6 *Buscar en Aspen*
*     C?lculo de densidades de l?quido por medio de la ecuaci?n IK-CAPI
*     C?lculo de densidades de l?quido por medio de la ecuaci?n DIPPR cr?tica
*     C?lculo de densidades de gas por medio de ecuaci?n de gas ideal corregida
*-------------------------------------------------------------------------------
*Constantes de la ecuaci?n DIPPR
parameters
MW(comp) "Peso molecular [kg/kmol]"
/
wat      18.01528
eth      46.07
gly      92.09382
/
Tcrit(comp) "Temperatura cr?tica [K]"
/
wat      647.0960000
eth      514.0000000
gly      850.0000000
/
Pcrit(comp) "Presi?n cr?tica [bar]"
/
wat     220.6351
eth     63.00
gly     75.00
/

C1r(comp) "kmol/cum"
/
wat     17.863
eth     1.6288
gly     0.92382
/
C2r(comp) "kmol/cum"
/
wat     58.606
eth     0.27469
gly     0.24386
/
C3r(comp) "kmol/cum"
/
wat     -95.396
eth     515
gly     850
/
C4r(comp) "kmol/cum"
/
wat     213.89
eth     0.23178
gly     0.22114
/
c5r "kmol/cum" /-141.26/
;
positive variable Tcritm(N,j,Net);

equation EqTcritm(N,j,Net);
EqTcritm(N,j,Net).. Tcritm(N,j,Net) =e= (sqr(sum(comp,(x(N,j,comp,Net)/100)*Tcrit(comp)/(Pcrit(comp)**0.5))))/(sum(comp,(x(N,j,comp,Net)/100)*Tcrit(comp)/Pcrit(comp)));

positive variables rho(N,j,comp,Net) "Densidad molar de l?quido [mol/m^3]";
equation Eqrho1(N,j,comp,Net), Eqrho2(N,j,Net);
Eqrho1(N,j,comp,Net)$(ord(comp) ge 2).. rho(N,j,comp,Net)=e=( C1r(comp)/(C2r(comp)**(1+((1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**C4r(comp)))) )*1000;
Eqrho2(N,j,Net).. rho(N,j,'wat',Net)=e=( C1r('wat')+(C2r('wat')*(1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**(0.35))+(C3r('wat')*(1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**(2/3))
                                            +(C4r('wat')*(1-(Temp(N,j,Net)/Tcritm(N,j,Net)))) + (C5r*(1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**(4/3)) )*1000;

rho.l(N,j,comp,Net)$(ord(comp) ge 2)=15000;
rho.l(N,j,comp,Net)$(ord(comp) eq 1)=50000;

positive variable rhoV(N,j,Net) "Densidad molar de vapor [mol/m^3]";
equation EqurhoV(N,j,Net);
EqurhoV(N,j,Net).. rhoV(N,j,Net)=e=P(N,j,Net)/(R/101325*Temp(N,j,Net));
rhoV.l(N,j,Net)=60;

positive variables
Qliq(N,j,Net) "Flujo volum�trico de l�quido (m3/hr)"
Qvap(N,j,Net) "Flujo volum�trico de vapor (m3/hr)"
;

equations EqQliq(N,j,Net), EqQvap(N,j,Net);
EqQliq(N,j,Net)..Qliq(N,j,Net)=e=L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100);
Qliq.l(N,j,Net)=6;

EqQvap(N,j,Net)..Qvap(N,j,Net)=e=V(N,j,Net)/rhoV(N,j,Net);
Qvap.l(N,j,Net)=4100;

*-------------------------------------------------------------------------------
*                                Secci?n 7
*     C?lculo de tensi?n superficial por medio de la ecuaci?n DIPPR cr?tica
*-------------------------------------------------------------------------------
*Constantes de la ecuaci?n DIPPR
parameters
C1sig(comp)
/
wat     0.17766
eth     0.0241005
gly     0.0645335
/
C2sig(comp)
/
wat     2.567
eth     -7.75658e-05
gly     -5.38024e-5
/
C3sig(comp)
/
wat     -3.3377
eth     -1.025e-07
gly     -2.1558e-7
/
C4sig(comp)
/
wat     1.9699
eth     0
gly     0
/
;
positive variables sigma(N,j,Net) "Tensi?n superficial l?quido vapor [N/m]";
equation Eqsigma(N,j,Net);
Eqsigma(N,j,Net).. sigma(N,j,Net)=e=sum(comp,(x(N,j,comp,Net)/100)*C1sig(comp)*(1-(Temp(N,j,Net)/Tcritm(N,j,Net)))**(C2sig(comp)+C3sig(comp)*(Temp(N,j,Net)/Tcritm(N,j,Net))+C4sig(comp)*((Temp(N,j,Net)/Tcritm(N,j,Net)))**2));
*------------------------------------------------------------------------------
*                                Secci?n 8
*          C?lculo de coeficientes de actividad por medio del modelo NRTL
*-------------------------------------------------------------------------------
table a_nrtl(comp,comp) Par?metro a de NRTL
        wat             eth             gly
wat                     3.4578          -1.2515
eth     -0.8009                         0.0
gly     -0.7318         0.0
;

table b_nrtl(comp,comp) Par?metro b de NRTL
        wat             eth             gly
wat                     -586.0809       272.6075
eth     246.18                          442.713
gly     170.9167        36.139
;

table c_nrtl(comp,comp) Par?metro c de NRTL
        wat             eth             gly
wat                     0.3             0.3
eth     0.3             0               0.3
gly     0.3             0.3             0
;
alias (comp,comp1);
parameter alfa_nrtl(comp,comp);
alfa_nrtl(comp,comp1)$(ord(comp) ne ord(comp1))=c_nrtl(comp,comp1);

*Par?metros G y Tao
variables tao_nrtl(N,j,comp,comp1,Net);
equations Eq_tao_nrtl(N,j,comp,comp1,Net);
tao_nrtl.fx(N,j,comp,comp1,Net)$(ord(comp) eq ord(comp1))=0;
Eq_tao_nrtl(N,j,comp,comp1,Net).. tao_nrtl(N,j,comp,comp1,Net)=e=a_nrtl(comp,comp1) + (b_nrtl(comp,comp1)/Temp(N,j,Net));
variables g_nrtl(N,j,comp,comp1,Net);
equations Eq_g_nrtl(N,j,comp,comp1,Net);
g_nrtl.fx(N,j,comp,comp1,Net)$(ord(comp) eq ord(comp1))=1;
Eq_g_nrtl(N,j,comp,comp1,Net).. g_nrtl(N,j,comp,comp1,Net)=e=exp( -alfa_nrtl(comp,comp1)*tao_nrtl(N,j,comp,comp1,Net));

*Coeficiente de actividad (gamma)

alias (comp,comp2,comp3);
variables gamma(N,j,comp,Net);
equations Eqgamma(N,j,comp,Net);
Eqgamma(N,j,comp,Net).. gamma(N,j,comp,Net)=e=
        exp(sum(comp1,x(N,j,comp1,Net)*tao_nrtl(N,j,comp1,comp,Net)*
        g_nrtl(N,j,comp1,comp,Net))/sum(comp1,x(N,j,comp1,Net)*
        g_nrtl(N,j,comp1,comp,Net))+sum(comp1,x(N,j,comp1,Net)*
        g_nrtl(N,j,comp,comp1,Net)/sum(comp2,x(N,j,comp2,Net)*
        g_nrtl(N,j,comp2,comp1,Net))*(tao_nrtl(N,j,comp,comp1,Net)-
        sum(comp2,x(N,j,comp2,Net)*tao_nrtl(N,j,comp2,comp1,Net)*
        g_nrtl(N,j,comp2,comp1,Net))/sum(comp3,x(N,j,comp3,Net)*
        g_nrtl(N,j,comp3,comp1,Net)))));

*$offtext
*-------------------------------------------------------------------------------
*                                Secci?n 10
*                           Ecuaci?n de estado (calculo de phi)
*-------------------------------------------------------------------------------
parameter
Omega(comp) "Factor ac�ntrico [-]"
/
wat     0.344861
eth     0.643558
gly     0.51269
/

positive variable phi(N,j,comp,Net);
equation EqPhi(N,j,comp,Net);

EqPhi(N,j,comp,Net).. phi(N,j,comp,Net) =e= 0.985;

*-------------------------------------------------------------------------------
*                                Secci?n 11
*                           C?lculo de entalp?as
*-------------------------------------------------------------------------------
*Constantes de Cp (kJ/kmol.K) gas ideal
parameters
C1c(comp)
/
wat      35.86
eth      6.638866920
gly      12.02602320
/
C2c(comp)
/
wat      -0.0229
eth      .2249913420
gly      .4250839170
/
C3c(comp)
/
wat      6e-5
eth      -1.0887479E-4
gly      -2.8316157E-4
/
C4c(comp)
/
wat      -4e-8
eth      1.83461616E-8
gly      7.58362604E-8
/
;

parameter
Tref "Temperatura de referencia [K]" /298.15/
Hscale "Escalar de las entalp�as" /1e3/
;

*Entalp?a de la fase vapor (kJ/mol) -- Int(CpdT)
variable HVi(N,j,comp,Net),HV(N,j,Net);
equations EqHVi(N,j,comp,Net),EqHV(N,j,Net);
EqHVi(N,j,comp,Net).. HVi(N,j,comp,Net)=e=( (C1c(comp)*(Temp(N,j,Net)-Tref)) + ((C2c(comp)/2)*((Temp(N,j,Net)**2)-(Tref**2)))
                                   + ((C3c(comp)/3)*((Temp(N,j,Net)**3)-(Tref**3))) + ((C4c(comp)/4)*((Temp(N,j,Net)**4)-(Tref**4))))/Hscale;
HVi.l(N,j,comp,Net)=-400;

EqHV(N,j,Net).. HV(N,j,Net)=e=sum(comp,HVi(N,j,comp,Net)*y(N,j,comp,Net)/100);
HV.l(N,j,Net)=-400;

*Constantes de entalpia de vaporizacion (kJ/kmol)
parameters
C1v(comp)
/
wat      51546.00000
eth      55789.00000
gly      1.10670000e+5
/
C2v(comp)
/
wat      0.2840200000
eth      0.3124500000
gly      0.4831900000
/
C3v(comp)
/
wat      -0.1584300000
eth      0
gly      0
/
C4v(comp)
/
wat      0.2375000000
eth      0
gly      0
/
Tb(comp)
/
wat      100
eth      78.29
gly      287.85
/
;

*Temperaturas reducidas
variable Tred(comp);
equation EqTred(comp);
EqTred(comp)..Tred(comp)=e=Tb(comp)/Tcrit(comp);

*Entalp?a de vaporizaci?n (kJ/mol)
variable DHvap(N,j,comp,Net);
equation EqdHvap(N,j,comp,Net);

EqdHvap(N,j,comp,Net)..DHvap(N,j,comp,Net)=e=(C1v(Comp)*((1-(Temp(N,j,Net)/Tcrit(comp))))**(C2v(Comp)+
         C3v(Comp)*(Temp(N,j,Net)/Tcrit(comp))+(C4v(Comp)*(Temp(N,j,Net)/Tcrit(comp)))))/Hscale;

*Entalpia de la fase liquida (kJ/mol)

variable HLi(N,j,comp,Net),HL(N,j,Net);
equation EqHLi(N,j,comp,Net),EqHL(N,j,Net);
EqHLi(N,j,comp,Net).. HLi(N,j,comp,Net)=e=Hvi(N,j,comp,Net)-DHVap(N,j,comp,Net);

EqHL(N,j,Net).. HL(N,j,Net)=e=sum(comp,HLi(N,j,comp,Net)*x(N,j,comp,Net)/100);

*-------------------------------------------------------------------------------
*                                Secci?n 12
*                  C?lculo de entalp?a de alimentaci?n
*-------------------------------------------------------------------------------
*Mezcla aze�tropica
*Entalp?a de la alimentaci?n de la mezcla azeotr�pica

parameter HV_Az(comp)    "Entalp?a de vapor de la alimentaci?n [kJ/mol]"
          Tred_Az(comp)  "Temperatura reducida alimentaci?n [-]"
          DHVap_Az(comp) "Entalp?a de vaporizaci?n alimentaci?n [kJ/mol]"
          HL_Az(comp)    "Entalp?a de l?quido de la alimentaci?n [kJ/mol]";

HV_Az(comp)=( (C1c(comp)*(Tali2-Tref)) + ((C2c(comp)/2)*((Tali2**2)-(Tref**2))) + ((C3c(comp)/3)*((Tali2**3)-(Tref**3))) + ((C4c(comp)/4)*((Tali2**4)-(Tref**4))))/Hscale;

Tred_Az(comp)=Tali2/Tcrit(comp);
DHVap_Az(comp)=( C1v(comp)*( (1-Tred_Az(comp))**( C2v(comp) + (C3v(comp)*Tred_Az(comp)) + (C4v(comp)*(Tred_Az(comp)**2)) ) ) )/Hscale;
HL_Az(comp)=HV_Az(comp)-DHVap_Az(comp);

*HFAz se calcula para todas las etapas internas
variable HFAz(N,j,Net) "Entalpia de la alimentacion de aze�tropo";
equation EqHFAz(N,j,Net);
EqHFAz(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. HFAz(N,j,Net) =e= sum(comp,(zAz(N,j,comp)/100)*(HL_Az(comp)));

*GLICEROL
*Entalp?a de la alimentaci?n de glicerol

parameter HV_g(comp)    "Entalp?a de vapor de la alimentaci?n [kJ/mol]"
          Tred_g(comp)  "Temperatura reducida alimentaci?n [K]"
          DHVap_g(comp) "Entalp?a de vaporizaci?n alimentaci?n [kJ/mol]"
          HL_g(comp)    "Entalp?a de l?quido de la alimentaci?n [kJ/mol]";

HV_g(comp)=( (C1c(comp)*(Tali1-Tref)) + ((C2c(comp)/2)*((Tali1**2)-(Tref**2))) + ((C3c(comp)/3)*((Tali1**3)-(Tref**3))) + ((C4c(comp)/4)*((Tali1**4)-(Tref**4))))/Hscale;

Tred_g(comp)=Tali1/Tcrit(comp);
DHVap_g(comp)=( C1v(comp)*( (1-Tred_g(comp))**( C2v(comp) + (C3v(comp)*Tred_g(comp)) + (C4v(comp)*(Tred_g(comp)**2)) ) ) )/Hscale;
HL_g(comp)=HV_g(comp)-DHVap_g(comp);

*HFG se alcula para todas las etapas internas
variable  HFG(N,j,Net)   "Entalp?a de la alimentaci?n de glicerol [kJ/mol]";
equation EqHFG(N,j,Net);
EqHFG(N,j,Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1)).. HFG(N,j,Net) =e= sum(comp,(zg(comp)/100)*(HL_g(comp)));

*-------------------------------------------------------------------------------
*                                Secci?n 13
*                Definicion de parametros binarios
*-------------------------------------------------------------------------------

*Existencia de reflujo
parameter yr(Net) "1 indica que en la etapa si hay reflujo";
yr(Net)=0;
yr('2')=1;

*Existencia de boil up
parameter yb(Net) "1 indica que en la etapa si hay boil up";
yb(Net)=0;

*Permite saber si la etapa es real o sobrante
parameter par(Net) "1 indica que la etapa es real fisicamente";

*Existencia de alimentacion
parameter yfG(Net) "1 indica que en la etapa hay alimentacion de Glicerol";
yfG(Net)=0;

parameter yfAz(Net) "1 indica que en la etapa hay alimentacion de Mezcla Az.";
yfAz(Net)=0;

*-------------------------------------------------------------------------------
*                                Secci?n 14
*                           Restricciones logicas-Realizarla caso hipot�tico
*-------------------------------------------------------------------------------
*$ontext

equation logic1(Net) "The boil up stage is below the reflux stage";
logic1(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=(yb(Net));

equation logic2 "There is one reflux stage";
logic2..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yr(Net1)))=e=1;

equation logic3 "There is one boil up stage";
logic3..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yb(Net1)))=e=1;

equation logic4 "There is one feed stage of Glycerol";
logic4..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yfG(Net1)))=e=1;

equation logic5 "There is one feed stage of Mezcla azeotr�pica";
logic5..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le (card(Net1)-1))),yfAz(Net1)))=e=1;

equation logic6(Net) "Glycerol feed stage is below the reflux";
logic6(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=yfG(Net);

equation logic7(Net) "Mezcla azeotr�pica feed stage is below the reflux";
logic7(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=yfAz(Net);

equation logic8(Net) "The boil up stage is below the glycerol feed stage";
logic8(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yfG(Net1)))=g=yb(Net);

equation logic9(Net) "The boil up stage is below the Mezcla azeotr�pica feed stage";
logic9(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yfAz(Net1)))=g=yb(Net);

equation logic10(Net)  "The Glycerol feed is above the Azeotropic feed";
logic10(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yfG(Net1)))=g=yfAz(Net);

equation logic11(Net) "Glycerol feed stage is below the reflux";
logic11(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=yfG(Net);

equation logic12(Net) "Azeotropic feed stage is below the reflux";
logic12(Net)$(ord(Net)>1 and ord(Net)<card(Net))..(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))=g=yfAz(Net);
*$offtext

*-------------------------------------------------------------------------------
*                                Secci?n 15
*                         Ecuaciones del condensador
*-------------------------------------------------------------------------------

*Condiciones iniciales (operaci?n en estado estable)
equation BalMasaC0,BalMasaParcialC0(comp),SumaC0,EquilibrioC0(comp),BalEnergiaC0;
BalMasaC0.. 0=e=V('1','1','2')-V('1','1','1')*(1+RR('1','1'));
BalMasaParcialC0(comp).. 0=e=V('1','1','2')*y('1','1',comp,'2')-V('1','1','1')*x('1','1',comp,'1')*(1+RR('1','1'));
SumaC0.. sum(comp,y('1','1',comp,'1')-x('1','1',comp,'1'))=e=0;
EquilibrioC0(comp).. y('1','1',comp,'1')*P('1','1','1')*phi('1','1',comp,'1')=e=Psat('1','1',comp,'1')*gamma('1','1',comp,'1')*x('1','1',comp,'1');
BalEnergiaC0.. 0=e=V('1','1','2')*HV('1','1','2')-V('1','1','1')*(1+RR('1','1'))*HL('1','1','1')-QC('1','1');

*Flujo de liquido fijo
equation fixedL(N,j);
fixedL(N,j)..L(N,j,'1')=e=0;
*-------------------------------------------------------------------------------
*                                Secci?n 16
*                    Ecuaciones de la columna (Punto Inicial)
*-------------------------------------------------------------------------------

*Condiciones iniciales (operacion en estado estable)
equations BalMasa0(Net,Net1),BalMasaParcial0(comp,Net,Net1),Suma0(Net),BalEnergia0(Net,Net1);

BalMasa0(Net,Net1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(Net1) eq card(Net1)))..0=e=yfG(Net)*FG*1000 + yfAz(Net)*FAz*1000 + RR('1','1')*V('1','1','1')*yr(Net) +BR('1','1')*L('1','1',Net1)*yb(Net)+ L('1','1',Net-1) +V('1','1',Net+1) -L('1','1',Net) -V('1','1',Net) ;
Suma0(Net)$(ord(Net)>1 and ord(Net)<card(Net)).. sum(comp,x('1','1',comp,Net)-y('1','1',comp,Net))=e=0;
BalEnergia0(Net,Net1)$(ord(Net)>1  and ord(Net)<card(Net) and ord(Net1) eq card(Net1))..0=e=yfG(Net)*FG*HFG('1','1',Net)*1000 + yfAz(Net)*FAz*HFAz('1','1',Net)*1000 + RR('1','1')*V('1','1','1')*yr(Net)*HL('1','1','1') + BR('1','1')*L('1','1',Net1)*yb(Net)*HV('1','1',Net1) +L ('1','1',Net-1)*HL('1','1',Net-1) + V('1','1',Net+1)*HV('1','1',Net+1) - L('1','1',Net)*HL('1','1',Net) - V('1','1',Net)*HV('1','1',Net);
BalMasaParcial0(comp,Net,Net1)$((ord(Net)>1 and ord(Net)<card(Net)) and (ord(Net1) eq card(Net1)))..0=e=yfG(Net)*FG*zg(comp)*1000 + yfAz(Net)*FAz*zAz('1','1',comp)*1000 + RR('1','1')*V('1','1','1')*yr(Net)*x('1','1',comp,'1') + BR('1','1')*L('1','1',Net1)*yb(Net)*y('1','1',comp,Net1) + L('1','1',Net-1)*x('1','1',comp,Net-1) + V('1','1',Net+1)*y('1','1',comp,Net+1) - L('1','1',Net)*x('1','1',comp,Net) - V('1','1',Net)*y('1','1',comp,Net);

*Relaciones de equilibrio
equations Equilibrio10(comp,Net);
Equilibrio10(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net))..0=e=par(net)*((y('1','1',comp,Net)*P('1','1',Net)*phi('1','1',comp,'1'))-
    (Psat('1','1',comp,Net)*gamma('1','1',comp,Net)*x('1','1',comp,Net)));

equations Equilibrio20(comp,Net);
Equilibrio20(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net) and ord(comp) ne 1)..0=e=(sum(Net1$(ord(Net1) ge 2 and ord(Net1) le ord(Net)),yr(Net1)))*(1-par(Net))*(y('1','1',comp,Net)-y('1','1',comp,Net+1));
equation Equilibrio30(comp,Net);
Equilibrio30(comp,Net)$(ord(Net)>1 and ord(Net)<card(Net) and ord(comp) ne 1)..0=e=
    (1-sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))*(1-par(Net))*(x('1','1',comp,Net)-x('1','1',comp,Net-1));
Equation Equilibrio40(Net,Net1);
Equilibrio40(Net,Net1)$(ord(Net)>1 and ord(Net)<card(Net) and ord(Net1) eq card(Net1))..0=e=
    (1-par(Net))*(V('1','1',Net)-V('1','1',Net+1)-BR('1','1')*L('1','1',Net1)*yb(Net));

*-------------------------------------------------------------------------------
*                                Secci?n 17
*                        Ecuaciones del rehervidor
*-------------------------------------------------------------------------------

*Condiciones iniciales (operaci?n en estado estable)
equation BalMasaR0(Net),BalMasaParcialR0(comp,Net),SumaR0(Net),EquilibrioR0(comp,Net),BalEnergiaR0(Net);

BalMasaR0(Net)$(ord(Net) eq card(Net))..0=e=L('1','1',Net-1)-L('1','1',Net)*(1+BR('1','1'));
BalMasaParcialR0(comp,Net)$(ord(Net) eq card(Net))..0=e=L('1','1',Net-1)*x('1','1',comp,Net-1)-L('1','1',Net)*(x('1','1',comp,Net)+BR('1','1')*y('1','1',comp,Net));
SumaR0(Net)$(ord(Net) eq card(Net)).. sum(comp,y('1','1',comp,Net)-x('1','1',comp,Net))=e=0;
EquilibrioR0(comp,Net)$(ord(Net) eq card(Net))..y('1','1',comp,Net)*P('1','1',Net)*phi('1','1',comp,'1')=e=Psat('1','1',comp,Net)*gamma('1','1',comp,Net)*x('1','1',comp,Net);
BalEnergiaR0(Net)$(ord(Net) eq card(Net))..0=e=QR('1','1')+L('1','1',Net-1)*HL('1','1',Net-1)-L('1','1',Net)*HL('1','1',Net)-BR('1','1')*L('1','1',Net)*HV('1','1',Net);

*Variable fija del flujo de vapor en la ultima etapa
equation fixedV(N,j,Net);
fixedV(N,j,Net)$(ord(Net) eq card(Net))..V(N,j,Net)=e=0;

*-------------------------------------------------------------------------------
*                                Secci?n 18
*               Relaciones hidr?ulicas para todas las etapas internas
*-------------------------------------------------------------------------------


*Definici?n de velocidad de vapor
positive variables far(N,j,Net) "Factor de areaciin [-]";
equations Eqfa(N,j,Net);
Eqfa(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(far(N,j,Net))=e=par(Net)*(0.981*exp(-0.411*(Qvap(N,j,Net)/3600)/At*(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)**(0.5)));
far.l(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))=0.35;

positive variable hD(N,j,Net)   "Altura del liquido por encima del divisor [m]";
equations EqhD(N,j,Net);
EqhD(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (hD(N,j,Net))=e=0.6*((Qliq(N,j,Net)/3600/Lw)**(2/3));
hD.l(N,j,Net)=1e-2;

positive variable uhv(N,j,Net) "Velocidad del vapor por los agujeros [m/s]";
equations Equhv(N,j,Net);
Equhv(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(uhv(N,j,Net)*A0)=e=par(Net)*((V(N,j,Net)/(rhoV(N,j,Net))))/3600;

positive variable unv(N,j,Net) "Velocidad del vapor por el plato [m/s]";
equations Equnv(N,j,Net);
Equnv(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*unv(N,j,Net)*At=e=par(Net)*((V(N,j,Net)/(rhoV(N,j,Net))))/3600;

*Definicion de velocidad del liquido
positive variable ul(N,j,Net) "Velocidad del l?quido en el derramadero [m/s]";
equations Equl(N,j,Net);
Equl(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*ul(N,j,Net)*Ad=e=par(Net)*((L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)))/3600;

*Carga de l�quido
positive variable hcl(N,j,Net)  "Altura del l?quido libre en r?gimen de spray [m]"
equation Eqhcl(N,j,Net);
scalar consmach /1e-20/;
Eqhcl(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*hcl(N,j,Net)=e=par(Net)*((0.157*((poro+consmach)**(-0.791))/(1+1.04E-4*(((((L(N,j,Net)+consmach)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)))/3600/Lw)**(-0.59))
                                                        *((poro+consmach)**(-1.791))))*(d_hole**0.833)
                                                        *(996/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))**(0.5*(1-0.91*d_hole/(poro+consmach))));

positive variable Csbf(N,j,Net);
equation EqCsbf(N,j,Net);
*EqCsbf(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(Csbf(N,j,Net))=e=par(Net)*0.15;
EqCsbf(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. par(Net)*(Csbf(N,j,Net))=e=par(Net)*(0.37*(((sqr(d_hole)*sigma(N,j,Net)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)))**0.125)
                                                        *((((rhoV(N,j,Net))*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)/(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))**0.1)
                                                        *((HS/hcl(N,j,Net))**0.5));

*Caida de presi�n
positive variables DPL(N,j,Net) "Ca?da de presi?n por la presencia de l?quido [atm]";
equations EqDPL(N,j,Net);
EqDPL(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DPL(N,j,Net))=e=((far(N,j,Net)*9.81*(sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)*(hD(N,j,Net)+hw))/101325);
DPL.l(N,j,Net)=9e-5;

positive variables DPS(N,j,Net) "Ca?da de presi?n debido a la presencia de los agujeros - seco [atm]";
equations EqDPS(N,j,Net);
EqDPS(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DPS(N,j,Net))=e=0.0001;

positive variable DPq(N,j,Net)  "Ca?da de presi?n en el derramadero [atm]";
equations EqDPq(N,j,Net);
EqDPq(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. DPq(N,j,Net)*101325=e=1.62*((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))/(sqr(Lw*hw))*(sqr((L(N,j,Net)/sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)/3600))+sqr((V(N,j,Net)/(rhoV(N,j,Net))/3600)));

positive variables DP(N,j,Net)  "Ca?da de presi?n total [atm]";
equations EqDP(N,j,Net),EqDPR(N,j,Net),EqP(N,j,Net),EqPC(N,j,Net),EqPR(N,j,Net) "Definici?n de presi?n por etapa [atm]";

EqDPR(N,j,Net)$(ord(Net) eq card(Net)).. DP(N,j,Net)=e=DP(N,j,Net-1);
EqDP(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. (DP(N,j,Net))=e=(DPS(N,j,Net)+DPL(N,j,Net));
EqP(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. P(N,j,Net)=e=P(N,j,Net-1)+par(Net)*DP(N,j,Net);
EqPC(N,j,Net)$(ord(Net) eq 1).. P(N,j,Net)=e=Pop;
EqPR(N,j,Net)$(ord(Net) eq card(Net)).. P(N,j,Net)=e=P(N,j,Net-1);



*Efectos indeseados en la columna
*Downflow flooding (inundaci?n en los derramaderos)
variable downf(N,j,Net);
equation DownFlood(N,j,Net)
;
DownFlood(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..downF(N,j,Net)*par(Net)=e=((HD(N,j,Net)+((DP(N,j,Net)*101325+DPq(N,j,Net)))
                                                        /(9.81*(((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))))-(HS))*par(Net);


*Entrainment flooding (inundaci?n por arrastre de l?quido)
variables
EntrainFloodLVar(N,j,Net)
EntrainFloodVVar(N,j,Net)
;
equations
EntrainFloodL(N,j,Net)
EntrainFloodV(N,j,Net)
;

EntrainFloodV(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. EntrainFloodVVar(N,j,Net)*par(Net) =e= par(Net)*((unv(N,j,Net))-(Csbf(N,j,Net)*(((((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                       -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)))
                                                        /(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))**0.5));

EntrainFloodL(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net))..EntrainFloodLVar(N,j,Net)*par(Net) =e= par(Net)*((ul(N,j,Net))-((sigma(N,j,Net)*9.81*(((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000))
                                                        -(rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000))
                                                        /((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)**2))**(1/4)));

*Weeping (lloriqueo)
equation Weep(N,j,Net);
Weep(N,j,Net)$(ord(Net)>1 and ord(Net)<card(Net)).. 0=g=(((0.68+0.12)/(((rhoV(N,j,Net)*sum(comp,MW(comp)*y(N,j,comp,Net)/100)/1000)
                                                                /((sum(comp,rho(N,j,comp,Net)*x(N,j,comp,Net)/100)*sum(comp,MW(comp)*x(N,j,comp,Net)/100)/1000)
                                                                *9.81*far(N,j,Net)*(hw+hd(N,j,Net))))**0.5))-(uhv(N,j,Net)))*par(Net);

*Construcci?n de la columna
equation Size "Tama�o del equipo";
Size.. Htotal =e= ((1+Sfactor)*sum(Net$(ord(Net)>1 and ord(Net)<card(Net)),HS*par(Net)));

equation DtoLratio "Relacion entre el diametro y la altura";
DtoLratio..htotal/D=l=20;

*===============================================================================
*                                Secci�n 19
*                     C�lculo de costos de infraestructura,
*                        Funci�n objetivo y soluci�n
*===============================================================================

scalars
CR_u   "Precio energ�a para el rehervidor $/kJ/h" /1e-5/
CC_u    "Precio energ�a para el condensador $/kJ/h"
CEt   "Precio Etanol $/mol" /0.03/
;
*CR_u = 146.8/3600;
CC_u = 24.5/3600/6000;
*===============================================================================
*Definici�n de los costos de infraestructura
*===============================================================================
*$ontext
scalars
*Indices CEPCI
I1 /397/
I2 /575.4/

*Par�metros Aspen Plus
sigma_d "Tensi�n superficial promedio del l�quido (dyn/cm)" /50/
sigma0 "Tensi�n superficial de referencia (dyn/cm)" /20/
*Valores de A. Energy Analyzer para los Intercambiadores
*Condensador
Ucond "Coeficiente de transferencia (kJ/hr K m2)"
lmtd_cond "LMTD promedio" /35/
*Rehervidor
Ureb "Coeficiente de transferencia (kJ/hr K m2)"
lmtd_reb "LMTD promedio" /35/
*Par�metros de Luyben para dimensionamiento de la columna
c0_fair "Par�metro 0 para la correlaci�n de Fair (ts=0.61m)" /439/
c1_fair "Par�metro 1 para la correlaci�n de Fair (ts=0.61m)" /2.5/
c2_fair "Par�metro 2 para la correlaci�n de Fair (ts=0.61m)" /1.2/
an_ova "An/A" /0.8/
ophour "Horas de operaci�n al a�o" /6000/
anfact "Factor anualizante" /0.25/
;

Ucond=0.95*3600;
Ureb=0.65*3600;

*-------------------------------------------------------------------------------
*Dimensionamiento de los intercambiadores
*-------------------------------------------------------------------------------
positive variables
Area_cond "�rea del condensador (m2)"
Area_reb  "�rea del rehervidor (m2)"
;
equations
def_a_cond "Definici�n del �rea del condensador (m2)"
def_a_reb  "Definici�n del �rea del rehervidor (m2)"
;

def_a_cond.. Area_cond*Ucond*lmtd_cond=e=1.25*Qc('1','1');
Area_cond.l=39.82;

def_a_reb.. Area_reb*Ureb*lmtd_reb=e=1.25*Qr('1','1');
Area_reb.l=73.94;

*-------------------------------------------------------------------------------
*Dimensionamiento de la columa
*-------------------------------------------------------------------------------
positive variable
dliqprom "Densidad promedio del l�quido (kmol/m3)"
dvapprom "Densidad promedio del vapor (kmol/m3)"
mwpromliq "MW promedio del l�quido"
mwpromvap "MW promedio del vapor"
FP_prom "FP promedio"
V_prom "Flujo de vapor promedio (kmol/hr)"
cap_par "Par�metro de capacidad"
f_flood "Par�metro de inundaci�n"
;

equation
dliqprom_def "Definici�n de la Densidad promedio del l�quido (kmol/m3)"
dvapprom_def "Definici�n de la Densidad promedio del vapor (kmol/m3)"
mwpromliq_def "Definici�n del MW promedio del l�quido (kg/mol)"
mwpromvap_def "Definici�n del MW promedio del vapor (kg/mol)"
*EqFP(N,j,Net)
FP_prom_def "Definici�n del FP promedio"
V_prom_def "Definici�n del flujo de vapor promedio"
cap_par_def "Definici�n del par�metro de capacidad"
f_flood_def "Definici�n del par�metro de inundaci�n"
A_col_def "Definici�n del �rea promedio de la columna"
;

dliqprom_def.. dliqprom=e=sum((Net,comp),par(Net)*rho('1','1',comp,Net)*x('1','1',comp,Net))/sum(Net,par(Net))/100;

dvapprom_def.. dvapprom=e=sum((Net),par(Net)*rhoV('1','1',Net))/sum(Net,par(Net));

mwpromliq_def.. mwpromliq=e=(sum((comp,Net),par(Net)*MW(comp)*x('1','1',comp,Net)/100)/1000)/sum(Net,par(Net));
mwpromliq.l=60/1000;

mwpromvap_def.. mwpromvap=e=(sum((comp,Net),par(Net)*MW(comp)*y('1','1',comp,Net)/100)/1000)/sum(Net,par(Net));
mwpromvap.l=30/1000;

V_prom_def.. v_prom=e=sum((Net),par(Net)*V('1','1',Net))/sum(Net,par(Net));

*FP_prom_def.. FP_prom=e=sum((Net)$(card(net) ge 2 and ord(net) le (card(net)-1)), L('1','1',Net)/V('1','1',Net))/sum(Net,par(Net))*mwpromliq/mwpromvap*((dvapprom/dliqprom*mwpromliq/mwpromvap)**(0.5));
fp_prom_def.. FP_prom=e=0.07;

cap_par_def.. cap_par*(1+c1_fair*(fp_prom**c2_fair))=e=c0_fair;
cap_par.l = 410.7;

scalar scaleflood "Escalar del factor de inundaci�n"  /1e3/
;

f_flood_def.. sqr(f_flood*scaleflood)=e=sqr(cap_par*an_ova*(sigma0/sigma_d)**0.2)*
         (dliqprom*mwpromliq-dvapprom*mwpromvap);

A_col_def.. sqr(A_col*0.6*f_flood*scaleflood)*dvapprom*mwpromvap
         =e=sqr(mwpromvap*V_prom);
A_col.l=0.373;
Area_cond.l=0.373;

*-------------------------------------------------------------------------------
*Costos de construcci�n: Turton
*-------------------------------------------------------------------------------
positive variables
cost_cond "Costo del condensador (1e5$)"
cost_reb  "Costo del rehervidor (1e5$)"
cost_col "Costo de la columna (1e5$)"
cost_sta "Costo de las etapas (1e5$)"
tot_inf "Costo total de infraestructura (1e5$)"
;

equations
cost_cond_def "Definici�n del costo del condensador"
cost_reb_def  "Definici�n del costo del rehervidor"
cost_col_def  "Definici�n del costo de la columna"
cost_sta_def "Definici�n del costo de las etapas"
tot_inf_def "Definici�n del costo total de infraestructura"
;

scalar
scale_cost "Escalar de los costos de infraestructura" /1e5/
;

cost_cond_def.. cost_cond=e=(I2/I1)*10**(3.7803+(0.8569*log10(Area_cond))+0.0349
         *sqr(log10(Area_cond)));
cost_cond.l=20523;

cost_reb_def.. cost_reb=e=(I2/I1)*10**(4.4646-(0.5277*log10(Area_reb))+0.3955*
         sqr(log10(Area_reb)));
cost_reb.l=10000;

cost_col_def.. cost_col=e=(I2/I1)*10**(3.4974+(0.4485*log10(htotal*A_col))
         +0.1074*sqr(log10(htotal*A_col)));
cost_col.l = 20000;

cost_sta_def.. cost_sta=e=(I2/I1)*(sum(Net,par(Net))-2)*10**(2.9949+(0.4465
         *log10(A_col))+0.3961*sqr(log10(A_col)));
cost_sta.l = 10000;

variable tot_inf "Costo total de infraestructura";
equation tot_inf_def "Definici�n del costo total de infraestructura";
tot_inf_def.. tot_inf=e=(cost_sta+cost_col+cost_cond+cost_reb)*anfact;
*-------------------------------------------------------------------------------
*Funci�n objetivo
*-------------------------------------------------------------------------------
variables
CostMP "Costo materia prima"
CostQr "Costo de operaci�n del rehervidor"
CostQc "Costo de operaci�n del condensador"
GanEth "Ganancia de etanol"
;
Equations EqCostMP "Costo materia prima", EqCostQr, EqCostQc, EqGanEth;
EqCostMP..CostMP =e= Cet/2*FAz*1000*ophour;
EqCostQr..CostQr =e= Qr('1','1')*CR_u*ophour;
EqCostQc..CostQc =e= Qc('1','1')*CC_u*ophour;
EqGanEth..GanEth =e= Cet*V('1','1','1')*ophour;

variables
zobj "Funci�n objetivo";
equations
Fobj "Funci�n objetivo";
*FObj..zobj=e=1;
FObj.. zobj=e=(CostMP + CostQr + CostQc - GanEth) +(tot_inf);

*-------------------------------------------------------------------------------
*                                Secci?n 20
*                            Cotas en las variables
*-------------------------------------------------------------------------------

Temp.lo(N,j,Net)= 300;
Temp.up(N,j,Net)= 510;

Tcritm.lo(N,j,Net)=510;
Tcritm.up(N,j,Net)=860;

x.lo(N,j,comp,Net)= 0;
x.up(N,j,comp,Net)= 100;

y.lo(N,j,comp,Net)= 0;
y.up(N,j,comp,Net)= 100;

rho.up(N,j,comp,Net)=70000;
rho.lo(N,j,comp,Net)=1000;

D.up=10;
D.lo=0.01;

Qr.up(N,j)=1e8;
Qr.lo(N,j)=1e5;
Qr.l(N,j)=6e6;

Qc.up(N,j)=2e8;
Qc.lo(N,j)=6e4;
Qc.l(N,j)=4e6;

RR.lo(N,j)=1e-4;
RR.up(N,j)=10;
RR.l(N,j)=0.01;

BR.up(N,j)=20;
Br.L(N,j)=2;

L.up(N,j,Net)=1e6;
V.up(N,j,Net)=1e6;

P.up(N,j,Net)= 2;
P.lo(N,j,Net)=Pop;

htotal.up=30;
htotal.lo=6;
htotal.l=10;

*execute_loadpoint "init.gdx";

*------INITIAL VALUE OF BINARY TERMS--------------------------------------------
*Existencia de reflujo

yfG('2')=1;
yfAz('11')=1;
yb('18')=1;

par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));


*----------------------MODEL AND SOLUTION---------------------------------------
*model Extractive_Distillation /all/;
model DSDA /all/;
option nlp=BARON;
*option nlp= IPOPT;
*option nlp= CONOPT;
*option nlp = CONOPT4;

$embeddedCode Python:
import csv
import os
###User params 1: independent & dependent binar terms---------------------------
Indep={"x1":["yb",None,"Net"],
"x2":["yfG",None,"Net"],
"x3":["yfAz",None,"Net"]}

Dep=["par(Net)$((ord(Net) ne card(Net)) and (ord(Net) ne 1))=(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yr(Net1)))+(yb(Net))-(sum(Net1$((ord(Net1) ge 2) and (ord(Net1) le ord(Net))),yb(Net1)));"]

###User params 2: neighborhood--------------------------------------------------
#use: Infinity,Separable,Mflat,Lflat
neighborhood="Infinity"

###User params 3:User defined reformulation,initialization----------------------
###and_ inequality constraints--------------------------------------------------


#use: Automatic,User_defined
Reform_type="Automatic"
Reform=[""]
Initial=[""]

#use: Automatic,User_defined
Ineq_type="Automatic"
#merge_constraints: 1 if_ user constraints must be merged with automatic constraints.
#merge_constraints: 0 if_ user defined constraints are going to be used alone.
merge_constraints="1"
Ineq=[""]

###user params 3: Decide if_ logical constranits are going to remain in the subProblems
#logic_sub:1 if_ logical constraints appear in the subproblems, gams is allowed to skip the solution stage for_ infeasible neighbors
#logic_sub:0 if_ logical constraints do_ not_ appear in the subproblems, the master problem oversees the logic_
logic_sub="1"

###Create list of external variabs (extvar)-------------------------------------
extvar=list()
for i in Indep:
 Indep[i][1]=str(sum(list(gams.get(Indep[i][0],KeyType.INT,KeyFormat.SKIP))))
 numb=int(float(Indep[i][1]))
 for j in range(1,numb+1):
  extvar.append(i+"_"+str(j))
  Indep[i].append(i+"_"+str(j))

#Total numb of external var
numb=len(extvar)


#print("number of external variables=",numb)
#print("set of external variables=",extvar)


#Total numb of neighbors
if neighborhood=="Infinity":
 num_neigh=(3**numb)-1
elif neighborhood=="Separable":
 num_neigh=2*numb
elif neighborhood=="Mflat":
 num_neigh=numb*(numb+1)
elif neighborhood=="Lflat":
 num_neigh=2**(numb+1)-2

#extvar set_ and_ subsets
Ext="set extvar /"+','.join(extvar)+"/; \n"
Ext_sub=""
for i in Indep:
 Ext_sub=Ext_sub+"set s_"+i+"(extvar) /"+','.join(Indep[i][3:None])+"/; \n"
#print(Ext)
#print(Ext_sub)

neigh="set neigh /d1*d"+str(num_neigh)+"/; \n"
#print(neigh)

#create de csv file_ with the niegihbors
os.chdir('./neighborhood_MATLAB')
file=str(numb)+'_'+neighborhood+'.csv'
with open(file) as csvfile:
 readCSV = csv.reader(csvfile, delimiter=',')
 neigh_values=","+','.join(extvar)+"\n"
 for row in readCSV:
  neigh_values=neigh_values+str(','.join(row))+"\n"
os.chdir('..')

#print(neigh_values)

#get and_ add ordered sets_in list format
for i in Indep:
 Indep[i].append(list(gams.get(Indep[i][2],KeyType.STRING)))

#create aux variables_ for_ reformulation
Ext_aux=""
for i in Indep:
 Ext_aux=Ext_aux+"parameter yaux_"+i+"(extvar,"+Indep[i][2]+");\nyaux_"+i+"(extvar,"+Indep[i][2]+")=0;\n"

#get and_ add value of intependente binary_ terms for_ initializations
for i in Indep:
 Indep[i].append(list(gams.get(Indep[i][0],KeyType.STRING)))

#print("\n",str(Indep))


#Creation of external variables_ functions_ and_ initialization
reform_string=""
initialization_string=""
if Reform_type=="User_defined":
 reform_string='\n'.join(Reform)
 initialization_string='\n'.join(Initial)
if Reform_type=="Automatic":
 cuent=0;
 for i in Indep:
  reform_string=reform_string+"fex(s_"+i+")=xvalue_p(s_"+i+");\n"
  for j in Indep[i][-1]:
   initialization_string=initialization_string+"xvalinit(\'"+extvar[cuent]+"\')="+str(j[0])+";\n"
   cuent=cuent+1

#print(reform_string)
#print(initialization_string)
#Total number of inequality constraints and_ definiton of ineqs
number_user_defined_ineq=len(Ineq)
numb_ineq=numb*2
ineq_string_Automatic=""
ineq_string_User=""
ineq_string=""
cuent1=1;
for i in Indep:
 cuent2=0;
 for j in Indep[i][-1]:
  cuent2=cuent2+1;
  ineq_string_Automatic=ineq_string_Automatic+"gval(\'g"+str(cuent1)+"\')=1-fex(\'"+i+"_"+str(cuent2)+"\');\n"
  ineq_string_Automatic=ineq_string_Automatic+"gval(\'g"+str(cuent1+1)+"\')=fex(\'"+i+"_"+str(cuent2)+"\')-"+str(len(Indep[i][-2]))+";\n"
  cuent1=cuent1+2
  if cuent2>=2:
   ineq_string_Automatic=ineq_string_Automatic+"gval(\'g"+str(cuent1)+"\')=1+fex(\'"+i+"_"+str(cuent2-1)+"\')-fex(\'"+i+"_"+str(cuent2)+"\');\n"
   numb_ineq=numb_ineq+1
   cuent1=cuent1+1

if Ineq_type=="Automatic":
 ineq_string=ineq_string_Automatic
elif Ineq_type=="User_defined" and merge_constraints=="1":
 numb_ineq=numb_ineq+int(number_user_defined_ineq)
 for i in Ineq:
  ineq_string_User=ineq_string_User+"gval(\'g"+str(cuent1)+"\')="+i
  cuent1=cuent1+1
 ineq_string=ineq_string_Automatic+ineq_string_User
elif Ineq_type=="User_defined" and merge_constraints=="0":
 numb_ineq=int(number_user_defined_ineq)
 cuent1=1;
 for i in Ineq:
  ineq_string_User=ineq_string_User+"gval(\'g"+str(cuent1)+"\')="+i
  cuent1=cuent1+1
 ineq_string=ineq_string_User

Ineq_set="set extineq /g1*g"+str(numb_ineq)+"/; \n"
#print(Ineq_set)
#print(ineq_string)

#Reformulation with external variables_
recalculation_string="\n"
for i in Indep:
 recalculation_string=recalculation_string+"yaux_"+i+"(extvar,"+Indep[i][2]+")=0;\n"
 recalculation_string=recalculation_string+"yaux_"+i+"(extvar,"+Indep[i][2]+")$(ord("+Indep[i][2]+") eq floor(fex(extvar)))=1-mod(fex(extvar),1);\n"
 recalculation_string=recalculation_string+"yaux_"+i+"(extvar,"+Indep[i][2]+")$(ord("+Indep[i][2]+") eq 1+floor(fex(extvar)))=mod(fex(extvar),1);\n\n"
 recalculation_string=recalculation_string+Indep[i][0]+"("+Indep[i][2]+")=sum(s_"+i+",yaux_"+i+"(s_"+i+","+Indep[i][2]+"));\n"
recalculation_string=recalculation_string+'\n'.join(Dep)

#print(recalculation_string)

#treatment of logical constraints
logic_string1=""
logic_string2=""

if logic_sub=="1":
 logic_string1=logic_string1+"maxExecError=100000;"
 logic_string2=logic_string2+"ExecError=0;"
###Writing txts-----------------------------------------------------------------
file = open("neigh_values.csv",'w')
file.write(neigh_values)
file.close()

file=open("param_sets.txt",'w')
file.write(Ext+Ext_sub+neigh+Ineq_set)
file.close()

file=open("param_extaux.txt",'w')
file.write(Ext_aux+initialization_string)
file.close()

file=open("recalculation.txt",'w')
file.write(recalculation_string)
file.close()

file=open("reformulation.txt",'w')
file.write(reform_string)
file.close()

file=open("Inequalities.txt",'w')
file.write(ineq_string)
file.close()

file=open("logics1.txt",'w')
file.write(logic_string1)
file.close()

file=open("logics2.txt",'w')
file.write(logic_string2)
file.close()
$endEmbeddedCode
*$offtext
