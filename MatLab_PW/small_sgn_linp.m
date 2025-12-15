function [A, B, pos_var, state_var]=small_sgn_linp(Rete,gendatname,avridx,pssidx)
%% Copyright Valentin Ilea
%for the moment
f=50;
[baseMVA,bus,gen,branch,area,gencost,nomnod,linee, nomare, lim_area, cpiu, cmeno, capa]=feval(Rete); 
out=feval(gendatname); gendat=out{1}; 

%% LOAD FLOW
warning('off','MATLAB:dispatcher:InexactMatch');
path(path,'Load_flow');
gauss=1; limitato=1; display=1;

[bus,gen,success,Yf,Yt,Ybus,V]=load_flow(baseMVA,bus,gen,branch,nomnod,gauss,limitato,display,capa);

%% READING GENERATORS DYYNAMIC DATA

% rearanging genidx depending on model

genidx=gen(:,1); model=gendat(:,19); loadidx=bus(:,1); loadidx(genidx)=[];
genidx2=diag(model==2)*genidx; genidx4=diag(model==4)*genidx; genidx51=diag(model==51)*genidx; genidx52=diag(model==52)*genidx;

km=1;
for m=1:4
    if m==1; for k=1:size(genidx,1); if genidx2(k)~=0; genidx(km)=genidx2(k); genidxp(km)=k; km=km+1; end; end; end
    if m==2; for k=1:size(genidx,1); if genidx4(k)~=0; genidx(km)=genidx4(k); genidxp(km)=k; km=km+1; end; end; end
    if m==3; for k=1:size(genidx,1); if genidx51(k)~=0; genidx(km)=genidx51(k); genidxp(km)=k; km=km+1; end; end; end
    if m==4; for k=1:size(genidx,1); if genidx52(k)~=0; genidx(km)=genidx52(k); genidxp(km)=k; km=km+1; end; end; end
end

% number of generators of certain model
True=(model==2);  a=trace(diag(True)); %model 2  generators no
True=(model==4);  b=trace(diag(True)); %model 4  generators no
True=(model==51); c=trace(diag(True)); %model 51 generators no
True=(model==52); d=trace(diag(True)); %model 52 generators no

% reading data
% generators:
H=gendat(genidxp,2); D=gendat(genidxp,3); 
Tpd=gendat(genidxp,14); Tpq=gendat(genidxp,15); Tppd=gendat(genidxp,16);
Tppq=gendat(genidxp,17); TAA=gendat(genidxp,18);
Pg=gen(genidxp,2)/baseMVA; Qg=gen(genidxp,3)/baseMVA; Sb=gen(genidxp,7); genstat=gen(genidxp,8); 
xd=gendat(genidxp,8).*baseMVA./gen(genidxp,7); xq=gendat(genidxp,9).*baseMVA./gen(genidxp,7);
xpd=gendat(genidxp,10).*baseMVA./gen(genidxp,7); xpq=gendat(genidxp,11).*baseMVA./gen(genidxp,7);  
xppd=gendat(genidxp,12).*baseMVA./gen(genidxp,7);  xppq=gendat(genidxp,13).*baseMVA./gen(genidxp,7);
Plgen=bus(genidx,3)./baseMVA; Qlgen=bus(genidx,4)./baseMVA;
Plload=bus(loadidx,3)./baseMVA; Qlload=bus(loadidx,4)./baseMVA;

% controllers:
if size(out,2)>=2&&avridx==1; avr=out{2}; type=avr(genidxp,2); KA=avr(genidxp,3); TA=avr(genidxp,4); KE=avr(genidxp,5); TE=avr(genidxp,6); KF=avr(genidxp,7); TF=avr(genidxp,8); A=avr(genidxp,9); B=avr(genidxp,10); clear avr; else type=0; end

if pssidx
    PSS=out{3};
    PSSgen=PSS(genidxp,1); PSSinput=PSS(genidxp,2); PSSblocks=PSS(genidxp,3); Tw=PSS(genidxp,4); KPSS=PSS(genidxp,5); T1=PSS(genidxp,6); T2=PSS(genidxp,7); T3=PSS(genidxp,8); T4=PSS(genidxp,9); clear PSS; else PSSgen=[];
end

clear gendat bus gen branch area gencost nomnod linee nomare lim_area cpiu cmeno capa success Yf Yt;
%% INITIAL DATA

%   MODEL 2
if a~=0
    Epq(1:a,1)=xpd(1:a)./abs(V(genidx(1:a))).*sqrt(Pg(1:a).^2+(Qg(1:a)+abs(V(genidx(1:a))).^2./xpd(1:a)).^2); 
    delta(1:a,1)=asin(xpd(1:a).*Pg(1:a)./Epq(1:a)./abs(V(genidx(1:a))))+angle(V(genidx(1:a)));   
end
%   MODEL 4
if b~=0
    delta(a+1:a+b,1)=angle(V(genidx(a+1:a+b)))+atan(Pg(a+1:a+b)./(Qg(a+1:a+b)+abs(V(genidx(a+1:a+b))).^2./xq(a+1:a+b))); 
    Vdg(1:b,1)=abs(V(genidx(a+1:a+b))).*sin(delta(a+1:a+b)-angle(V(genidx(a+1:a+b)))); 
    Vqg(1:b,1)=abs(V(genidx(a+1:a+b))).*cos(delta(a+1:a+b)-angle(V(genidx(a+1:a+b)))); 
    Iq(1:b,1)=(Vqg(1:b)./Vdg(1:b).*Pg(a+1:a+b)-Qg(a+1:a+b))./(Vqg(1:b).^2./Vdg(1:b)+Vdg(1:b)); 
    Id(1:b,1)=(Pg(a+1:a+b)-Vqg(1:b).*Iq(1:b))./Vdg(1:b); 
    Ef(1:b,1)=Vqg(1:b)+xd(a+1:a+b).*Id(1:b); 
    Epd(1:b,1)=Vdg(1:b)-xpq(a+1:a+b).*Iq(1:b); 
    Epq(a+1:a+b,1)=Vqg(1:b)+xpd(a+1:a+b).*Id(1:b); 
end
%   MODEL 51
if c~=0
    delta(a+b+1:a+b+c,1)=angle(V(genidx(a+b+1:a+b+c)))+atan(Pg(a+b+1:a+b+c)./(Qg(a+b+1:a+b+c)+abs(V(genidx(a+b+1:a+b+c))).^2./xq(a+b+1:a+b+c))); 
    Vdg(b+1:b+c,1)=abs(V(genidx(a+b+1:a+b+c))).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); 
    Vqg(b+1:b+c,1)=abs(V(genidx(a+b+1:a+b+c))).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); 
    Iq(b+1:b+c,1)=(Vqg(b+1:b+c)./Vdg(b+1:b+c).*Pg(a+b+1:a+b+c)-Qg(a+b+1:a+b+c))./(Vqg(b+1:b+c).^2./Vdg(b+1:b+c)+Vdg(b+1:b+c)); 
    Id(b+1:b+c,1)=(Pg(a+b+1:a+b+c)-Vqg(b+1:b+c).*Iq(b+1:b+c))./Vdg(b+1:b+c); 
    Epq(a+b+1:a+b+c,1)=Vqg(b+1:b+c)+xpd(a+b+1:a+b+c).*Id(b+1:b+c);    
    Eppd(1:c,1)=Vdg(b+1:b+c)-xppq(a+b+1:a+b+c).*Iq(b+1:b+c);
    Epd(b+1:b+c,1)=Eppd(1:c)-(xpq(a+b+1:a+b+c)-xppq(a+b+1:a+b+c)+Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c))).*Iq(b+1:b+c);
    Ef(b+1:b+c,1)=Epq(a+b+1:a+b+c)+(xd(a+b+1:a+b+c)-xpd(a+b+1:a+b+c)).*Id(b+1:b+c); 
end
%   MODEL 52
if d~=0
    delta(a+b+c+1:a+b+c+d,1)=angle(V(genidx(a+b+c+1:a+b+c+d)))+atan(Pg(a+b+c+1:a+b+c+d)./(Qg(a+b+c+1:a+b+c+d)+abs(V(genidx(a+b+c+1:a+b+c+d))).^2./xq(a+b+c+1:a+b+c+d))); 
    Vdg(b+c+1:b+c+d,1)=abs(V(genidx(a+b+c+1:a+b+c+d))).*sin(delta(a+b+c+1:a+b+c+d)-angle(V(genidx(a+b+c+1:a+b+c+d)))); 
    Vqg(b+c+1:b+c+d,1)=abs(V(genidx(a+b+c+1:a+b+c+d))).*cos(delta(a+b+c+1:a+b+c+d)-angle(V(genidx(a+b+c+1:a+b+c+d)))); 
    Iq(b+c+1:b+c+d,1)=(Vqg(b+c+1:b+c+d)./Vdg(b+c+1:b+c+d).*Pg(a+b+c+1:a+b+c+d)-Qg(a+b+c+1:a+b+c+d))./(Vqg(b+c+1:b+c+d).^2./Vdg(b+c+1:b+c+d)+Vdg(b+c+1:b+c+d)); 
    Id(b+c+1:b+c+d,1)=(Pg(a+b+c+1:a+b+c+d)-Vqg(b+c+1:b+c+d).*Iq(b+c+1:b+c+d))./Vdg(b+c+1:b+c+d); 
    Eppq(1:d,1)=Vqg(b+c+1:b+c+d)+xppd(a+b+c+1:a+b+c+d).*Id(b+c+1:b+c+d);    
    Eppd(c+1:c+d,1)=Vdg(b+c+1:b+c+d)-xppq(a+b+c+1:a+b+c+d).*Iq(b+c+1:b+c+d);
    Ef(b+c+1:b+c+d,1)=Eppq(1:d)+(xd(k)-xppd(a+b+c+1:a+b+c+d)).*Id(b+c+1:b+c+d);
    Epq(a+b+c+1:a+b+c+d,1)=(1-TAA(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d)).*Ef(b+c+1:b+c+d)-(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d)-Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d))).*Id(b+c+1:b+c+d);     
end
%% K1~K10 CONSTANTS

if a~=0; K1(1:a,1)=Epq(1:a).*abs(V(genidx(1:a)))./xpd(1:a).*cos(delta(1:a)-angle(V(genidx(1:a)))); end % MODEL 2
if b~=0; K1(a+1:a+b,1)=Epd(1:b).*abs((V(genidx(a+1:a+b))))./(xpq(a+1:a+b)).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))... % MODEL 4
    +Epq(a+1:a+b).*abs((V(genidx(a+1:a+b))))./(xpd(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))...
    -abs((V(genidx(a+1:a+b)))).^2.*(1./(xpd(a+1:a+b))-1./(xpq(a+1:a+b))).*cos(2*(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))));
end
if c~=0; K1(a+b+1:a+b+c,1)=Epq(a+b+1:a+b+c).*abs(V(genidx(a+b+1:a+b+c)))./xpd(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))... % MODEL 51
    -abs(V(genidx(a+b+1:a+b+c))).^2.*(1./xpd(a+b+1:a+b+c)-1./xppq(a+b+1:a+b+c)).*cos(2*(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))))...
    +Eppd(1:c).*abs(V(genidx(a+b+1:a+b+c)))./xppq(a+b+1:a+b+c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); 
end
if d~=0; K1(a+b+c+1:a+c+b+d,1)=Eppq(1:d).*abs(V(genidx(a+b+c+1:a+c+b+d)))./xppd(a+b+c+1:a+c+b+d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))... % MODEL 52
    -abs(V(genidx(a+b+c+1:a+c+b+d))).^2.*(1./xppd(a+b+c+1:a+c+b+d)-1./xppq(a+b+c+1:a+c+b+d)).*cos(2*(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))))...
    +Eppd(c+1:c+d).*abs(V(genidx(a+b+c+1:a+c+b+d)))./xppq(a+b+c+1:a+c+b+d).*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))); 
end
   
if a~=0; K2(1:a,1)=zeros(a,1); end % MODEL 2
if b~=0; K2(a+1:a+b,1)=-abs((V(genidx(a+1:a+b))))./(xpq(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))); end % MODEL 4
if c~=0; K2(a+b+1:a+b+c,1)=-abs(V(genidx(a+b+1:a+b+c)))./xppq(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); end % MODEL 51 
if d~=0; K2(a+b+c+1:a+c+b+d,1)=-abs(V(genidx(a+b+c+1:a+c+b+d)))./xppq(a+b+c+1:a+c+b+d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))); end % MODEL 52 

if a~=0; K3(1:a,1)=zeros(a,1); end % MODEL 2
if b~=0; K3(a+1:a+b,1)=abs((V(genidx(a+1:a+b))))./(xpd(a+1:a+b)).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))); end % MODEL 4
if c~=0; K3(a+b+1:a+b+c,1)=abs(V(genidx(a+b+1:a+b+c)))./xpd(a+b+1:a+b+c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); end % MODEL 51 
if d~=0; K3(a+b+c+1:a+c+b+d,1)=abs(V(genidx(a+b+c+1:a+c+b+d)))./xppd(a+b+c+1:a+c+b+d).*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))); end % MODEL 52

K4=-K1; % MODEL 2_4_51_52

if a~=0; K5(1:a,1)=Epq(1:a)./(xpd(1:a)).*sin(delta(1:a)-angle((V(genidx(1:a))))); end % MODEL 2
if b~=0; K5(a+1:a+b,1)=Epq(a+1:a+b)./(xpd(a+1:a+b)).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))... % MODEL 4
    -Epd(1:b)./(xpq(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))-abs((V(genidx(a+1:a+b))))...
    .*(1./(xpd(a+1:a+b))-1./(xpq(a+1:a+b))).*sin(2*(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))));
end
if c~=0; K5(a+b+1:a+b+c,1)=Epq(a+b+1:a+b+c)./xpd(a+b+1:a+b+c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))... % MODEL 51
    -Eppd(1:c)./xppq(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))-abs(V(genidx(a+b+1:a+b+c)))...
    .*(-1./xppq(a+b+1:a+b+c)+1./xpd(a+b+1:a+b+c)).*sin(2*(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))); 
end
if d~=0; K5(a+b+c+1:a+c+b+d,1)=Eppq(1:d)./xppd(a+b+c+1:a+c+b+d).*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))... % MODEL 52
    -Eppd(c+1:c+d)./xppq(a+b+c+1:a+c+b+d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))...
    -abs(V(genidx(a+b+c+1:a+c+b+d))).*(-1./xppq(a+b+c+1:a+c+b+d)+1./xppd(a+b+c+1:a+c+b+d)).*sin(2*(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))); 
end

if a~=0; K6(1:a,1)=-Epq(1:a).*abs((V(genidx(1:a))))./(xpd(1:a)).*sin(delta(1:a)-angle((V(genidx(1:a))))); end % MODEL 2
if b~=0; K6(a+1:a+b,1)=abs((V(genidx(a+1:a+b)))).^2.*sin(2*(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))))... % MODEL 4
    .*(1./(xpd(a+1:a+b))-1./(xpq(a+1:a+b)))-Epq(a+1:a+b).*abs((V(genidx(a+1:a+b))))./(xpd(a+1:a+b))...
    .*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))+Epd(1:b).*abs((V(genidx(a+1:a+b))))./(xpq(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))); end
if c~=0; K6(a+b+1:a+b+c,1)=abs(V(genidx(a+b+1:a+b+c))).^2.*sin(2*(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))).*(1./xpd(a+b+1:a+b+c)-1./xppq(a+b+1:a+b+c))... % MODEL 51
    -Epq(a+b+1:a+b+c).*abs(V(genidx(a+b+1:a+b+c)))./xpd(a+b+1:a+b+c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))...
    +Eppd(1:c).*abs(V(genidx(a+b+1:a+b+c)))./xppq(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); end
if d~=0; K6(a+b+c+1:a+c+b+d,1)=abs(V(genidx(a+b+c+1:a+c+b+d))).^2.*sin(2*(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))))... % MODEL 52
    .*(1./xppd(a+b+c+1:a+c+b+d)-1./xppq(a+b+c+1:a+c+b+d))-Eppq(1:d).*abs(V(genidx(a+b+c+1:a+c+b+d)))./xppd(a+b+c+1:a+c+b+d)...
    .*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))+Eppd(c+1:c+d).*abs(V(genidx(a+b+c+1:a+c+b+d)))./xppq(a+b+c+1:a+c+b+d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))); end

if a~=0; K7(1:a,1)=zeros(a,1); end % MODEL 2
if b~=0; K7(a+1:a+b,1)=abs((V(genidx(a+1:a+b))))./(xpq(a+1:a+b)).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))); end % MODEL 4
if c~=0; K7(a+b+1:a+b+c,1)=abs(V(genidx(a+b+1:a+b+c)))./xppq(a+b+1:a+b+c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); end % MODEL 51
if d~=0; K7(a+b+c+1:a+c+b+d,1)=abs(V(genidx(a+b+c+1:a+c+b+d)))./xppq(a+b+c+1:a+c+b+d).*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))); end % MODEL 52 

if a~=0; K8(1:a,1)=zeros(a,1); end % MODEL 2
if b~=0; K8(a+1:a+b,1)=abs((V(genidx(a+1:a+b))))./(xpd(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))); end % MODEL 4
if c~=0; K8(a+b+1:a+b+c,1)=abs(V(genidx(a+b+1:a+b+c)))./xpd(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))); end % MODEL 51
if d~=0; K8(a+b+c+1:a+c+b+d,1)=abs(V(genidx(a+b+c+1:a+c+b+d)))./xppd(a+b+c+1:a+c+b+d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))); end % MODEL 52

K9=-K6; % MODEL 2_4_51_52

if a~=0; K10(1:a,1)=Epq(1:a)./xpd(1:a).*cos(delta(1:a)-angle((V(genidx(1:a)))))-2*abs((V(genidx(1:a))))./(xpd(1:a)); end % MODEL 2
if b~=0; K10(a+1:a+b,1)=Epq(a+1:a+b)./(xpd(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))... % MODEL 4
    +Epd(1:b)./(xpq(a+1:a+b)).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))-2*abs((V(genidx(a+1:a+b))))...
    ./(xpd(a+1:a+b)).*cos((delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))).^2-2*abs((V(genidx(a+1:a+b))))...
    ./(xpq(a+1:a+b)).*sin((delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))).^2; end
if c~=0; K10(a+b+1:a+b+c,1)=Epq(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))./xpd(a+b+1:a+b+c)... % MODEL 51
    +Eppd(1:c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c))))./xppq(a+b+1:a+b+c)-2*abs(V(genidx(a+b+1:a+b+c)))./xpd(a+b+1:a+b+c).*cos(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))).^2-...
    2.*abs(V(genidx(a+b+1:a+b+c)))./xppq(a+b+1:a+b+c).*sin(delta(a+b+1:a+b+c)-angle(V(genidx(a+b+1:a+b+c)))).^2; end 
if d~=0; K10(a+b+c+1:a+c+b+d,1)=Eppq(1:d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))./xppd(a+b+c+1:a+c+b+d)...
    +Eppd(c+1:c+d).*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d))))./xppq(a+b+c+1:a+c+b+d)... % MODEL 52
    -2*abs(V(genidx(a+b+c+1:a+c+b+d)))./xppd(a+b+c+1:a+c+b+d).*cos(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))).^2-...
    2*abs(V(genidx(a+b+c+1:a+c+b+d)))./xppq(a+b+c+1:a+c+b+d).*sin(delta(a+b+c+1:a+c+b+d)-angle(V(genidx(a+b+c+1:a+c+b+d)))).^2; end

%% LINEARIZED SYSTEM BLOCKS
%% A1 submatrices

A1_delta_omega=diag(ones(1,size(genidx,1))*2*pi*f);
A1_omega_delta=(diag(Sb)*diag(H))^-1*diag(K1)*(-baseMVA/2);
A1_omega_omega=-diag(D);

%   A1_omega_E

A1_omega_E_4_1=[(diag(Sb(a+1:a+b))*diag(H(a+1:a+b)))^-1*diag(K2(a+1:a+b))*(-baseMVA/2) zeros(b,c)];  
A1_omega_E_4_3=[(diag(Sb(a+1:a+b))*diag(H(a+1:a+b)))^-1*diag(K3(a+1:a+b))*(-baseMVA/2) zeros(b,c+d)];
    
A1_omega_E_51_2=[(diag(Sb(a+b+1:a+b+c))*diag(H(a+b+1:a+b+c)))^-1*diag(K2(a+b+1:a+b+c))*(-baseMVA/2) zeros(c,d)];
A1_omega_E_51_3=[zeros(c,b) (diag(Sb(a+b+1:a+b+c))*diag(H(a+b+1:a+b+c)))^-1*diag(K3(a+b+1:a+b+c))*(-baseMVA/2) zeros(c,d)]; 

A1_omega_E_52_2=[zeros(d,c) (diag(Sb(a+b+c+1:a+b+c+d))*diag(H(a+b+c+1:a+b+c+d)))^-1*diag(K2(a+b+c+1:a+b+c+d))*(-baseMVA/2)]; 
A1_omega_E_52_4=(diag(Sb(a+b+c+1:a+b+c+d))*diag(H(a+b+c+1:a+b+c+d)))^-1*diag(K3(a+b+c+1:a+b+c+d))*(-baseMVA/2);  
           
A1_omega_E=[zeros(a,2*b+3*c+3*d); A1_omega_E_4_1 zeros(b,c+d) A1_omega_E_4_3 zeros(b,d);...
   zeros(c,b+c) A1_omega_E_51_2 A1_omega_E_51_3 zeros(c,d); zeros(d,b+c) A1_omega_E_52_2 zeros(d,b+c+d) A1_omega_E_52_4];
 
%   A1_E_delta

if b~=0
A1_Epd_delta_42=diag(1./(Tpq(a+1:a+b)).*((xq(a+1:a+b))-(xpq(a+1:a+b))).*abs((V(genidx(a+1:a+b)))) ...
    ./(xpq(a+1:a+b)).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))));
A1_Epd_delta_42=[A1_Epd_delta_42; zeros(c,b)]; 
    
A1_Epq_delta_42=diag(-1./(Tpd(a+1:a+b)).*((xd(a+1:a+b))-(xpd(a+1:a+b))).*abs((V(genidx(a+1:a+b)))) ...
    ./(xpd(a+1:a+b)).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b))))));
A1_Epq_delta_42=[A1_Epq_delta_42; zeros(c+d,b)]; 
    
A1_E_delta_4=[A1_Epd_delta_42; zeros(c+d,b); A1_Epq_delta_42; zeros(d,b)];
else
    A1_E_delta_4=zeros(2*b+3*c+3*d,0);
end

if c~=0
A1_Epd_delta_513=diag(1./(Tpq(a+b+1:a+b+c)).*((xq(a+b+1:a+b+c))-(xpq(a+b+1:a+b+c)) ...
    -Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c)))...
    .*abs((V(genidx(a+b+1:a+b+c))))./(xppq(a+b+1:a+b+c)).*cos(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c))))));
A1_Epd_delta_513=[zeros(b,c); A1_Epd_delta_513];
    
A1_Eppd_delta_513=diag(1./(Tppq(a+b+1:a+b+c)).*((xpq(a+b+1:a+b+c))-(xppq(a+b+1:a+b+c)) ...
    +Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c)))...
    .*abs((V(genidx(a+b+1:a+b+c))))./(xppq(a+b+1:a+b+c)).*cos(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c))))));
A1_Eppd_delta_513=[A1_Eppd_delta_513; zeros(d,c)]; 
     
A1_Epq_delta_513=diag(-1./(Tpd(a+b+1:a+b+c)).*((xd(a+b+1:a+b+c))-(xpd(a+b+1:a+b+c))) ...
    .*abs((V(genidx(a+b+1:a+b+c))))./(xpd(a+b+1:a+b+c)).*sin(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c))))));
A1_Epq_delta_513=[zeros(b,c); A1_Epq_delta_513; zeros(d,c)]; 
 
A1_E_delta_51=[A1_Epd_delta_513; A1_Eppd_delta_513; A1_Epq_delta_513; zeros(d,c)];
else
    A1_E_delta_51=zeros(2*b+3*c+3*d,0);
end

if d~=0
    A1_Eppd_delta_524=diag(1./(Tppq(a+b+c+1:a+b+c+d)).*((xq(a+b+c+1:a+b+c+d))-(xppq(a+b+c+1:a+b+c+d))) ...
        .*abs((V(genidx(a+b+c+1:a+b+c+d))))./(xppq(a+b+c+1:a+b+c+d)).*cos(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d))))));
    A1_Eppd_delta_524=[zeros(c,d); A1_Eppd_delta_524];

    A1_Epq_delta_524=diag(-1./(Tpd(a+b+c+1:a+b+c+d)).*((xd(a+b+c+1:a+b+c+d))-(xpd(a+b+c+1:a+b+c+d)) ...
        -Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d)))...
        .*abs((V(genidx(a+b+c+1:a+b+c+d))))./(xppd(a+b+c+1:a+b+c+d)).*sin(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d))))));
    A1_Epq_delta_524=[zeros(b+c,d); A1_Epq_delta_524];
    
    A1_Eppq_delta_524=diag(-1./(Tppd(a+b+c+1:a+b+c+d)).*((xpd(a+b+c+1:a+b+c+d))-(xppd(a+b+c+1:a+b+c+d)) ...
        +Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d)))...
        .*abs((V(genidx(a+b+c+1:a+b+c+d))))./(xppd(a+b+c+1:a+b+c+d)).*sin(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d))))));

    A1_E_delta_52=[zeros(b+c,d); A1_Eppd_delta_524; A1_Epq_delta_524; A1_Eppq_delta_524];
else
        A1_E_delta_52=zeros(2*b+3*c+3*d,0);
end    

A1_E_delta=[zeros(2*b+3*c+3*d,a) A1_E_delta_4 A1_E_delta_51 A1_E_delta_52];

%   A1_E_E

if b~=0
    A1_E_E_1=[diag(-1./(Tpq(a+1:a+b)).*((xq(a+1:a+b))-(xpq(a+1:a+b)))./(xpq(a+1:a+b))-1./(Tpq(a+1:a+b))) zeros(b,b+3*c+3*d)];
    A1_E_E_5=[zeros(b,b+2*c+d) diag(-1./(Tpd(a+1:a+b)).*((xd(a+1:a+b))-(xpd(a+1:a+b)))./(xpd(a+1:a+b))-1./(Tpd(a+1:a+b))) zeros(b,c+2*d)];
else
    A1_E_E_1=zeros(0,2*b+3*c+3*d); A1_E_E_5=zeros(0,2*b+3*c+3*d);
end

if c~=0
    A1_E_E_2=[zeros(c,b) diag(-1./(Tpq(a+b+1:a+b+c)))...
        diag(-1./(Tpq(a+b+1:a+b+c)).*((xq(a+b+1:a+b+c))-(xpq(a+b+1:a+b+c))...
        -Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c)))./(xppq(a+b+1:a+b+c))) zeros(c,b+c+3*d)];
    
    A1_E_E_3=[zeros(c,b) diag(1./(Tppq(a+b+1:a+b+c)))...
        diag(-1./(Tppq(a+b+1:a+b+c)).*((xpq(a+b+1:a+b+c))-(xppq(a+b+1:a+b+c))...
        +Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c)))./(xppq(a+b+1:a+b+c))-1./(Tppq(a+b+1:a+b+c))) zeros(c,b+c+3*d)];
    
    A1_E_E_6=[zeros(c,2*b+2*c+d) diag(-1./(Tpd(a+b+1:a+b+c)).*((xd(a+b+1:a+b+c))-(xpd(a+b+1:a+b+c)))./(xpd(a+b+1:a+b+c))-1./(Tpd(a+b+1:a+b+c))) zeros(c,2*d)];
else
    A1_E_E_2=zeros(0,2*b+3*c+3*d); A1_E_E_3=zeros(0,2*b+3*c+3*d); A1_E_E_6=zeros(0,2*b+3*c+3*d); 
end

if d~=0
    A1_E_E_4=[zeros(d,b+2*c) diag(-1./(Tppq(a+b+c+1:a+b+c+d)).*((xq(a+b+c+1:a+b+c+d))-(xppq(a+b+c+1:a+b+c+d)))./(xppq(a+b+c+1:a+b+c+d))-1./(Tppq(a+b+c+1:a+b+c+d))) zeros(d,b+c+2*d)];

    A1_E_E_7=[zeros(d,2*b+3*c+d) diag(-1./(Tpd(a+b+c+1:a+b+c+d))) ...
        diag(-1./(Tpd(a+b+c+1:a+b+c+d)).*((xd(a+b+c+1:a+b+c+d))-(xpd(a+b+c+1:a+b+c+d))...
        -Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d)))./(xppd(a+b+c+1:a+b+c+d)))];
    
    A1_E_E_8=[zeros(d,2*b+3*c+d) diag(1./(Tppd(a+b+c+1:a+b+c+d))) ...
        diag(-1./(Tppd(a+b+c+1:a+b+c+d)).*((xpd(a+b+c+1:a+b+c+d))-(xppd(a+b+c+1:a+b+c+d))...
        +Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d)))./(xppd(a+b+c+1:a+b+c+d))-1./(Tppd(a+b+c+1:a+b+c+d)))];
    
else
    A1_E_E_4=zeros(0,2*b+3*c+3*d); A1_E_E_7=zeros(0,2*b+3*c+3*d); A1_E_E_8=zeros(0,2*b+3*c+3*d); 
end

A1_E_E=[A1_E_E_1; A1_E_E_2; A1_E_E_3; A1_E_E_4; A1_E_E_5; A1_E_E_6; A1_E_E_7; A1_E_E_8];

%% B1 matrix

B1_omega=[(diag(Sb)*diag(H))^-1*diag(K5)*(-baseMVA/2) (diag(Sb)*diag(H))^-1*diag(K4)*(-baseMVA/2)];

if b~=0
    B1_E_1=[zeros(b,a) diag(1./(Tpq(a+1:a+b)).*((xq(a+1:a+b))-(xpq(a+1:a+b))).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))./(xpq(a+1:a+b))) ...
        zeros(b,c+d+a) diag(-1./(Tpq(a+1:a+b)).*((xq(a+1:a+b))-(xpq(a+1:a+b))).*abs((V(genidx(a+1:a+b)))).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))./(xpq(a+1:a+b))) zeros(b,c+d)];
    B1_E_5=[zeros(b,a) diag(1./(Tpd(a+1:a+b)).*((xd(a+1:a+b))-(xpd(a+1:a+b))).*cos(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))./(xpd(a+1:a+b))) ...
        zeros(b,c+d+a) diag(1./(Tpd(a+1:a+b)).*((xd(a+1:a+b))-(xpd(a+1:a+b))).*abs((V(genidx(a+1:a+b)))).*sin(delta(a+1:a+b)-angle((V(genidx(a+1:a+b)))))./(xpd(a+1:a+b))) zeros(b,c+d)];
else
    B1_E_1=zeros(0,2*(a+b+c+d)); B1_E_5=zeros(0,2*(a+b+c+d));
end

if c~=0
    B1_E_2=[zeros(c,a+b) diag(1./(Tpq(a+b+1:a+b+c)).*((xq(a+b+1:a+b+c))-(xpq(a+b+1:a+b+c))...
        -Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c))).*sin(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c)))))./(xppq(a+b+1:a+b+c))) ...
        zeros(c,b+d+a) diag(-1./(Tpq(a+b+1:a+b+c)).*((xq(a+b+1:a+b+c))-(xpq(a+b+1:a+b+c))...
        -Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c))).*abs((V(genidx(a+b+1:a+b+c)))).*cos(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c)))))./(xppq(a+b+1:a+b+c))) zeros(c,d)];
    
    B1_E_3=[zeros(c,a+b) diag(1./(Tppq(a+b+1:a+b+c)).*((xpq(a+b+1:a+b+c))-(xppq(a+b+1:a+b+c))...
        +Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c))).*sin(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c)))))./(xppq(a+b+1:a+b+c))) ...
        zeros(c,b+d+a) diag(-1./(Tppq(a+b+1:a+b+c)).*((xpq(a+b+1:a+b+c))-(xppq(a+b+1:a+b+c))...
        +Tppq(a+b+1:a+b+c)./Tpq(a+b+1:a+b+c).*xppq(a+b+1:a+b+c)./xpq(a+b+1:a+b+c).*(xq(a+b+1:a+b+c)-xpq(a+b+1:a+b+c))).*abs((V(genidx(a+b+1:a+b+c)))).*cos(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c)))))./(xppq(a+b+1:a+b+c))) zeros(c,d)];
    
    B1_E_6=[zeros(c,a+b) diag(1./(Tpd(a+b+1:a+b+c)).*((xd(a+b+1:a+b+c))-(xpd(a+b+1:a+b+c))).*cos(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c)))))./(xpd(a+b+1:a+b+c))) ...
        zeros(c,b+d+a) diag(1./(Tpd(a+b+1:a+b+c)).*((xd(a+b+1:a+b+c))-(xpd(a+b+1:a+b+c))).*abs((V(genidx(a+b+1:a+b+c)))).*sin(delta(a+b+1:a+b+c)-angle((V(genidx(a+b+1:a+b+c)))))./(xpd(a+b+1:a+b+c))) zeros(c,d)];

else
    B1_E_2=zeros(0,2*(a+b+c+d)); B1_E_3=zeros(0,2*(a+b+c+d)); B1_E_6=zeros(0,2*(a+b+c+d));
end

if d~=0
    B1_E_4=[zeros(d,a+b+c) diag(1./(Tppq(a+b+c+1:a+b+c+d)).*((xq(a+b+c+1:a+b+c+d))-(xppq(a+b+c+1:a+b+c+d))).*sin(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d)))))./(xppq(a+b+c+1:a+b+c+d))) ...
        zeros(d,b+c+a) diag(-1./(Tppq(a+b+c+1:a+b+c+d)).*((xq(a+b+c+1:a+b+c+d))-(xppq(a+b+c+1:a+b+c+d))).*abs((V(genidx(a+b+c+1:a+b+c+d)))).*cos(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d)))))./(xppq(a+b+c+1:a+b+c+d)))];

    B1_E_7=[zeros(d,a+b+c) diag(1./(Tpd(a+b+c+1:a+b+c+d)).*((xd(a+b+c+1:a+b+c+d))-(xpd(a+b+c+1:a+b+c+d))...
        -Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d))).*cos(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d)))))./(xppd(a+b+c+1:a+b+c+d))) ...
        zeros(d,b+c+a) diag(1./(Tpd(a+b+c+1:a+b+c+d)).*((xd(a+b+c+1:a+b+c+d))-(xpd(a+b+c+1:a+b+c+d))...
        -Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d))).*abs((V(genidx(a+b+c+1:a+b+c+d)))).*sin(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d)))))./(xppd(a+b+c+1:a+b+c+d)))];
    
    B1_E_8=[zeros(d,a+b+c) diag(1./(Tppd(a+b+c+1:a+b+c+d)).*((xpd(a+b+c+1:a+b+c+d))-(xppd(a+b+c+1:a+b+c+d))...
        +Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d))).*cos(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d)))))./(xppd(a+b+c+1:a+b+c+d))) ...
        zeros(d,b+c+a) diag(1./(Tppd(a+b+c+1:a+b+c+d)).*((xpd(a+b+c+1:a+b+c+d))-(xppd(a+b+c+1:a+b+c+d))...
        +Tppd(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d).*xppd(a+b+c+1:a+b+c+d)./xpd(a+b+c+1:a+b+c+d).*(xd(a+b+c+1:a+b+c+d)-xpd(a+b+c+1:a+b+c+d))).*abs((V(genidx(a+b+c+1:a+b+c+d)))).*sin(delta(a+b+c+1:a+b+c+d)-angle((V(genidx(a+b+c+1:a+b+c+d)))))./(xppd(a+b+c+1:a+b+c+d)))];
   
else
    B1_E_4=zeros(0,2*(a+b+c+d)); B1_E_7=zeros(0,2*(a+b+c+d)); B1_E_8=zeros(0,2*(a+b+c+d));
end
    
B1_E=[B1_E_1; B1_E_2; B1_E_3; B1_E_4; B1_E_5; B1_E_6; B1_E_7; B1_E_8];

%% C1 matrix

C1_delta=[ diag(K1); diag(K6)];

if b~=0; C1_P_E4=[diag(K2(a+1:a+b)) zeros(b,2*c+d) diag(K3(a+1:a+b)) zeros(b,c+2*d)]; ...
         C1_Q_E4=[diag(K7(a+1:a+b)) zeros(b,2*c+d) diag(K8(a+1:a+b)) zeros(b,c+2*d)];...
         else C1_P_E4=zeros(0,2*b+3*c+3*d); C1_Q_E4=zeros(0,2*b+3*c+3*d); 
end
if c~=0; C1_P_E51=[zeros(c,b+c) diag(K2(a+b+1:a+b+c)) zeros(c,b+d) diag(K3(a+b+1:a+b+c)) zeros(c,2*d)];...
         C1_Q_E51=[zeros(c,b+c) diag(K7(a+b+1:a+b+c)) zeros(c,b+d) diag(K8(a+b+1:a+b+c)) zeros(c,2*d)];...
         else C1_P_E51=zeros(0,2*b+3*c+3*d); C1_Q_E51=zeros(0,2*b+3*c+3*d); 
end
if d~=0; C1_P_E52=[zeros(d,b+2*c) diag(K2(a+b+c+1:a+b+c+d)) zeros(d,b+d+c) diag(K3(a+b+c+1:a+b+c+d))];...
         C1_Q_E52=[zeros(d,b+2*c) diag(K7(a+b+c+1:a+b+c+d)) zeros(d,b+d+c) diag(K8(a+b+c+1:a+b+c+d))];...
         else C1_P_E52=zeros(0,2*b+3*c+3*d); C1_Q_E52=zeros(0,2*b+3*c+3*d); 
end

C1_E=[zeros(a,2*b+3*c+3*d); C1_P_E4; C1_P_E51; C1_P_E52; zeros(a,2*b+3*c+3*d); C1_Q_E4; C1_Q_E51; C1_Q_E52]; 

%% D1 matrix && D2 matrix
D1_P_teta=zeros(size(genidx,1),size(genidx,1)); D1_P_V=D1_P_teta; D1_Q_teta=D1_P_teta; D1_Q_V=D1_P_teta;
D2_P_teta=zeros(size(genidx,1),size(loadidx,1)); D2_P_V=D2_P_teta; D2_Q_teta=D2_P_teta; D2_Q_V=D2_P_teta;
for k=1:size(genidx,1)
    for m=1:size(genidx,1)
        if k~=m
            D1_P_teta(k,m)=-abs(V(genidx(k)))*abs(V(genidx(m)))*abs(Ybus(genidx(k),genidx(m)))*sin(angle(V(genidx(k)))-angle(V(genidx(m)))-angle(Ybus(genidx(k),genidx(m))));
            D1_P_V(k,m)=-abs(V(genidx(k)))*abs(Ybus(genidx(k),genidx(m)))*cos(angle(V(genidx(k)))-angle(V(genidx(m)))-angle(Ybus(genidx(k),genidx(m))));
            D1_Q_teta(k,m)=abs(V(genidx(k)))*abs(V(genidx(m)))*abs(Ybus(genidx(k),genidx(m)))*cos(angle(V(genidx(k)))-angle(V(genidx(m)))-angle(Ybus(genidx(k),genidx(m))));
            D1_Q_V(k,m)=-abs(V(genidx(k)))*abs(Ybus(genidx(k),genidx(m)))*sin(angle(V(genidx(k)))-angle(V(genidx(m)))-angle(Ybus(genidx(k),genidx(m))));
        else
            D1_P_teta(k,m)=abs(V(genidx(k)))^2*imag(Ybus(genidx(k),genidx(k)))-Qlgen(k)+Qg(k)+K4(k);
            D1_P_V(k,m)=(Plgen(k)-Pg(k))/abs(V(genidx(k)))-abs(V(genidx(k)))*real(Ybus(genidx(k),genidx(k)))+K5(k);
            D1_Q_teta(k,m)=abs(V(genidx(k)))^2*real(Ybus(genidx(k),genidx(k)))+Plgen(k)-Pg(k)+K9(k);
            D1_Q_V(k,m)=(Qlgen(k)-Qg(k))/abs(V(genidx(k)))+abs(V(genidx(k)))*imag(Ybus(genidx(k),genidx(k)))+K10(k);
        end
    end
    for m=1:size(loadidx,1)
        D2_P_teta(k,m)=-abs(V(genidx(k)))*abs(V(loadidx(m)))*abs(Ybus(genidx(k),loadidx(m)))*sin(angle(V(genidx(k)))-angle(V(loadidx(m)))-angle(Ybus(genidx(k),loadidx(m))));
        D2_P_V(k,m)=-abs(V(genidx(k)))*abs(Ybus(genidx(k),loadidx(m)))*cos(angle(V(genidx(k)))-angle(V(loadidx(m)))-angle(Ybus(genidx(k),loadidx(m))));
        D2_Q_teta(k,m)=abs(V(genidx(k)))*abs(V(loadidx(m)))*abs(Ybus(genidx(k),loadidx(m)))*cos(angle(V(genidx(k)))-angle(V(loadidx(m)))-angle(Ybus(genidx(k),loadidx(m))));
        D2_Q_V(k,m)=-abs(V(genidx(k)))*abs(Ybus(genidx(k),loadidx(m)))*sin(angle(V(genidx(k)))-angle(V(loadidx(m)))-angle(Ybus(genidx(k),loadidx(m))));
    end
end

D1=[D1_P_V D1_P_teta; D1_Q_V D1_Q_teta]; D2=[D2_P_V D2_P_teta; D2_Q_V D2_Q_teta];

%% D3 && D4 matrix
D3_P_teta=zeros(size(loadidx,1),size(genidx,1)); D3_P_V=D3_P_teta; D3_Q_teta=D3_P_teta; D3_Q_V=D3_P_teta;
D4_P_teta=zeros(size(loadidx,1),size(loadidx,1)); D4_P_V=D4_P_teta; D4_Q_teta=D4_P_teta; D4_Q_V=D4_P_teta;

for k=1:size(loadidx,1)
    for m=1:size(genidx,1)
        D3_P_teta(k,m)=-abs(V(loadidx(k)))*abs(V(genidx(m)))*abs(Ybus(loadidx(k),genidx(m)))*sin(angle(V(loadidx(k)))-angle(V(genidx(m)))-angle(Ybus(loadidx(k),genidx(m))));
        D3_P_V(k,m)=-abs(V(loadidx(k)))*abs(Ybus(loadidx(k),genidx(m)))*cos(angle(V(loadidx(k)))-angle(V(genidx(m)))-angle(Ybus(loadidx(k),genidx(m))));
        D3_Q_teta(k,m)=abs(V(loadidx(k)))*abs(V(genidx(m)))*abs(Ybus(loadidx(k),genidx(m)))*cos(angle(V(loadidx(k)))-angle(V(genidx(m)))-angle(Ybus(loadidx(k),genidx(m))));
        D3_Q_V(k,m)=-abs(V(loadidx(k)))*abs(Ybus(loadidx(k),genidx(m)))*sin(angle(V(loadidx(k)))-angle(V(genidx(m)))-angle(Ybus(loadidx(k),genidx(m))));
    end
    for m=1:size(loadidx,1)
        if k~=m
            D4_P_teta(k,m)=-abs(V(loadidx(k)))*abs(V(loadidx(m)))*abs(Ybus(loadidx(k),loadidx(m)))*sin(angle(V(loadidx(k)))-angle(V(loadidx(m)))-angle(Ybus(loadidx(k),loadidx(m))));
            D4_P_V(k,m)=-abs(V(loadidx(k)))*abs(Ybus(loadidx(k),loadidx(m)))*cos(angle(V(loadidx(k)))-angle(V(loadidx(m)))-angle(Ybus(loadidx(k),loadidx(m))));
            D4_Q_teta(k,m)=abs(V(loadidx(k)))*abs(V(loadidx(m)))*abs(Ybus(loadidx(k),loadidx(m)))*cos(angle(V(loadidx(k)))-angle(V(loadidx(m)))-angle(Ybus(loadidx(k),loadidx(m))));
            D4_Q_V(k,m)=-abs(V(loadidx(k)))*abs(Ybus(loadidx(k),loadidx(m)))*sin(angle(V(loadidx(k)))-angle(V(loadidx(m)))-angle(Ybus(loadidx(k),loadidx(m))));
        else
            D4_P_teta(k,m)=abs(V(loadidx(k)))^2*imag(Ybus(loadidx(k),loadidx(k)))-Qlload(k);
            D4_P_V(k,m)=Plload(k)/abs(V(loadidx(k)))-abs(V(loadidx(k)))*real(Ybus(loadidx(k),loadidx(k)));
            D4_Q_teta(k,m)=abs(V(loadidx(k)))^2*real(Ybus(loadidx(k),loadidx(k)))+Plload(k);
            D4_Q_V(k,m)=Qlload(k)/abs(V(loadidx(k)))+abs(V(loadidx(k)))*imag(Ybus(loadidx(k),loadidx(k)));
        end
    end
end
D3=[D3_P_V D3_P_teta; D3_Q_V D3_Q_teta]; D4=[D4_P_V D4_P_teta; D4_Q_V D4_Q_teta];

%% Adding controls
%% AVR

if avridx
if b~=0||c~=0; A1_E_xv_451=[diag(1./Tpd(a+1:a+b+c)) zeros(b+c,d+2*trace(diag(type==2)))]; else A1_E_xv_451=zeros(0,b+c+d+2*trace(diag(type==2))); end
if d~=0; A1_E_xv_52=[zeros(d,a+b+c) diag(1./Tpd(a+b+c+1:a+b+c+d).*(1-TAA(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d))) zeros(d,2*trace(diag(type==2)));...
        [zeros(d,a+b+c) diag(TAA(a+b+c+1:a+b+c+d)./Tpd(a+b+c+1:a+b+c+d)./Tppd(a+b+c+1:a+b+c+d)) zeros(d,2*trace(diag(type==2)))]]; else A1_E_xv_52=zeros(0,b+c+d+2*trace(diag(type==2)));
end
A1_E_xv=[zeros(b+2*c+d,b+c+d+2*trace(diag(type==2))); A1_E_xv_451; A1_E_xv_52];

if trace(diag(type==2))==0
    A1_xv_xv=diag(-1./TA); B1_xv_V=[zeros(b+c+d,a) diag(-KA(a+1:a+b+c+d)./TA(a+1:a+b+c+d))]; B1_xv=[B1_xv_V zeros(size(B1_xv_V))];
else
    elim=[];
    for k=1:size(genidx,1)
        if type(k)~=2; elim=[elim k]; end
    end
    
    A1_xv_xv_Ef_Ef=diag(type==1)*diag(-1./TA)+diag(type==2)*diag((-KE-A.*exp(B.*(Ef))-(Ef).*A.*B.*exp(B.*(Ef)))./TE);
    A1_xv_xv_Ef_VR=diag(1./TE); A1_xv_xv_Ef_VR(:,elim)=[]; 
    A1_xv_xv_Ef_RF=zeros(b+c+d,trace(diag(type==2)));
    
    A1_xv_xv_VR_Ef=diag(-KA.*KF./TA./TF)*diag(type==2); A1_xv_xv_VR_Ef(elim,:)=[];
    A1_xv_xv_VR_VR=diag(-1./TA); A1_xv_xv_VR_VR(:,elim)=[]; A1_xv_xv_VR_VR(elim,:)=[]; 
    A1_xv_xv_VR_RF=diag(KA./TA); A1_xv_xv_VR_RF(:,elim)=[]; A1_xv_xv_VR_RF(elim,:)=[];

    A1_xv_xv_RF_Ef=diag(KF./TF.^2)*diag(type==2); A1_xv_xv_RF_Ef(elim,:)=[];
    A1_xv_xv_RF_VR=zeros(trace(diag(type==2)),trace(diag(type==2))); 
    A1_xv_xv_RF_RF=diag(-1./TF); A1_xv_xv_RF_RF(:,elim)=[]; A1_xv_xv_RF_RF(elim,:)=[];
         
    A1_xv_xv=[[A1_xv_xv_Ef_Ef A1_xv_xv_Ef_VR A1_xv_xv_Ef_RF]; [A1_xv_xv_VR_Ef A1_xv_xv_VR_VR A1_xv_xv_VR_RF]; [A1_xv_xv_RF_Ef A1_xv_xv_RF_VR A1_xv_xv_RF_RF]];

    B1_xv_V_Ef=[zeros(b+c+d,a) diag(-KA(a+1:a+b+c+d)./TA(a+1:a+b+c+d))*diag(type==1)];
    B1_xv_V_VR=diag(-KA(a+1:a+b+c+d)./TA(a+1:a+b+c+d)); B1_xv_V_VR(elim,:)=[]; B1_xv_V_VR=[zeros(size(B1_xv_V_VR,1),a) B1_xv_V_VR];
    B1_xv_V_RF=zeros(size(B1_xv_V_VR,1),a+b+c+d);
    
    B1_xv_V=[B1_xv_V_Ef; B1_xv_V_VR; B1_xv_V_RF]; B1_xv=[B1_xv_V zeros(size(B1_xv_V))];
end

% control matrix E
    E_xv_Ef=diag(type==1)*diag(KA./TA);
    if trace(diag(type==2))~=0; E_xv_VR=diag(type==2)*diag(KA./TA); E_xv_VR(elim,:)=[]; E_xv_RF=zeros(size(E_xv_VR)); else E_xv_VR=zeros(0,b+c+d); E_xv_RF=zeros(0,b+c+d); end
    E_xv=[E_xv_Ef; E_xv_VR; E_xv_RF];
    E=[zeros( 2*a+4*b+5*c+5*d,size(E_xv,2)); E_xv];

%% Adding PSS aditional blocks
    if pssidx
        kpss=1; elim1=[]; for k=1:size(genidx,1); if ~(PSSgen(kpss)==genidx(k)); elim1=[elim1 k]; elseif kpss<size(PSSgen,1); kpss=kpss+1; end; end
        PSSindicator=diag(ones(1,d+b+c));
        PSSindicator(elim1,:)=[]; 
        
        A1_xv_xPSS1=zeros(b+c+d,size(PSSgen,1)); %%%%
        A1_xv_xPSS2=zeros(b+c+d,size(PSSgen,1));
        A1_xv_xPSS3=PSSindicator.'*(diag(KA(PSSgen)./TA(PSSgen))*PSSindicator)*diag(type==1); A1_xv_xPSS3(:,elim1)=[];
        A1_xv_xPSSI=[A1_xv_xPSS1 A1_xv_xPSS2 A1_xv_xPSS3];
        
        A1_xPSS_delta=[A1_omega_delta(PSSgen,PSSgen)*diag(KPSS)*PSSindicator; A1_omega_delta(PSSgen,PSSgen)*diag(KPSS.*T1./T2)*PSSindicator;...
            A1_omega_delta(PSSgen,PSSgen)*diag(KPSS.*T1.*T3./T2./T4)*PSSindicator]; %%%%%
         
        A1_xPSS_Omega=[A1_omega_omega(PSSgen,PSSgen)*diag(KPSS)*PSSindicator; A1_omega_omega(PSSgen,PSSgen)*diag(KPSS.*T1./T2)*PSSindicator;...
            A1_omega_omega(PSSgen,PSSgen)*diag(KPSS.*T1.*T3./T2./T4)*PSSindicator]; %%%%%
        
        for k=1:size(PSSgen,1); A1_xPSS_E(k,:)=A1_omega_E(PSSgen(k),:)*KPSS(k); A1_xPSS_E(k+size(PSSgen,1),:)=A1_omega_E(PSSgen(k),:)*KPSS(k)*T1(k)/T2(k);
        A1_xPSS_E(k+2*size(PSSgen,1),:)=A1_omega_E(PSSgen(k),:)*KPSS(k)*T1(k)/T2(k)*T3(k)/T4(k); end %%%%     
               
        if trace(diag(type==2))~=0
            A1_xv_xPSS4=zeros(trace(diag(type==2)),size(PSSgen,1)); %%%%
            A1_xv_xPSS5=zeros(trace(diag(type==2)),size(PSSgen,1));
            A1_xv_xPSS6=PSSindicator.'*(diag(KA(PSSgen)./TA(PSSgen))*PSSindicator)*diag(type==2); A1_xv_xPSS6(elim,:)=[]; A1_xv_xPSS6(:,elim1)=[];
            A1_xv_xPSSII=[[A1_xv_xPSS4 A1_xv_xPSS5 A1_xv_xPSS6]; zeros(size(A1_xv_xPSS6,1),size(A1_xv_xPSSI,2))]; %%%%
        else
            A1_xv_xPSSII=zeros(0,size(A1_xv_xPSSI,2));
        end
        A1_xv_xPSS=[A1_xv_xPSSI; A1_xv_xPSSII];      
        A1_xPSS_xPSS=[[diag(-1./Tw)  zeros(size(PSSgen,1),2*size(PSSgen,1))]; [diag((Tw-T1)./Tw./T2) diag(-1./T2) zeros(size(PSSgen,1),size(PSSgen,1))];...
            [diag((Tw-T1).*T3./Tw./T2./T4) diag((T2-T3)./T2./T4) diag(-1./T4)]];
        
        B1_xPSS_omega=[[B1_omega(PSSgen,PSSgen)*diag(KPSS)*PSSindicator B1_omega(PSSgen,PSSgen+size(genidx,1))*diag(KPSS)*PSSindicator];...
            [B1_omega(PSSgen,PSSgen)*diag(KPSS.*T1./T2)*PSSindicator B1_omega(PSSgen,PSSgen+size(genidx,1))*diag(KPSS.*T1./T2)*PSSindicator];...
            [B1_omega(PSSgen,PSSgen)*diag(KPSS.*T1.*T3./T2./T4)*PSSindicator B1_omega(PSSgen,PSSgen+size(genidx,1))*diag(KPSS.*T1.*T3./T2./T4)*PSSindicator]];
    else
         A1_xv_xPSS=zeros(size(A1_xv_xv,1),0); A1_xPSS_xPSS=[]; A1_xPSS_Omega=zeros(0,a+b+c+d); A1_xPSS_delta=zeros(0,a+b+c+d); A1_xPSS_E=zeros(0,2*b+3*c+3*d);
         B1_xPSS_omega=zeros(0,2*(a+b+c+d));
    end   
else 
    A1_E_xv=zeros(size(A1_E_E,1),0);
    A1_xv_xv=[]; A1_xv_xPSS=[];
    A1_xPSS_xPSS=[]; A1_xPSS_Omega=zeros(0,a+b+c+d); A1_xPSS_delta=zeros(0,a+b+c+d); A1_xPSS_E=zeros(0,2*b+3*c+3*d);
    B1_xv=zeros(0,2*(a+b+c+d)); B1_xPSS_omega=zeros(0,2*(a+b+c+d)); E=[];
end

%% Composing matrices
A1=[                   zeros(a+b+c+d,a+b+c+d)                                A1_delta_omega                       zeros(a+b+c+d,2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+pssidx*(size(PSSgen,1)*3))
                           A1_omega_delta                                    A1_omega_omega                                   A1_omega_E        zeros(a+b+c+d,avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+pssidx*(size(PSSgen,1)*3))
                             A1_E_delta                                 zeros(2*b+3*c+3*d,a+b+c+d)                           A1_E_E A1_E_xv     zeros(size(A1_E_E,1),pssidx*(size(PSSgen,1)*3))
    zeros(avridx*(trace(diag(type==1))+3*trace(diag(type==2))),a+b+c+d)                     zeros(avridx*(trace(diag(type==1))+3*trace(diag(type==2))),size(A1_E_E,1)+size(A1_omega_omega,1))         A1_xv_xv A1_xv_xPSS
                            A1_xPSS_delta                                      A1_xPSS_Omega                                 A1_xPSS_E zeros(pssidx*(size(PSSgen,1)*3),avridx*(trace(diag(type==1))+3*trace(diag(type==2)))) A1_xPSS_xPSS];

B1=[        zeros(a+b+c+d,2*(a+b+c+d))
                     B1_omega
                      B1_E 
                     B1_xv
                B1_xPSS_omega];            

C1=[C1_delta zeros(2*(a+b+c+d),a+b+c+d) C1_E zeros(2*(a+b+c+d),avridx*(b+c+d+2*trace(diag(type==2)))) zeros(2*(a+b+c+d),pssidx*(size(PSSgen,1)*3))];

clear K1 K2 K3 K4 K5 K6 K7 K8 K9 K10 A B;
A=A1-B1*inv(D1-D2*inv(D4)*D3)*C1;
B=E; 

%% Creating state variable list
pos_var=cell(2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+pssidx*(size(PSSgen,1)*3),2);

for k=1:size(genidx,1); pos_var{k,1}=k; pos_var{k,2}=strcat('delta_G',num2str(genidx(k))); pos_var{size(genidx,1)+k,1}=size(genidx,1)+k; pos_var{size(genidx,1)+k,2}=strcat('Omega_G',num2str(genidx(k))); end
if b~=0; for k=1:b; pos_var{2*size(genidx,1)+k,1}=2*size(genidx,1)+k; pos_var{2*size(genidx,1)+k,2}=strcat('Epd_G',num2str(genidx(k)));
        pos_var{2*size(genidx,1)+b+2*c+d+k,1}=2*size(genidx,1)+2*c+d+b+k; pos_var{2*size(genidx,1)+2*c+b+d+k,2}=strcat('Epq_G',num2str(genidx(k))); end; end
if c~=0; for k=1:c; pos_var{2*size(genidx,1)+b+k,1}=2*size(genidx,1)+b+k; pos_var{2*size(genidx,1)+b+k,2}=strcat('Epd_G',num2str(genidx(k+a+b)));
        pos_var{2*size(genidx,1)+b+c+k,1}=2*size(genidx,1)+c+b+k; pos_var{2*size(genidx,1)+c+b+k,2}=strcat('Eppd_G',num2str(genidx(k+a+b)));
        pos_var{2*size(genidx,1)+2*b+2*c+k,1}=2*size(genidx,1)+2*c+2*b+k; pos_var{2*size(genidx,1)+2*c+2*b+k,2}=strcat('Epq_G',num2str(genidx(k+a+b)));end; end
if d~=0; for k=1:d; pos_var{2*size(genidx,1)+b+2*c+k,1}=2*size(genidx,1)+b+2*c+k; pos_var{2*size(genidx,1)+b+2*c+k,2}=strcat('Eppd_G',num2str(genidx(k+a+b+c)));
        pos_var{2*size(genidx,1)+2*b+2*c+d+k,1}=2*size(genidx,1)+2*c+d+2*b+k; pos_var{2*size(genidx,1)+2*c+d+2*b+k,2}=strcat('Epq_G',num2str(genidx(k+a+b+c)));
        pos_var{2*size(genidx,1)+2*b+3*c+2*d+k,1}=2*size(genidx,1)+3*c+2*d+2*b+k; pos_var{2*size(genidx,1)+3*c+2*d+2*b+k,2}=strcat('Eppq_G',num2str(genidx(k+a+b+c)));end; end
if avridx; for k=1:size(genidx,1); pos_var{2*size(genidx,1)+2*b+3*c+3*d+k,1}=2*size(genidx,1)+2*b+3*c+3*d+k; pos_var{2*size(genidx,1)+2*b+3*c+3*d+k,2}=strcat('Ef_G',num2str(genidx(k))); 
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+k+3,1}=2*size(genidx,1)+2*b+3*c+3*d+k+3; pos_var{2*size(genidx,1)+2*b+3*c+3*d+k+3,2}=strcat('VR_G',num2str(genidx(k))); 
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+k+6,1}=2*size(genidx,1)+2*b+3*c+3*d+k+6; pos_var{2*size(genidx,1)+2*b+3*c+3*d+k+6,2}=strcat('Rf_G',num2str(genidx(k)));end; end
if pssidx; for k=1:size(PSSgen,1); pos_var{2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+k,1}=2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+k;
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+k,2}=strcat('V1_G',num2str(PSSgen(k)));
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+size(PSSgen,1)+k,1}=2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+size(PSSgen,1)+k;
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+size(PSSgen,1)+k,2}=strcat('V2_G',num2str(PSSgen(k)));
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+2*size(PSSgen,1)+k,1}=2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+2*size(PSSgen,1)+k;
        pos_var{2*size(genidx,1)+2*b+3*c+3*d+avridx*(trace(diag(type==1))+3*trace(diag(type==2)))+2*size(PSSgen,1)+k,2}=strcat('V3_G',num2str(PSSgen(k)));
    end
end

state_var=diag((model==2)*2)+diag((model==4)*4)+diag((model==51)*5)+diag((model==52)*5);
if avridx; state_var=state_var++avridx*diag([zeros(a,1); (type==1)*1])+avridx*diag([zeros(a,1); (type==2)*3]); end
if pssidx; for k=1:size(PSSgen,1); state_var=state_var+diag((genidx==PSSgen(k))*3); end; end
state_var=[genidx state_var*ones(size(genidx,1),1)];
end