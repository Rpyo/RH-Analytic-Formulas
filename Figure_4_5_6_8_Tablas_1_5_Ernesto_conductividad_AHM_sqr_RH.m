clearvars;close all;
clc 
clear all
% Conductivity gain predictions for multiscale fibrous composites
% with interfacial thermal barrier resistance
% Ernesto Iglesias-Rodríguez*1
% | Julián Bravo-Castillero1
% | Manuel E. Cruz2
% | Raul Guinovart-Díaz



%% Figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bi=[1000, 100];
xas=500;
cm=1;
% p=[0:.001:.1, 0.1:.01:1];
% alpha=(1-p)./(p*(.1-1)+1);
alpha=[0:.001:.1, 0.1:.01:1];
phi1=[.2:.1:.6, .65,.68:.01:.74];
 for pp=1:length(Bi)
  
  K=Bi(pp); % imperfect parameter
  cf=xas;
  for ii=1:length(alpha)
      for  jj=1:length(phi1)
      V=phi1(jj);
      lamCH = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,V,K);
%       lamCH1(jj)=lamCH;
      Vc=alpha(ii)*V;
      Vf=V*(1-alpha(ii))/(1-Vc);
      lamInt = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,Vf,K);
      lamRH = AHMelastico_IMPERFECT_sqr_orden3a(lamInt,cf,Vc,K);
%       lamRH1(ii,jj)=lamRH;
      lamGain1(ii,jj,pp)=lamRH/lamCH;
  
  
      end
    end
 end

colores=[[0,0,1];[.6,.6,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];

    figure
 for pp=1:length(Bi)
subplot(2,2,2*(pp-1)+1) % 1 y 3
hold on
for jj=1:length(phi1(1:5))
plot(alpha,lamGain1(:,jj,pp),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
 
end
plot(alpha,ones(size(alpha)),':k','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
for jj=1:length(phi1(1:5))                   
text(.4,lamGain1(136,jj,pp),['\textbf{$\phi=$}',num2str(phi1(jj),3)],'fontsize',12,'interpreter','latex')
end
hold off
ylim([1, 1.1])
xlim([0, 1])
box on
grid on
grid minor
ylabel('$\kappa_{gain}$','fontsize',20,'interpreter','latex')
if pp==2
xlabel('aggregation parameter ($\alpha$)','fontsize',16,'interpreter','latex')
end
title( ['$Bi=$',num2str(Bi(pp)),',    ','$\rho=500$'],'fontsize',12,'interpreter','latex')
 end

  for pp=1:length(Bi)
subplot(2,2,2*(pp-1)+2) % 2 y 4
hold on
for jj=5:length(phi1)
plot(alpha,lamGain1(:,jj,pp),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
 
end
plot(alpha,ones(size(alpha)),':k','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
for jj=5:length(phi1)                   
text(.4,lamGain1(136,jj,pp),['\textbf{$\phi=$}',num2str(phi1(jj),3)],'fontsize',12,'interpreter','latex')
end
hold off
ylim([0.82, 1.1])
yticks([0.82:.04: 1.1])
xlim([0, 1])
box on
grid on
grid minor
% ylabel('$\kappa_{gain}$','fontsize',20,'interpreter','latex')
if pp==2
xlabel('aggregation parameter ($\alpha$)','fontsize',16,'interpreter','latex')
end
title( ['$Bi=$',num2str(Bi(pp)),',    ','$\rho=500$'],'fontsize',12,'interpreter','latex')
 end

%% END figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bi=[1e10];
xas=[ 2,6 , 50, 500];
% cf=50; %%% propiedad de conductividad de las fibras
cm=1;
% p=[0:.001:.1, 0.1:.01:1];
% alpha=(1-p)./(p*(.1-1)+1);
alpha=[0:.001:.1, 0.1:.01:1];
 phi2=[.5,.6, 0.65,.67,.68 , .7];

 for mm=1:length(xas)
  phi=phi2;
  K=Bi; % parámetro de imperfeccion
  cf=xas(mm);
  for ii=1:length(alpha)
      for   jj=1:length(phi)
      V=phi(jj);
      lamCH = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,V,K);
%       lamCH1(jj)=lamCH;
      Vc=alpha(ii)*V;
      Vf=V*(1-alpha(ii))/(1-Vc);
      lamInt = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,Vf,K);
      lamRH = AHMelastico_IMPERFECT_sqr_orden3a(lamInt,cf,Vc,K);
%       lamRH1(ii,jj)=lamRH;
      lamGain2(ii,jj,mm)=lamRH/lamCH;
  
  
      end
    end
 end


%%%% Figura 5 del paper Julian Ernesto
marca1=[':b';':y';':g';':r';':c';':m';':b';':y';':m'];
marca=['-b';'-y';'-g';'-r';'-c';'-m';'-b';'-y';'-m'];

% colores=[[0,0,1];[1,1,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];

    figure
 for mm=1:length(xas)
subplot(2,2,mm)
hold on
for jj=1:length(phi)
plot(alpha,lamGain2(:,jj,mm),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
 
end
plot(alpha,ones(size(alpha)),':k','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
for jj=1:length(phi)                   
text(.4,lamGain2(136,jj,mm),['\textbf{$\phi=$}',num2str(phi(jj),3)],'fontsize',12,'interpreter','latex')
end
hold off
% ylim([0.98, 1.08])
xlim([0, 1])
box on
grid on
grid minor
if mm==1 | mm==3
ylabel('$\kappa_{gain}$','fontsize',20,'interpreter','latex')
end
if mm==3 | mm==4
xlabel('aggregation parameter ($\alpha$)','fontsize',16,'interpreter','latex')
end
title( ['$\rho=$',num2str(xas(mm),3),',    ','$Bi=10^{10}$'],'fontsize',12,'interpreter','latex')
 end

%% END figure 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Figure 6
Bi=[.1,.01, .00001];
xas=500;
% cf=50; %%% propiedad de conductividad de las fibras
cm=1;
% p=[0:.001:.1, 0.1:.01:1];
% alpha=(1-p)./(p*(.1-1)+1);
alpha=[0:.001:.1, 0.1:.01:1];
phi1=[.2:.1:.6, .65,.68:.01:.74];
 for pp=1:length(Bi)
  
  K=Bi(pp); % parámetro de imperfeccion
  cf=xas;
  for ii=1:length(alpha)
      for  jj=1:length(phi1)
      V=phi1(jj);
      lamCH = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,V,K);
%       lamCH1(jj)=lamCH;
      Vc=alpha(ii)*V;
      Vf=V*(1-alpha(ii))/(1-Vc);
      lamInt = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,Vf,K);
      lamRH = AHMelastico_IMPERFECT_sqr_orden5a(lamInt,cf,Vc,K);
%       lamRH1(ii,jj)=lamRH;
      lamGain3(ii,jj,pp)=lamRH/lamCH;
  
  
      end
    end
 end

% colores=[[0,0,1];[1,1,0];[0,1,0];[1,0,0];[0,1,1];[1,0,1];[0,0,0];[0.5,0.5,0.5];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];

    figure % 6
 for pp=1:length(Bi)
subplot(3,2,2*(pp-1)+1) % 1,  3 y 5
hold on
for jj=1:length(phi1(1:5))
plot(alpha,lamGain3(:,jj,pp),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
 
end
plot(alpha,ones(size(alpha)),':k','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
for jj=1:length(phi1(1:5))                   
text(.6,lamGain3(152,jj,pp),['\textbf{$\phi=$}',num2str(phi1(jj),3)],'fontsize',12,'interpreter','latex')
end
hold off

ylim([0.91, 1.])
yticks([0.91:.03: 1.])
xlim([0, 1])
if pp==1
ylim([0.96, 1.])
yticks([0.96:.02: 1.])
end
box on
grid on
grid minor
ylabel('$\kappa_{gain}$','fontsize',20,'interpreter','latex')
if pp==3
xlabel('aggregation parameter ($\alpha$)','fontsize',16,'interpreter','latex')
end
title( ['$Bi=$',num2str(Bi(pp),5),',    ','$\rho=500$'],'fontsize',12,'interpreter','latex')
 end

  for pp=1:length(Bi)
subplot(3,2,2*(pp-1)+2) % 2 , 4 6
hold on
for jj=5:length(phi1)
plot(alpha,lamGain3(:,jj,pp),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
 
end
plot(alpha,ones(size(alpha)),':k','LineWidth',1,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
for jj=5:length(phi1)                   
text(.6,lamGain3(152,jj,pp),['\textbf{$\phi=$}',num2str(phi1(jj),3)],'fontsize',12,'interpreter','latex')
end
hold off

ylim([0.916, 1.26])
yticks([0.92:.08: 1.26])
xlim([0, 1])
if pp==1
ylim([0.96, 1.15])
yticks([0.96:.03: 1.15])
end
box on
grid on
grid minor
% ylabel('$\kappa_{gain}$','fontsize',20,'interpreter','latex')
if pp==3
xlabel('aggregation parameter ($\alpha$)','fontsize',16,'interpreter','latex')
end
title( ['$Bi=$',num2str(Bi(pp),6),',    ','$\rho=500$'],'fontsize',12,'interpreter','latex')
 end

%% END figure 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 
 
 %% FIGURE 8  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
Bi=[.0005,5, 50000];
xas=[2,50];
% cf=50; %%% propiedad de conductividad de las fibras
cm=1;
% p=[0:.001:.1, 0.1:.01:1];
% alpha=(1-p)./(p*(.1-1)+1);
load Sk_Tk_30_mas_5_675_hasta_90_ENoP
No=15;
ang=14;
H1=(conj(delta1(ang))*conj(w2(ang))-conj(delta2(ang))*conj(w1(ang)))/(w1(ang)*conj(w2(ang))-w2(ang)*conj(w1(ang)));
H2=((delta1(ang))*conj(w2(ang))-(delta2(ang))*conj(w1(ang)))/(w1(ang)*conj(w2(ang))-w2(ang)*conj(w1(ang)));
S1=Sk(:,ang);

phi=[.73:.001:0.785,.7854];
for kk=1:length( Bi)
   K=Bi(kk); % parámetro de imperfeccion
    for mm=1:length(xas)
       cf=xas(mm);
          for jj=1:length(phi)
           V=phi(jj);
           lamCH1(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden1a(cm,cf,V,K);
           lamCH2(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden2a(cm,cf,V,K);
           lamCH3(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,V,K);
           lamCH4(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden4a(cm,cf,V,K);
           lamCH4b(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden4b(cm,cf,V,K);
           lamCH5(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,V,K);
           lamCH6(jj,mm,kk) = AHMelastico_IMPERFECT_sqr_orden6a(cm,cf,V,K);

          end
       end
end

lamCH7=AHMelastico_IMPERFECT_Realpozo_2011; % 
%  colores=[[0,0,1];[.7,.7,0];[1,0,0];[0,1,0];[0,.1,1];[1,0,1];[0,0,0];[0.9,0.5,0.5];[0.6,0.6,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
 colores=[[0,0,1];[1,0,0];[0,1,1];[.8,0.8,0];[0,0,1];[1,0,1];[1.,0,0];[0.,0.,0];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
%***************************************************************
Rocha_Cruz=[0.07903,0.07903;1.3048,2.8948;1.6766,9.5331];% "Rocha & Cruz 2001" Tabla 4, Realpozo 2011

    figure % 8
    jc=0;
 for pp=1:length(Bi)
     for mm=1:length(xas)
subplot(3,2,2*(pp-1)+mm) % 1,  3 y 5
hold on
plot(phi,lamCH1(:,mm,pp),'-','Color',colores(1,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
plot(phi,lamCH2(:,mm,pp),'-','Color',colores(2,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
plot(phi,lamCH3(:,mm,pp),'-','Color',colores(3,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
plot(phi,lamCH4(:,mm,pp),'-','Color',colores(4,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
plot(phi,lamCH4b(:,mm,pp),':','Color',colores(5,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
plot(phi,lamCH5(:,mm,pp),'-','Color',colores(6,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
plot(phi,lamCH6(:,mm,pp),':','Color',colores(7,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4) 
                   jc=jc+1;
plot(phi,lamCH7(:,jc),':','Color',colores(8,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4) 
plot(.75,Rocha_Cruz(pp,mm),'o','Color',colores(9,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)   
 
% for jj=1:length(phi1(1:5))                   
% text(.6,lamGain3(152,jj,pp),['\textbf{$\phi=$}',num2str(phi1(jj),3)],'fontsize',12,'interpreter','latex')
% end
hold off

% ylim([0.91, 1.])
% yticks([0.91:.03: 1.])
 xlim([0.73, 0.7854])
xticks([0.73:.01:.78, .785])
box on
grid on
grid minor
ylabel('$\widehat{\kappa}$','fontsize',16,'interpreter','latex')
if pp==3
xlabel('Fiber volume fraction ($\phi$)','fontsize',12,'interpreter','latex')
end
title( ['$\rho=$',num2str(xas(mm),5),', ','$Bi=$',num2str(Bi(pp),5)],'fontsize',12,'interpreter','latex')
 end
 end
 leg={'1st order','2nd order','3rd order','4th order','4th order(*)','5th order','6th order(*)','Realpozo et al. 2011','Rocha and Cruz, 2001'};
 legend(leg{:},'location','NorthWest','fontsize',12,'NumColumns',9,'interpreter','latex');

%% END figure 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Table 1
phi=[0.1:0.1:0.6];
cm=1;
cf=1/3*10^(-6);
K=1e10;
lam1=zeros(length(phi),7);
for jj=1:length(phi)
    V=phi(jj);
           lam1(jj,1) = AHMelastico_IMPERFECT_sqr_orden1a(cm,cf,V,K);
           lam1(jj,2) = AHMelastico_IMPERFECT_sqr_orden2a(cm,cf,V,K);
           lam1(jj,3) = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,V,K);
           lam1(jj,4) = AHMelastico_IMPERFECT_sqr_orden4a(cm,cf,V,K);
           lam1(jj,5) = AHMelastico_IMPERFECT_sqr_orden4b(cm,cf,V,K);
           lam1(jj,6) = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,V,K);
           lam1(jj,7) = AHMelastico_IMPERFECT_sqr_orden6a(cm,cf,V,K);
end 
 disp('Tabla 1 paper Ernesto, solo las cuatro primeras columnas')
 format long
%% %Tabla
disp('Tabla 1 de A1, con AHM')
phi=phi';
F_1=lam1(:,1);
F_2=lam1(:,2);
F_3=lam1(:,3);
F_4=lam1(:,4);

T = table(phi,F_1,F_2,F_3,F_4);
disp(T)


%% Table 2, Tabla3

phi=[0.3,0.75];
cm=1;
ro=[2,50];

lam1=zeros(length(Bi),7,length(ro),length(phi));
for kk=1:length(phi)
for ii=1:length(ro)
 Bi=[1e-7,1e-5,1e-3,1e-2,1e-1,1,ro(ii)/(ro(ii)-1),10,100,1000,1e5,1e7];
for jj=1:length(Bi)
    K=Bi(jj);cf=ro(ii); V=phi(kk);
           lam1(jj,1,ii,kk) = AHMelastico_IMPERFECT_sqr_orden1a(cm,cf,V,K);
           lam1(jj,2,ii,kk) = AHMelastico_IMPERFECT_sqr_orden2a(cm,cf,V,K);
           lam1(jj,3,ii,kk) = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,V,K);
           lam1(jj,4,ii,kk) = AHMelastico_IMPERFECT_sqr_orden4a(cm,cf,V,K);
           lam1(jj,5,ii,kk) = AHMelastico_IMPERFECT_sqr_orden4b(cm,cf,V,K);
           lam1(jj,6,ii,kk) = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,V,K);
           lam1(jj,7,ii,kk) = AHMelastico_IMPERFECT_sqr_orden6a(cm,cf,V,K);
end 
end
end
 disp('Tabla 2 paper Ernesto, solo las cuatro primeras columnas')
 format long
%% %Tabla 2 y 3, ro=2, phi=.3
for kk=1:length(phi)
   for ii=1:length(ro)
disp(['\rho=',num2str(ro(ii)),'   ','\phi=',num2str(phi(kk))] )
disp('*******************************************************************************************************')  
Bi=[1e-7,1e-5,1e-3,1e-2,1e-1,1,ro(ii)/(ro(ii)-1),10,100,1000,1e5,1e7]';
F_1=lam1(:,1,ii,kk);
F_2=lam1(:,2,ii,kk);
F_3=lam1(:,3,ii,kk);
F_4=lam1(:,4,ii,kk);

T = table(Bi,F_1,F_2,F_3,F_4);
disp(T)
disp('*******************************************************************************************************') 
   end
end

%% Tabla 5  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modelo analitico del paper
lamf=19; %%% propiedad de conductividad de las fibras
lamm=1;
cm=lamm;cf=lamf;%% propiedad de conductividad de la fibra
% ro=(lamf-1)/(lamf+1);
ro=(cf-cm)/(cf+cm);
xas=cf/cm;
% ro=1/5;
r = [0.10, .13, .16, .19]; % radio de las fibras
R=[5:30]; % parámetro de imperfeccion
Bi=[xas./R(1:end)];

% Bi=[0.143, 0.133,0.125,0.118,0.111,0.105,0.1,0.095,0.091]
n=4;
alpha=.4;
for ii=1:length(r)
   for jj=1:length(Bi)
      K=Bi(jj); % imperfect parameter
      V(ii)=n*pi*(r(ii))^2;
      lamCH(jj,ii) = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,V(ii),K);
%       lamCH1(jj)=lamCH;
      Vc=alpha*V(ii);
      Vf=V(ii)*(1-alpha)/(1-Vc);
      lamInt(jj,ii) = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,Vf,K);
      lamRH(jj,ii) = AHMelastico_IMPERFECT_sqr_orden5a(lamInt(jj,ii),cf,Vc,K);
end
end


%%%% Figura 3 del paper
 colores=[[0,0,1];[1,0,0];[0,1,1];[.8,0.8,0];[0,0,1];[1,0,1];[1.,0,0];[0.,0.,0];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
figure
hold on
for jj=1:length(r)
plot(R,lamRH(:,jj),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
end
hold off
% ylim([0.6, 1.4])
 xlim([5, 30])
box on
grid on
grid minor
% textoA=['R= ',num2str(R(1)),'R= ',num2str(R(2)), 'R= ',num2str(R(3)), 'R= ',num2str(R(4)),'R= ',num2str(R(5)),...
%     'R= ',num2str(R(6)),  'R= ',num2str(R(7))];
lgd=legend(['r= ',num2str(r(1))],['r= ',num2str(r(2))], ['r= ',num2str(r(3))], ['r= ',num2str(r(4))], 'location','northwest');
title(lgd,'AHM-RH')
xlabel('R=\rho/Bi')
ylabel('$\widehat{\kappa}$','fontsize',16,'interpreter','latex')

%% %Tabla  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Tabla 1 de A1, con AHM')

disp('Tabla 1 de A1, con AHM')
R_=(R');
r_1=[0.1;lamRH(:,1)];
r_2=[0.13;(lamRH(:,2))];
r_3=[0.16; (lamRH(:,3))];
r_4=[0.19;(lamRH(:,4))];
lamAHM=[r_1,r_2,r_3,r_4];

% T = table(categorical({'R';'14';'15';'16';'17';'18';'19';'20';'21';'22'}),lamAHM);
% disp(T)
 
% Rc=2*ro/(1-ro)
%%%%%%%%%%%%%%%%%%%%%


%%% Modelo analitico del paper
lamf=19; %%% propiedad de conductividad de las fibras
lamm=1;
cm=lamm;
cf=lamf;%% propiedad de conductividad de la fibra
% ro=(lamf-1)/(lamf+1);
ro=(cf-cm)/(cf+cm);
xas=cf/cm;
% ro=1/5;
r = [0.10, .13, .16, .19]; % radio de las fibras
 R=[5:30]; % parámetro de imperfeccion
% Bi=[1e12, xas./R(2:end)];
R1=[25:-1:5];
% Bi=[0:.1:(xas/(xas-1)-.1),xas/(xas-1),(xas/(xas-1)+.1):.1:(2*xas-1)/(xas-1)];
 Bi=[0:.1:.7, xas./R1];

n=4;
alpha=.6;
for ii=1:length(r)
   for jj=1:length(Bi)
      K=Bi(jj); % imperfect parameter
      V(ii)=n*pi*(r(ii))^2;
      lamCH(jj,ii) = AHMelastico_IMPERFECT_sqr_orden2a(cm,cf,V(ii),K);
%       lamCH1(jj)=lamCH;
      Vc=alpha*V(ii);
      Vf=V(ii)*(1-alpha)/(1-Vc);
      lamInt(jj,ii) = AHMelastico_IMPERFECT_sqr_orden2a(cm,cf,Vf,K);
      lamRH1(jj,ii) = AHMelastico_IMPERFECT_sqr_orden2a(lamInt(jj,ii),cf,Vc,K);
end
end


%%%% Figura 3 del paper
 colores=[[0,0,1];[1,0,0];[0,1,1];[.8,0.8,0];[0,0,1];[1,0,1];[1.,0,0];[0.,0.,0];[0.5,0.5,1];[0.8,0.8,0];[0.1,0.8,0];[.8,0,0];[0,.8,.8]];
figure
hold on
for jj=1:length(r)
plot(Bi,lamRH1(:,jj),'-','Color',colores(jj,:),'LineWidth',2,...
                       'MarkerEdgeColor','k',...
                       'MarkerFaceColor','w',...
                       'MarkerSize',4)
end
hold off
% ylim([0.6, 1.4])
xlim([min(Bi), max(Bi)])
box on
grid on
grid minor
% textoA=['R= ',num2str(R(1)),'R= ',num2str(R(2)), 'R= ',num2str(R(3)), 'R= ',num2str(R(4)),'R= ',num2str(R(5)),...
%     'R= ',num2str(R(6)),  'R= ',num2str(R(7))];
lgd=legend(['\phi= ',num2str(V(1))],['\phi= ',num2str(V(2))], ['\phi= ',num2str(V(3))], ['\phi= ',num2str(V(4))], 'location','northwest');
title(lgd,'AHM-RH')
xlabel('Bi')
ylabel('$\widehat{\kappa}$','fontsize',16,'interpreter','latex')




%% %Tabla
disp('Tabla 1 de A1, con AHM')
 Bi=Bi';
V_1=lamRH1(:,1);
V_2=lamRH1(:,2);
V_3=lamRH1(:,3);
V_4=lamRH1(:,4);
V_i={'\phi= ',num2str(V(1)),'\phi= ',num2str(V(2)), '\phi= ',num2str(V(3)),'\phi= ',num2str(V(4))};
T = table(Bi,V_1,V_2,V_3,V_4);
disp(V_i)
disp(T)

% Functions used

function RR = AHMelastico_IMPERFECT_sqr_orden1a(cm,cf,phi,K)
Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
RR=cm*(1-phi*Beta1)/(1+phi*Beta1);
end


function RR = AHMelastico_IMPERFECT_sqr_orden2a(cm,cf,phi,K)

Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
Beta3=((1-cf/cm)*K + 3*cf)/((1+cf/cm)*K+3*cf);
C4=0.305827833;
Tau =C4*Beta1*Beta3*phi^4;
RR=cm*(1-phi*Beta1-Tau)/(1+phi*Beta1-Tau);
end


function RR = AHMelastico_IMPERFECT_sqr_orden3a(cm,cf,phi,K)

Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
Beta3=((1-cf/cm)*K + 3*cf)/((1+cf/cm)*K+3*cf);
Beta5=((1-cf/cm)*K + 5*cf)/((1+cf/cm)*K+5*cf);
C4=0.305827833;  b8=1.40295995;
Tau =C4*Beta1*Beta3*phi^4/(1-b8*Beta3*Beta5*phi^8);
RR=cm*(1-phi*Beta1-Tau)/(1+phi*Beta1-Tau);
end

function RR = AHMelastico_IMPERFECT_sqr_orden4a(cm,cf,phi,K)

Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
Beta3=((1-cf/cm)*K + 3*cf)/((1+cf/cm)*K+3*cf);
Beta5=((1-cf/cm)*K + 5*cf)/((1+cf/cm)*K+5*cf);
Beta7=((1-cf/cm)*K + 7*cf)/((1+cf/cm)*K+7*cf);
C4=0.305827833;  b8=1.402959953646316; b12=2.55915216; 
Delta1=1-b8*Beta3*Beta5*phi^8; Delta2=Delta1-b12*Beta5*Beta7*phi^12;
 C8=0.013361523; %Delta2=(1-5041/169*b8*Beta3*Beta5*phi^8);
Tau =C4*Beta1*Beta3*phi^4/Delta2 + C8*Beta1*Beta7*phi^8*Delta1/Delta2;
RR=cm*(1-phi*Beta1-Tau)/(1+phi*Beta1-Tau);
end


function RR = AHMelastico_IMPERFECT_sqr_orden5a(cm,cf,phi,K)

Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
Beta3=((1-cf/cm)*K + 3*cf)/((1+cf/cm)*K+3*cf);
Beta5=((1-cf/cm)*K + 5*cf)/((1+cf/cm)*K+5*cf);
Beta7=((1-cf/cm)*K + 7*cf)/((1+cf/cm)*K+7*cf);
Beta9=((1-cf/cm)*K + 9*cf)/((1+cf/cm)*K+9*cf);
C4=0.305827833;  b8=1.402959953646316; b12=2.55915216; 
Delta1=1-b8*Beta3*Beta5*phi^8;b16=5.76867559;
Delta0=1-b12*Beta5*Beta7*phi^12-b16*Beta7*Beta9*phi^16 ; Delta2=Delta1-b12*Beta5*Beta7*phi^12;
C8=0.013361523;C24= 4.93057350; b12p=0.55233049;
Delta1p=Delta1-b12*Beta5*Beta7*phi^12;
Delta2p=Delta1p+168/13*b8*Beta3*Beta5*phi^8+1001/17*b12p*Beta3*Beta9*phi^12;
Delta3=Delta2-b12p*Beta3*Beta9*phi^12-b16*Beta7*Beta9*phi^16-C24*Beta3*Beta5*Beta7*Beta9*phi^24;
Tau =C4*Beta1*Beta3*phi^4*Delta0/Delta3 + C8*Beta1*Beta7*phi^8*Delta2p/Delta3;
RR=cm*(1-phi*Beta1-Tau)/(1+phi*Beta1-Tau);
end


function RR =AHMelastico_IMPERFECT_sqr_orden4b(cm,cf,phi,K)
Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
Beta3=((1-cf/cm)*K + 3*cf)/((1+cf/cm)*K+3*cf);
Beta5=((1-cf/cm)*K + 5*cf)/((1+cf/cm)*K+5*cf);
Beta7=((1-cf/cm)*K + 7*cf)/((1+cf/cm)*K+7*cf);
C4=0.305827833;
C8=0.013361523;
b8=1.402959953646316;  
D1=1-b8*Beta3*Beta5*phi^8;
Tau4MT=C4*Beta1*Beta3*phi^4/D1 +...
       C8*Beta1*Beta7*phi^8;
  RR=cm*(1-phi*Beta1-Tau4MT)/(1+phi*Beta1-Tau4MT);
end


function RR =AHMelastico_IMPERFECT_sqr_orden6a(cm,cf,phi,K)
%%%% orden 6to tabla
Beta1=((1-cf/cm)*K+cf)/((1+cf/cm)*K+cf);
Beta3=((1-cf/cm)*K + 3*cf)/((1+cf/cm)*K+3*cf);
Beta5=((1-cf/cm)*K + 5*cf)/((1+cf/cm)*K+5*cf);
Beta7=((1-cf/cm)*K + 7*cf)/((1+cf/cm)*K+7*cf);
Beta9=((1-cf/cm)*K + 9*cf)/((1+cf/cm)*K+9*cf);
C4=0.305827833;
C8=0.013361523;
C12=0.000184643;C16 = 0.242252;
C20 = 0.0341942;C24= 4.93057350;
C28 = 0.0479731;
b8=1.402959953646316; b12=2.55915216; 
bb12=0.15233; b20 = 3.59039;
b24 = 0.389837; bb24 = 6.54926;
b32 = 9.18835; b36 = 0.997652;
D0= 1 - b12*Beta5*Beta7*phi^12;
D2= 1 - b8*Beta3*Beta5*phi^8/D0 - bb12*Beta3*Beta9*phi^12;
D3= 1 - b8*Beta3*Beta5*phi^8/D0 - (b12*Beta5*Beta7+ bb12*Beta3*Beta9)*phi^12 +... 
   b20*Beta3*Beta5^2*Beta7*phi^20/D0 + b24*Beta3*Beta5*Beta7*Beta9*phi^24;
D4=1 - b8*Beta3*Beta5*phi^8/D0-(2*b12*Beta5*Beta7 + bb12*Beta3*Beta9)*phi^12+... 
   2*b20*Beta3*Beta5^2*Beta7*phi^20/D0+(2*b24*Beta3*Beta5*Beta7*Beta9+... 
    bb24*Beta5*Beta7^2)*phi^24-b32*Beta3*Beta5^3*Beta7^2*phi^32/D0-... 
   b36*Beta3*Beta5^2*Beta7^2*Beta9^2*phi^36;
Tau6MT=C4*Beta1*Beta3*phi^4/D2 +...
       C8*Beta1*Beta7*phi^8 +... 
       C12*Beta1^2*phi^12 +...
       C16*Beta1*Beta3*Beta5*Beta7*phi^16/D3 +...
       C20*Beta1*Beta5*Beta7^2*phi^20/D0+...
       C28*Beta1*Beta3*Beta5^2*Beta7^2*phi^28/D4;
  RR=cm*(1-phi*Beta1-Tau6MT)/(1+phi*Beta1-Tau6MT);
end

%%% Imperfecto sqr Realpozo 2011
function lamCH7 =AHMelastico_IMPERFECT_Realpozo_2011

%     Bi=0.0005		Bi=5		Bi=50000	
% ro=2	ro=50	ro=2	ro=50	ro=2	ro=50
	
lamCH7=[0.1031	0.1031	1.2955	2.7938	1.6517	7.8525
	   0.1019	0.1019	1.2959	2.7987	1.6529	7.919
	   0.1008	0.1008	1.2964	2.8036	1.6541	7.9869
	   0.0996	0.0996	1.2969	2.8085	1.6554	8.0563
	   0.0985	0.0985	1.2973	2.8135	1.6566	8.1272
	   0.0973	0.0973	1.2978	2.8185	1.6579	8.1998
	   0.0962	0.0962	1.2983	2.8234	1.6591	8.2739
	   0.095	0.095	1.2987	2.8284	1.6603	8.3498
	   0.0938	0.0938	1.2992	2.8335	1.6616	8.4275
	   0.0926	0.0926	1.2997	2.8385	1.6628	8.5071
	   0.0914	0.0914	1.3001	2.8435	1.6641	8.5886
	   0.0902	0.0902	1.3006	2.8486	1.6653	8.6721
	   0.089	0.089	1.301	2.8537	1.6666	8.7577
	   0.0878	0.0878	1.3015	2.8587	1.6678	8.8456
	   0.0866	0.0866	1.302	2.8638	1.6691	8.9358
	   0.0854	0.0854	1.3024	2.869	1.6703	9.0284
	   0.0841	0.0841	1.3029	2.8741	1.6716	9.1236
	   0.0829	0.0829	1.3034	2.8792	1.6729	9.2214
	   0.0816	0.0816	1.3038	2.8844	1.6741	9.322
	   0.0803	0.0803	1.3043	2.8896	1.6754	9.4256
	   0.079	0.079	1.3048	2.8948	1.6766	9.5323
	   0.0777	0.0777	1.3052	2.9	    1.6779	9.6422
	   0.0764	0.0764	1.3057	2.9052	1.6792	9.7556
	   0.0751	0.0751	1.3062	2.9105	1.6805	9.8726
	   0.0737	0.0737	1.3066	2.9157	1.6817	9.9935
	   0.0724	0.0724	1.3071	2.921	1.683	10.1184
	   0.071	0.071	1.3076	2.9263	1.6843	10.2477
	   0.0696	0.0696	1.308	2.9316	1.6855	10.3816
	   0.0682	0.0682	1.3085	2.9369	1.6868	10.5204
	   0.0668	0.0668	1.309	2.9423	1.6881	10.6644
	   0.0654	0.0654	1.3094	2.9476	1.6894	10.8141
	   0.0639	0.0639	1.3099	2.953	1.6907	10.9697
	   0.0624	0.0624	1.3104	2.9584	1.6919	11.1317
	   0.0609	0.0609	1.3109	2.9638	1.6932	11.3007
	   0.0594	0.0594	1.3113	2.9692	1.6945	11.4771
	   0.0578	0.0578	1.3118	2.9747	1.6958	11.6617
	   0.0563	0.0563	1.3123	2.9801	1.6971	11.855
	   0.0546	0.0546	1.3127	2.9856	1.6984	12.0578
	   0.053	0.053	1.3132	2.9911	1.6997	12.2711
	   0.0513	0.0513	1.3137	2.9966	1.701	12.4959
	   0.0496	0.0496	1.3141	3.0021	1.7023	12.7333
	   0.0478	0.0478	1.3146	3.0077	1.7036	12.9848
	   0.046	0.046	1.3151	3.0132	1.7049	13.2518
	   0.0441	0.0441	1.3155	3.0188	1.7062	13.5362
	   0.0422	0.0422	1.316	3.0244	1.7075	13.8403
	   0.0402	0.0402	1.3165	3.03	1.7088	14.1668
	   0.0381	0.0381	1.317	3.0356	1.7101	14.5188
	   0.0359	0.0359	1.3174	3.0413	1.7114	14.9004
	   0.0336	0.0336	1.3179	3.047	1.7127	15.3166
	   0.0312	0.0312	1.3184	3.0527	1.7141	15.7737
	   0.0286	0.0286	1.3188	3.0584	1.7154	16.2801
	   0.0259	0.0259	1.3193	3.0641	1.7167	16.8468
	   0.0229	0.0229	1.3198	3.0698	1.718	17.4889
	   0.0196	0.0196	1.3203	3.0756	1.7193	18.2281
	   0.0158	0.0158	1.3207	3.0814	1.7207	19.096
	   0.0115	0.0115	1.3212	3.0872	1.722	20.1424
	   0.0095	0.0095	1.3214	3.0895	1.7225	20.6278];

end