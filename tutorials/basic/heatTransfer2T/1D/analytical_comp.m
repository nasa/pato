clear all
clc

fiche = fopen('postProcessing/singleGraph/15/line_T_g_T_s.xy','r') ;

% Lecture des données qui se trouvent sur 3 colonnes
% Quelque soit la quantité de données

a = fscanf(fiche,'%f ',[3 inf]);
a = a.';
xf = a(:,1);
Tf1 = a(:,2);

% Fermeture du fichier texte
fclose(fiche);

T0 = 293;
Ts = 323;
u_x = 1e-2;

hsf = 6000;
eps = 0.8;
rho_g = 8.248e-01;
cp_g = 1676;

a = hsf / (eps*rho_g*cp_g*u_x);

T = (T0 - Ts) * exp(-a*xf) + Ts;

figure
plot(xf,T,xf,Tf1,'o')

save('AnalyticalT.txt','T')