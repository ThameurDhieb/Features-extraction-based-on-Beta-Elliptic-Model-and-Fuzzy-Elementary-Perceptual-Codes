function [points,V,DPV,T,t,X_init,Y_init,T_init, pas_t] = copie_pre_traitement_seg_H_p(X,Y,T, rayon_filtre_V, sigma_p_filtre_V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction pré_trait retourne les signaux de position X Y filtrés                                     %%%%%
%%%%% en utilisant la méthode FTRPB                                                                          %%%%%
%%%%% ainsi que les signaux de vitesse d'accélération ,le signal de vitesse curviligne                       %%%%%
%%%%% et  la vitesse angulaire filtrés en utilisant la méthode IIRPB                                         %%%%%
%%%%% avec  B et A sont les coefficients des filtres                                                         %%%%%
%%%%% les paramètres d'entrée sont X, Y et Z sont les données originales (tablette)                          %%%%%
%%%%% les paramètres de sortie points:(X,Y) et vitesse curviligne V et la dirivée seconde du signal vitesse  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

points = [];
Vfst = [];
longt = length(T);

%%% pas_t = 0.01;
pas_t = T(2) - T(1);

%%%%%%%%%%%%%%%%%%%% 
temp_depart = T(1);
for i = 1 : longt
    T(i) = T(i) - temp_depart;
end
%%%%%%%%%%%%%%%%%%%%
j = 1;
k = 0;
i = 1;
M = length(X);
while i < M + 1    
      k = k + 1;
      points(k,j) = X(i);
      points(k,j+1) = Y(i);       
      i = i + 1;
end
t = T;

%%%%%%% long : taille de la matrice initaile
long = length(t);

%%%%%%% il faut que la longueur min de la matrice "points"=46 pour ne pas avoir le
%%%%%%% problème du filtre; l'idée est de remplire le reste d'une autre matrice
%%%%%%% "points46" par les coordonnées du dérnier point de la matrice "points"
%%%%%%% traitement pour X
points46 = [];
t46 = [];
if long <= 20
   for i = 1 : long - 1          
       points46 = [points46;points(i,:);(points(i,:)+points(i+1,:))/2];
       t46 = [t46;t(i);(t(i)+t(i+1))/2];
   end
   i = long;
   points46 = [points46;points(i,:)];
   t46 = [t46;t(i)];
   points = points46;
   t = t46;
   T = t;
   pp = size(points);
   
   %%% pas_t = 0.005;
   pas_t = T(2) - T(1);
end


% T_init = T;
% X_init = points(:,1);
% Y_init = points(:,2);


%------------------détermination des paramètres du filtre utilisé  
%------------------A et B les coefficients du filtre

%[B,A] = iirpb1(5,100,16,20);
if size(points,1) > (3 * 8)
   [B,A] = iirpb1(8,80,16,20); % résérvé au mots obtenus par le système de re-échantillonnage
elseif size(points,1) > (3 * 5)
   [B,A] = iirpb1(5,80,16,20);
else
   lar_fil = fix( size(points,1) / 3 );
   [B,A] = iirpb1(lar_fil,100,16,20);
end

%------------------filtrer les signaux de position X et Y par la méthode time reversal 

%points(:,1) = filtfilt(B,A,points(:,1));
%points(:,2) = filtfilt(B,A,points(:,2));
rayon = 3;%4;
sigma_p = 2;
if ((rayon * 2) + 1) > size(points,1)
%   rayon = max( round( size(points,1) / 2 ) , 1);
   rayon = size(points,1);
end

[points,XYpoint] = filtre_lineaire_1(rayon, sigma_p, 1, size(points,1), points);

T_init = T;
X_init = points(:,1);
Y_init = points(:,2);
 
%------------------filtrer les signaux de vitesse VX et VY et d'acceleration AX et AY
%------------------1-détermination des paramètres du filtre utilisé  A et B pour filtrer les signaux de vitesse VX et VY

VX = vitessex(points(:,1),t);
VY = vitessey(points(:,2),t);

%------------------filtrer les signaux de vitesse VX et VY et d'acceleration AX et AY


N = 6; Fs = 100; Fc = 12; Flag = 0; beta = 6; ripple = 35;

if length(VX) < ( N * 3)
   N = fix( length(VX) / 3);
end

%Wind : la fenetre de poids (windows) qui peut etre:
%            lanczos(N), hamming(N), hanning(N), blackman(N), bartlett(N)
%            triang(N), kaiser(N,beta), chebwim(N,ripple)
Wind = hamming(N);
%Wind = kaiser(N,beta);
%Wind = chebwin(N,ripple)
[B] = firpb(N,Fs,Fc,Wind,Flag);
A = 1;
VX = filtfilt(B,A,VX);
%%%%%%%%VY = filtfilt(B,A,VY);
V = sqrt(VX.^2+VY.^2);

%rayon_filtre_V = 2;
%sigma_p_filtre_V = 0.7;

[V] = filtre_lineaire(rayon_filtre_V, sigma_p_filtre_V, V);

%----------------------- La dérivée première de la vitesse V
dt = diff(t);
dV = diff(V);
DPV = (dV)./(dt);
DPV = [0; DPV];


%[B,A] = iirpb1(10,100,10,20);
%DPV = filtfilt(B,A,DPV);

%HBHB% N = 6; Fs = 100; Fc = 16; Flag = 0; beta = 6; ripple = 34;
%HBHB% N = 3; Fs = 100; Fc = 16; Flag = 0; beta = 6; ripple = 34;
 N = 2; Fs = 100; Fc = 16; Flag = 0; beta = 6; ripple = 34;  % HBH : meilleur compromis entre detection 
                                                             %       des double inflexions de la vitesse
                                                             %       (pas de filtrage de DPV) 
                                                             %       et diminution des segments excedentaires
                                                             %       (filtrage accentiué de DPV au debut at à la fin du tracé) 

if length(VX) < ( N * 3)
   N = fix( length(VX) / 3);
end
                                                             
%Wind : la fenetre de poids (windows) qui peut etre:
%            lanczos(N), hamming(N), hanning(N), blackman(N), bartlett(N)
%            triang(N), kaiser(N,beta), chebwim(N,ripple)
%Wind = hamming(N);
%Wind = kaiser(N,beta)
Wind = chebwin(N,ripple);
[B] = firpb(N,Fs,Fc,Wind,Flag);
A = 1;

DPV = filtfilt(B,A,DPV);%HBHBHBHB%

AX = accelerationx(VX,t);
AY = accelerationy(VY,t);

%-----------------2-filtrer les signaux d'accéleration AX et AY
AXf = filtfilt(B,A,AX);
AYf = filtfilt(B,A,AY);
AXY = sqrt(AX.^2+AY.^2);
AXYf = filtfilt(B,A,AXY);
 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

