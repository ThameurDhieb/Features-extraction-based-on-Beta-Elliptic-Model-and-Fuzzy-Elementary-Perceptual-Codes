function [points,V,DPV,T,t,X_init,Y_init,T_init, pas_t] = copie_pre_traitement_seg_H_Re_p(X, Y, T, rayon_filtre_V, sigma_p_filtre_V)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction pré_trait retourne les signaux de position X Y filtrés                                   %%%%%
%%%%% en utilisant la méthode FTRPB                                                                        %%%%%
%%%%% ainsi que les signaux de vitesse d'accélération ,le signal de vitesse curviligne                     %%%%%
%%%%% et  la vitesse angulaire filtrés en utilisant la méthode IIRPB                                       %%%%%
%%%%% avec  B et A sont les coefficients des filtres                                                       %%%%%
%%%%% les paramètres d'entrée sont X, Y et Z sont les données originales (tablette)                        %%%%%
%%%%% les paramètres de sortie points:(X,Y) et vitesse curviligne V et la dirivée seconde du signal vitesse%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
hh = size(points(:,:));
t = T;

%%%% long : taille de la matrice initaile
long = length(t);

%%%%%%%%%%%%%%%% il faut que la longueur min de la matrice "points"=46 pour ne pas avoir le
%%%%%%%%%%%%%%%% problème du filtre; l'idée est de remplire le reste d'une autre matrice
%%%%%%%%%%%%%%%% "points46" par les coordonnées de la dérnière point de la matrice "points"
%%%%%%%%%%%%%%%% traitement pour X
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


T_init = T;
X_init = points(:,1);
Y_init = points(:,2);



 
%------------------filtrer les signaux de vitesse VX et VY

VX = vitessex(points(:,1),t);
VY = vitessey(points(:,2),t);

%------------------détermination des paramètres du filtre utilisé A et B les coefficients du filtre
N = 6; Fs = 100; Fc = 12; Flag = 0; beta = 6; ripple = 35;
if length(VX) < ( N * 3)
   N = fix( length(VX) / 3);
end

%Wind : la fenetre de poids (windows) qui peut etre:
%            lanczos(N), hamming(N), hanning(N), blackman(N), bartlett(N)
%            triang(N), kaiser(N,beta), chebwim(N,ripple)
Wind = hamming(N);
%Wind = kaiser(N,beta);
%Wind = chebwin(N,ripple);
[B] = firpb(N,Fs,Fc,Wind,Flag);
A = 1;
VX = filtfilt(B,A,VX);
V = sqrt(VX.^2+VY.^2);

%rayon_filtre_V = 2;
%sigma_p_filtre_V = 0.7;
[V] = filtre_lineaire(rayon_filtre_V, sigma_p_filtre_V, V);

%----------------------- La dérivée première de la vitesse V
dt = diff(t);
dV = diff(V);
DPV = (dV)./(dt);
DPV = [0; DPV];

N = 2; Fs = 100; Fc = 16; Flag = 0; beta = 6; ripple = 34;

if length(VX) < ( N * 3)
   N = fix( length(VX) / 3);
end

%Wind : la fenetre de poids (windows) qui peut etre:
%            lanczos(N), hamming(N), hanning(N), blackman(N), bartlett(N)
%            triang(N), kaiser(N,beta), chebwim(N,ripple);
%Wind = hamming(N);
%Wind = kaiser(N,beta);
Wind = chebwin(N,ripple);
[B] = firpb(N,Fs,Fc,Wind,Flag);
A = 1;

DPV = filtfilt(B,A,DPV);%HBHBHBHB%

AX = accelerationx(VX,t);
AY = accelerationy(VY,t);

%-----------------2-filtrer les signaux d'acc%eleration AX et AY

AXf = filtfilt(B,A,AX);
AYf = filtfilt(B,A,AY);
AXY = sqrt(AX.^2+AY.^2);
AXYf = filtfilt(B,A,AXY);



return;