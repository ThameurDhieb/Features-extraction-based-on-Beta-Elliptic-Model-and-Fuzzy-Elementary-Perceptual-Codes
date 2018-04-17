function [points,V,DPV,T,t,X_init,Y_init,T_init] = copie_pre_traitement_seg_H(X,Y,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction pr�_trait retourne les signaux de position X Y filtr�s                                     %%%%%
%%%%% en utilisant la m�thode FTRPB                                                                          %%%%%
%%%%% ainsi que les signaux de vitesse d'acc�l�ration ,le signal de vitesse curviligne                       %%%%%
%%%%% et  la vitesse angulaire filtr�s en utilisant la m�thode IIRPB                                         %%%%%
%%%%% avec  B et A sont les coefficients des filtres                                                         %%%%%
%%%%% les param�tres d'entr�e sont X, Y et Z sont les donn�es originales (tablette)                          %%%%%
%%%%% les param�tres de sortie points:(X,Y) et vitesse curviligne V et la diriv�e seconde du signal vitesse  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

points = [];
Vfst = [];
longt = length(T);

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
%%%%%%% probl�me du filtre; l'id�e est de remplire le reste d'une autre matrice
%%%%%%% "points46" par les coordonn�es du d�rnier point de la matrice "points"
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
end


% T_init = T;
% X_init = points(:,1);
% Y_init = points(:,2);


%------------------d�termination des param�tres du filtre utilis�  
%------------------A et B les coefficients du filtre

[B,A] = iirpb1(5,100,16,20);

%------------------filtrer les signaux de position X et Y par la m�thode time reversal 

points(:,1) = filtfilt(B,A,points(:,1));
points(:,2) = filtfilt(B,A,points(:,2));

T_init = T;
X_init = points(:,1);
Y_init = points(:,2);
 
%------------------filtrer les signaux de vitesse VX et VY et d'acceleration AX et AY
%------------------1-d�termination des param�tres du filtre utilis�  A et B pour filtrer les signaux de vitesse VX et VY

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

%----------------------- La d�riv�e premi�re de la vitesse V
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
                                                             %       (filtrage accentiu� de DPV au debut at � la fin du trac�) 
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

%-----------------2-filtrer les signaux d'acc�leration AX et AY
AXf = filtfilt(B,A,AX);
AYf = filtfilt(B,A,AY);
AXY = sqrt(AX.^2+AY.^2);
AXYf = filtfilt(B,A,AXY);
 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

