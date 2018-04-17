function [Co]=firpb(N,Fs,Fc,Wind,Flag);
% FIRPB est un programme qui retourne les parametres de la reponse
% en frequence du filtre de N point de type passe-bas a phase lineaire
% ainsi que les coefficients du filtre.
% LES PARAMETRES D'ENTREE SONT :
%     N    : nombre de point du filtre
%     Fs   : la frequence d'echantillonnage en Hz
%     Fc   : la frequence de coupure en Hz
%     Wind : la fenetre de poids (windows) qui peut etre :
%            lanczos(N), hamming(N), hanning(N), blackman(N), bartlett(N)
%            triang(N), kaiser(N,beta), chebwim(N,ripple)
%     Flag : si Flag=1, il y aura un fichier de sortie BUFFER.MAT
%            valeur par defaut = 0
%
% LES PARAMETRES DE SORTIE SONT :
%        reponse en frequence : H    (complexe)
%        vecteur de frequence : W    (Hz) de 0 a Fs/2
%        erreur relative      : Er   (%)
%        coefficient          : Co
%
% Pour les graphiques, voir la fonction : GRAFPB

flag=0;
if nargin < 4
    win = boxcar(N+1);
elseif nargin > 4
    win = Wind;
    flag = Flag;
else
    win = Wind;
end;

Wn=Fc/(Fs/2);

%------------- COEFFICIENT DU FILTRE -------------

Co = fir1(N-1,Wn,win);

return;
