function [B,A]=iirpb1(N,Fs,Fc,Rs);
%function [ H,W,Er,B,A]=iirpb(N,Fs,Fc,Sorte,Flag);
% IIRPB est un programme qui retourne les parametres de la reponse
% en frequence du filtre de type passe-bas de N points en
% utilisant un prototype analogique.
%
% [H,W,Er,B,A]=iirpb(N,Fs,Fc,Sorte,Flag)
%
% LES PARAMETRES D'ENTREE SONT :
%     N     : nombre de points du filtre
%     Fs    : la frequence d'echantillonnage en Hz
%     Fc    : la frequence de coupure en Hz
%     Sorte : type de filtre     = chebychev2
%                           
%             valeur par defaut = 1
%     Flag  : si Flag=1, il y aura un fichier de sortie BUFFER.MAT
%             valeur par defaut = 0
%
% LES PARAMETRES DE SORTIES SONT :
%        reponse en frequence : H    (complexe)
%        vecteur de frequence : W    (Hz) de 0 a Fs/2
%        erreur relative      : Er   (%)
%        coefficient          : B    numerateur
%                             : A    denominateur
%
% Pour les graphique, voir la fonction : GRAFPB
%Fc=16;
%Fs=100;
%N=6;
%Rs=30;

flag=0;              % valeur par defaut
Wn=Fc/(Fs)   ;
%  _за))азз__ии-('"й&!!mщ% entre 0 et 1 ou 1 = Fs/2
[B,A]=cheby2(N-1,Rs,Wn);
%[B,A]=hamming(N-1);

%if nargin < 4
 %   sort=1;         % valeur par defaut
 %elseif nargin > 4
 %   sort = Sorte;
   % flag = Flag;
   %else
   % sort = Sorte;
   %end;

%------------- COEFFICIENT DU FILTRE -------------

return; 

