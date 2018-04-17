function [data_norm, rapport_normalisation] = normalisation_H2_hauteur(data,n)
%function norm= normalisation_H(x,y,pression,n)

% methode norm_128
% cette fonction permet de normaliser un caract�re d�fini par les coordonn�s des points qui les constituent
% x:cha�ne des abscisses
% y:cha�ne des ordonn�es
% n:la longueur de la cha�ne x ou y, c'est le nombre des points d'un caract�re
% n=length (x);

X_pression = [];
Y_pression = [];

x = data(:,1);
y = data(:,2);
pression = data(:,3);

hauteur_norm = 102.5;    %128;
ordonnee_min = 25.5;     %20;

for i = 1 : n
    if pression(i) ~= 0
       X_pression = [X_pression; x(i)];
       Y_pression = [Y_pression; y(i)];
    end
end



maxi_x = max(X_pression);% maxi_x est le maximum des abscisses
maxi_y = max(Y_pression);% maxi_y est le maximum des ordonn�es
mini_x = min(X_pression);% mini_x est le minimum des abscisses
mini_y = min(Y_pression);% mini_y est le minimum des ordonn�es

%m = max(maxi_x - mini_x , maxi_y  - mini_y);
m = maxi_y  - mini_y; % Normalisation suivant la hauteur

ye = [];
xe = [];

if(m ~= 0)
   
   for i = 1 : n
      
      if (pression(i) ~= 0)
         xn_i = hauteur_norm*((x(i)-mini_x)/m);
         xe = [xe,xn_i];
         yn_i = (hauteur_norm*((y(i)-mini_y)/m)) + ordonnee_min;
         ye = [ye,yn_i];
      else
         xn_i = x(i);
         xe = [xe,xn_i];
         yn_i = y(i);
         ye = [ye,yn_i];
      end
 
   end

end

data_norm = []; % norm est une matrice de n lignes et de trois colonnes
                % la premi�re colonne repr�sente les abscisses et la deuxi�me repr�sente les ordonn�es
                % la troisi�me la pression
for i = 1 : n
   coord = [xe(i),ye(i)];
   data_norm = [data_norm; coord pression(i)];
end 


rapport_normalisation = hauteur_norm / m;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % % % % function [data_norm] = normalisation_H2_hauteur(data,n)
% % % % % %function norm= normalisation_H(x,y,pression,n)
% % % % % 
% % % % % % methode norm_128
% % % % % % cette fonction permet de normaliser un caract�re d�fini par les coordonn�s des points qui les constituent
% % % % % % x:cha�ne des abscisses
% % % % % % y:cha�ne des ordonn�es
% % % % % % n:la longueur de la cha�ne x ou y, c'est le nombre des points d'un caract�re
% % % % % % n=length (x);
% % % % % 
% % % % % rapport_norm = 128;
% % % % % 
% % % % % X = data(:,1);
% % % % % Y = data(:,2);
% % % % %        
% % % % %        
% % % % % maxi_x = max(X);% maxi_x est le maximum des abscisses
% % % % % maxi_y = max(Y);% maxi_y est le maximum des ordonn�es
% % % % % mini_x = min(X);% mini_x est le minimum des abscisses
% % % % % mini_y = min(Y);% mini_y est le minimum des ordonn�es
% % % % % 
% % % % % %m=max(maxi_x - mini_x , maxi_y  - mini_y);
% % % % % m = maxi_y  - mini_y;  % Normalisation suivant la hauteur
% % % % % 
% % % % % xs = [];
% % % % % ys = [];
% % % % % 
% % % % % if(m ~= 0)
% % % % %    
% % % % %    for i = 1 : n
% % % % %       
% % % % %       xn(i) = rapport_norm * ((X(i)-mini_x)/m);
% % % % %       xs = [xs; xn(i)];
% % % % %       yn(i) = rapport_norm * ((Y(i)-mini_y)/m);
% % % % %       ys = [ys; yn(i)];
% % % % %  
% % % % %    end
% % % % % else
% % % % %    xs = X;
% % % % %    ys = Y;    
% % % % % end 
% % % % % 
% % % % % % data_norm est une matrice de n lignes et de deux colonnes
% % % % % % la premi�re colonne repr�sente les abscisses et la deuxi�me repr�sente les ordonn�es
% % % % % 
% % % % % 
% % % % % data_norm = [xs , ys];
