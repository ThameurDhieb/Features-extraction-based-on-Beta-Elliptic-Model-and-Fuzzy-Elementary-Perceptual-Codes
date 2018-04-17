function [N_seg,data_mot] = copie_segmentation_tablet_2(data)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction copie_segmentation_tablet permet décomposer le signa enligne    %%%%%
%%%%% en segements séparés par une levée suivie par une posée du stylo            %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


levee1 = [];
levee = [];
posee = [];

Tlevee1 = [];
Plevee1 = [];
indice_seg = [];
X_seg = [];
Y_seg = [];
P_seg = [];
T_seg = [];
data_mot = [];



X = data(:,1);
Y = data(:,2);
P = data(:,3);
T = data(:,4);

%------------------ dimension du vecteur X
M = length(X);
%------------------ définir les "levee"s et les "posee"s du stylo
N_seg = 0;
for i = 2 : M
    
      if ( P(i-1) == 0 ) & ( P(i) ~= 0 )
         posee=[posee;data(i,1) data(i,2) P(i) T(i) i ];       
      elseif ( P(i-1) ~= 0 )& ( P(i) == 0 )        
         levee=[levee;data(i,1) data(i,2) P(i) T(i) i ];
         N_seg = N_seg + 1;
      end 
   
end

levee1 = levee;

%------------------reconstruction des données :segmentation
 j = 1;


 if ( N_seg >= 1 )
     
    for j = 1 : N_seg

      for i = 1 : M
          if ( i >= posee(j,5) ) & ( i < levee(j,5) )
            indice_seg(i,j) = i;
            X_seg(i,j) = data(i,1);
            Y_seg(i,j) = data(i,2);
            P_seg(i,j) = data(i,3);
            T_seg(i,j) = data(i,4);
          elseif  ( i < posee(j,5) )
            indice_seg(i,j) = 0;        % les points d'indice_seg = 0 à la colonne j des vecteurs X_seg et Y_seg
            X_seg(i,j) = posee(j,1);    % n'appartiennent pas au segment j
            Y_seg(i,j) = posee(j,2);
            P_seg(i,j) = posee(j,3);
            T_seg(i,j) = posee(j,4);
          else  %i>levee(j,5)          
            indice_seg(i,j) = 0;        % les points d'indice_seg = 0 à la colonne j des vecteurs X_seg et Y_seg
            X_seg(i,j) = levee(j,1);    % n'appartiennent pas au segment j
            Y_seg(i,j) = levee(j,2);
            P_seg(i,j) = levee(j,3);
            T_seg(i,j) = levee(j,4);
          end
      end
    end

 else
      j = 1
      for i = 1 : M
          if  ( i >= posee(j,5) )
          
            indice_seg(i,j) = i;
            X_seg(i,j) = data(i,1);
            Y_seg(i,j) = data(i,2);
            P_seg(i,j) = data(i,3);
            T_seg(i,j) = data(i,4);
          else
            indice_seg(i,j) = 0;         % les points d'indice_seg = 0 à la colonne j des vecteurs X_seg et Y_seg
            X_seg(i,j) = posee(j,1);     % n'appartiennent pas au segment j
            Y_seg(i,j) = posee(j,2);
            P_seg(i,j) = posee(j,3);
            T_seg(i,j) = posee(j,4);
           
          end
      end
    
 end

 
 data_mot = [indice_seg X_seg Y_seg P_seg T_seg];
 

 return 