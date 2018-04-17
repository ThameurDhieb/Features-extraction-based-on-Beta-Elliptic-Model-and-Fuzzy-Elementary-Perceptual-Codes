function [matpseudomots, nbr_pseudomots, taille_pseudomot] = extract_pseudo_mot_p(a , longueur_nbr_points_min)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction extractcaracters permet de delimiter les données (x,y,p) de chaque prototype                %%%%%
%%%%% dans la matrice de données (a) représentant les inscriptions successives des pseudo-mots                %%%%%
%%%%% d'un même phrase ou proposition saisie sur la tablette                                                                      %%%%%
%%%%% Cette delimitation est effectuée en transformant (dévisant) la matrice 2Dim (a) (a=[(x,y,p)*nbr point]) %%%%%
%%%%% en une matrice de 3Dim (matpseudomots) (matpseudomots=[ [(x,y,p)*nbr point/prototype]*nbr prototype]    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


di = 0;
meh = 0;
 
for i = size(a,1) : -1 : 1
   
%    if ((a(i,3)) == 0)
    if ((a(i,1)) == 0) & ((a(i,2)) == 0)

       
        a(i,:) = [];
        if i > 1
%           if (a(i-1,3)) == 1
           if ((a(i-1,1)) ~= 0) | ((a(i-1,2)) ~= 0)
               break;   % fin boucle for
           end
        end
    end
end

data_G = a;


indicateur_premier_passage = 0;
indicateur_atente_debut_caractere = 1;
%longueur_nbr_points_min = 9;

k = 1; % indice du pseudo_mot courant
j = 1; % indice dans la nouvelle matrice 'matpseudomots'
taille_pseudomot = [];

size_a = size(a,1);
if size_a > 1
for i = 1 : size_a

%    if (a(i,3)) == 1
    if ( ((a(i,1)) ~= 0) | ((a(i,2)) ~= 0) )
                 
       matpseudomots(j,:,k) = a(i,:);
         
       if ( i < length(a) )
%          if ( j < longueur_nbr_points_min) & ( a(i+1,3) == 0 )
          if ( j < longueur_nbr_points_min) & ( (a(i+1,1) == 0)&(a(i+1,2) == 0) )
             matpseudomots(:,:,k) = [];
             k = k - 1;
          end
%          if ( j >= longueur_nbr_points_min) & ( a(i+1,3) == 0 )
          if ( j >= longueur_nbr_points_min) & ( (a(i+1,1) == 0)&(a(i+1,2) == 0) )
             taille_pseudomot = [ taille_pseudomot; k , j];
          end
       elseif ( i == length(a) )
          if ( j < longueur_nbr_points_min )
             matpseudomots(:,:,k) = [];
             k = k - 1;
          end
          if ( j >= longueur_nbr_points_min)
             taille_pseudomot = [ taille_pseudomot; k , j];
          end           
       end

       j = j + 1;
       indicateur_premier_passage = 1;
       indicateur_atente_debut_caractere = 0;

    else  
%       if(a(i,3) == 0)&(a(i+1,3) == 1)&(indicateur_premier_passage == 1)
       if( (a(i,1) == 0)&(a(i,2) == 0) )&( ((a(i+1,1)) ~= 0) | ((a(i+1,2)) ~= 0) )&(indicateur_premier_passage == 1)
         k = k + 1;
         j = 1;
         indicateur_atente_debut_caractere = 1;
       end

%     else
%         if ((a(i,3)) == 0)
%         a(i,:) = [];
%     end
    end

end
end

nbr_pseudomots = k;

% if size_a == 1
%    matpseudomots(1,:,1) = a(1,:);
%    nbr_pseudomots = 1;
%    taille_pseudomot = [1];
% elseif size_a == 0
if ( size_a == 0 ) | ( size_a == 1 )
   matpseudomots = [];
   nbr_pseudomots = 0;
   taille_pseudomot = [];

end
