function [Err_moy_d, Err_Max_d, Err_moy_p, Err_Max_p] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction mesure_erreur_reconstruction_arc permet de calculer                                %%%%%
%%%%% les erreurs de reconstruction d'un segment de la trajectoire par un                            %%%%%
%%%%% arc d'ellipse défini par les parmètres de son équation. le calcul                              %%%%%
%%%%% de l'erreur repose sur la mesure de différence entre les rayons de                             %%%%%
%%%%% divergant du centre de l'ellipse et atteignant l'arc et ceux                                   %%%%%
%%%%% divergant du même centre et atteignant la trajectoire.                                         %%%%%
%%%%% quatre types d'erreurs sont fournies : erreur moyenne en distance, erreur maximale en distance,%%%%%
%%%%% erreur moyenne en pourcent et erreur maximale en pourcent                                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nbr_points = size(liste_points,1);

a_quart = ortho_dist(1);
b_quart = ortho_dist(2);

rapport_grandeur = sqrt( (a_quart^2) + (b_quart^2) );

somme_dist_Erreur = 0;
dist_Max_Erreur = 0;
if nbr_points >= 3
   for i = 2 : (nbr_points - 1)
       xM = liste_points(i,1);
       yM = liste_points(i,2);
       [dist_Erreur] = mesure_erreur_reconstruction_point(xM,yM,parametre_ellipse);
       somme_dist_Erreur = somme_dist_Erreur + dist_Erreur;
       if dist_Erreur > dist_Max_Erreur
          dist_Max_Erreur = dist_Erreur;
       end
   end
   
   Err_moy_d = somme_dist_Erreur / (nbr_points);
   Err_Max_d = dist_Max_Erreur;

   Err_moy_p = ( Err_moy_d / rapport_grandeur ) * 100;
   Err_Max_p = ( Err_Max_d / rapport_grandeur ) * 100;
else
   Err_moy_d = 0;
   Err_Max_d = 0;

   Err_moy_p = 0;
   Err_Max_p = 0;
end

        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
