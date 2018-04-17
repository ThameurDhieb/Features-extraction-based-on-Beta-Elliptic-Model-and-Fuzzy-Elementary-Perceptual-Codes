function [reponse, X_rest, Y_rest, T_rest] = detection_segments_excedentaires_beta_chevauchees(param_grandeur, tableau_affectation_intervalles, points, X_init, Y_init, T_init)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction detection_segments_excedentaires_beta_chevauchees permet             %%%%%
%%%%% de détécter l'existence detraits excédentaires au extrimités de (début et fin)   %%%%%
%%%%% de la trajectoire.                                                               %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 reponse = 1;          % pas de segments excedentaires detectés
%reponse = 0;          % detection de segments excedentaires
pas = 0.01;
nombre_points_min = 18;

nbr_intervalles = size(tableau_affectation_intervalles, 1);


T_Deb = tableau_affectation_intervalles(:,1);
T_Fin = tableau_affectation_intervalles(:,2);
A1 = param_grandeur(:,1);
B1 = param_grandeur(:,2);
Longueur = param_grandeur(:,3);

Longueur_Segment = [];
for i = 1 : nbr_intervalles
    longueur_segment_i = sqrt( ( (A1(i))^2 ) + ( (B1(i))^2 ) )  ;
    longueur_i = Longueur(i);
    Longueur_Segment = [Longueur_Segment; (longueur_segment_i + longueur_i)/2];
end

if nbr_intervalles <= 4
   reponse = 1;
   X_rest = X_init;
   Y_rest = Y_init;
   T_rest = T_init;

elseif nbr_intervalles >= 5
   t_deb = T_Deb(1);
   t_fin = T_Fin(nbr_intervalles);
   
   K1 = tableau_affectation_intervalles(2,5);
   K2 = tableau_affectation_intervalles(2,6);
   if K1 > K2
      num_dernier_intervalle_considere_deb = 2;
   else
      num_dernier_intervalle_considere_deb = 1;
   end

   K1 = tableau_affectation_intervalles(nbr_intervalles - 1,5);
   K2 = tableau_affectation_intervalles(nbr_intervalles - 1,6);
   if K1 < K2
      num_premier_intervalle_considere_fin = nbr_intervalles - 1;
   else
      num_premier_intervalle_considere_fin = nbr_intervalles;
   end
   
   longueur_totale_intra = 0;
   nbr_seg_intra = 0;
   for i = num_dernier_intervalle_considere_deb + 1 : (num_premier_intervalle_considere_fin - 1)
       longueur_totale_intra = longueur_totale_intra + Longueur_Segment(i);
       nbr_seg_intra = nbr_seg_intra + 1;
   end
   longueur_moy = longueur_totale_intra / nbr_seg_intra;
   somme_longueur_deb = 0;
   for i = 1 : num_dernier_intervalle_considere_deb
       somme_longueur_deb = somme_longueur_deb + Longueur_Segment(i);
   end
   somme_longueur_fin = 0;
   for i = num_premier_intervalle_considere_fin : nbr_intervalles
       somme_longueur_fin = somme_longueur_fin + Longueur_Segment(i);
   end

   rapport_segment_deb = ( 2 * longueur_moy ) / somme_longueur_deb;
   rapport_segment_fin = ( 2 * longueur_moy ) / somme_longueur_fin;
   
   if rapport_segment_deb > 5.5%4.5%3
      reponse = 0;
      t_deb = T_Fin(num_dernier_intervalle_considere_deb) + pas;
   end
   if rapport_segment_fin > 6%5%4
      reponse = 0;
      t_fin = T_Deb(num_premier_intervalle_considere_fin) - pas;
   end

   if reponse == 0
       
      indice_vect_t_deb = find( T_init >= t_deb);
      indice_t_deb = indice_vect_t_deb(1);
      indice_vect_t_fin = find( T_init >= t_fin);
      indice_t_fin = indice_vect_t_fin(1);
      
      if indice_t_fin > indice_t_deb + nombre_points_min
         X_rest = [];
         Y_rest = [];
         T_rest = [];
      
         for jr = indice_t_deb : indice_t_fin
             X_rest = [X_rest; X_init(jr)];
             Y_rest = [Y_rest; Y_init(jr)];
             T_rest = [T_rest; T_init(jr)];
         end

      else
         reponse = 1;
         X_rest = X_init;
         Y_rest = Y_init;
         T_rest = T_init;
      end
   else
      reponse = 1;
      X_rest = X_init;
      Y_rest = Y_init;
      T_rest = T_init;
   end
      
end


        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%





