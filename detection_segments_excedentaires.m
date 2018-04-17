function [reponse, X_rest, Y_rest, T_rest] = detection_segments_excedentaires(param_grandeur, IMPULSION_BETA, points, X_init, Y_init, T_init);



 reponse = 1;          % pas de segments excedentaires detectés
%reponse = 0;          % detection de segments excedentaires
pas = 0.01;
nombre_points_min = 18;

nbr_seg = size(IMPULSION_BETA, 1);


T0 = IMPULSION_BETA(:,6);
T1 = IMPULSION_BETA(:,7);
A1 = param_grandeur(:,1);
B1 = param_grandeur(:,2);
A2 = param_grandeur(:,3);
B2 = param_grandeur(:,4);
Longueur = param_grandeur(:,5);

Longueur_Segment = [];
for i = 1 : nbr_seg
    longueur_segment_i = ( sqrt( ( (A1(i))^2 ) + ( (B1(i))^2 ) ) ) + (  sqrt( ( (A2(i))^2 ) + ( (B2(i))^2 ) ) );
    longueur_i = Longueur(i);
    Longueur_Segment = [Longueur_Segment; (longueur_segment_i + longueur_i)/2];
end

if nbr_seg <= 2
   reponse = 1;
   X_rest = X_init;
   Y_rest = Y_init;
   T_rest = T_init;

elseif nbr_seg >= 3
   t_deb = T0(1);
   t_fin = T1(nbr_seg);
   longueur_totale_intra = 0;
   for i = 2 : (nbr_seg - 1)
        longueur_totale_intra = longueur_totale_intra + Longueur_Segment(i);
   end
   longueur_moy = longueur_totale_intra / (nbr_seg - 2);
   rapport_segment_deb = longueur_moy / Longueur_Segment(1);
   rapport_segment_fin = longueur_moy / Longueur_Segment(nbr_seg);
   

   if rapport_segment_deb > 5.5%4.5%3
      reponse = 0;
      t_deb = T1(1) + pas;
   end
   if rapport_segment_fin > 6%5%4
      reponse = 0;
      t_fin = T0(nbr_seg) - pas;
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
