function [data_prim] = reparation_de_l_echantillonnage_faible_des_pseudo_mots( data )

%%%%HHHHHHHHHHHHHH%BBBBBBBBBB%%%%%%
%%%%HHHHHHHHHHHHHH%BBBBBBBBBB%%%%%%
%%%%HHHHHHHHHHHHHH%BBBBBBBBBB%%%%%%

  %%%%% initialisation par défaut
  nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture = 0;
  nbr_de_points_dans_le_pseudo_mot_en_cours = 0;
  
  parcourt_d_un_pseudo_mot = 0;
  passage_d_un_pseudo_mot_a_un_autre = 1;
  
  vect_indices_points_limites_des_pseudo_mots = [];

  for k = 1 : size(data,1)
      x_data = data(k,1);
      y_data = data(k,2);

      indice_point_courant = k;
      if ( x_data ~= 0 ) | ( y_data ~= 0 )
         x_point_data = x_data;
         y_point_data = y_data;
         
         nbr_de_points_dans_le_pseudo_mot_en_cours = nbr_de_points_dans_le_pseudo_mot_en_cours + 1;
         if ( passage_d_un_pseudo_mot_a_un_autre == 1 ) % & (parcourt_d_un_pseudo_mot == 0)
            nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture = nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture + 1;
            indice_premier_point_du_pseudo_mot_en_cours = indice_point_courant;
            
            nbr_de_points_dans_le_pseudo_mot_en_cours = 1;
         end

         parcourt_d_un_pseudo_mot = 1;
         passage_d_un_pseudo_mot_a_un_autre = 0;
         indice_dernier_point_du_pseudo_mot_en_cours = indice_point_courant;
      else
         if ( parcourt_d_un_pseudo_mot == 1 ) & (indice_point_courant > 1)
            parcourt_d_un_pseudo_mot = 0;
            passage_d_un_pseudo_mot_a_un_autre = 1;
            
            %indice_dernier_point_du_pseudo_mot_en_cours = indice_point_courant - 1;
            
            vect_indices_points_limites_des_pseudo_mots = [vect_indices_points_limites_des_pseudo_mots; indice_premier_point_du_pseudo_mot_en_cours , indice_dernier_point_du_pseudo_mot_en_cours];
         end
         
      end
  end

  
  
  data_prim = [];
  
  taille_minimale_pseudo_mot = 7;
  taille_minimale_points_reechantillones = 4;
  
  nbr_de_pseudo_mots_dans_la_ligne_d_ecriture = size(vect_indices_points_limites_des_pseudo_mots , 1);
  for k = 1 : nbr_de_pseudo_mots_dans_la_ligne_d_ecriture
      indice_premier_point_du_pseudo_mot_en_cours = vect_indices_points_limites_des_pseudo_mots(k,1);
      indice_dernier_point_du_pseudo_mot_en_cours = vect_indices_points_limites_des_pseudo_mots(k,2);
      nbr_de_points_dans_pseudo_mot = max( ( indice_dernier_point_du_pseudo_mot_en_cours - indice_premier_point_du_pseudo_mot_en_cours + 1) , 0 );
      
      if ( nbr_de_points_dans_pseudo_mot <= taille_minimale_pseudo_mot ) & ( nbr_de_points_dans_pseudo_mot >= taille_minimale_points_reechantillones - 1 )
         troncon_data_prim_pseudo_mot = [];
         if ( nbr_de_points_dans_pseudo_mot > taille_minimale_points_reechantillones )
            for k_hh = indice_premier_point_du_pseudo_mot_en_cours : (indice_dernier_point_du_pseudo_mot_en_cours - 1)
                x_point_data_Mk = data(k_hh , 1);
                y_point_data_Mk = data(k_hh , 2);
                x_point_data_Mkplus1 = data(k_hh + 1 , 1);
                y_point_data_Mkplus1 = data(k_hh + 1 , 2);
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; x_point_data_Mk , y_point_data_Mk];
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; ( ( x_point_data_Mk + x_point_data_Mkplus1 ) / 2 ) , ( ( y_point_data_Mk + y_point_data_Mkplus1 ) / 2 ) ];
            end
         elseif ( nbr_de_points_dans_pseudo_mot == taille_minimale_points_reechantillones )
            for k_hh = indice_premier_point_du_pseudo_mot_en_cours : (indice_dernier_point_du_pseudo_mot_en_cours - 1)
                x_point_data_Mk = data(k_hh , 1);
                y_point_data_Mk = data(k_hh , 2);
                x_point_data_Mkplus1 = data(k_hh + 1 , 1);
                y_point_data_Mkplus1 = data(k_hh + 1 , 2);
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; x_point_data_Mk , y_point_data_Mk];
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; ( ( (1 * x_point_data_Mk) + (2 * x_point_data_Mkplus1) ) / 3 ) , ( ( (1 * y_point_data_Mk) + (2 * y_point_data_Mkplus1) ) / 3 ) ];
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; ( ( (2 * x_point_data_Mk) + (1 * x_point_data_Mkplus1) ) / 3 ) , ( ( (2 * y_point_data_Mk) + (1 * y_point_data_Mkplus1) ) / 3 ) ];
            end
         elseif ( nbr_de_points_dans_pseudo_mot <= taille_minimale_points_reechantillones - 1)
            for k_hh = indice_premier_point_du_pseudo_mot_en_cours : (indice_dernier_point_du_pseudo_mot_en_cours - 1)
                x_point_data_Mk = data(k_hh , 1);
                y_point_data_Mk = data(k_hh , 2);
                x_point_data_Mkplus1 = data(k_hh + 1 , 1);
                y_point_data_Mkplus1 = data(k_hh + 1 , 2);
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; x_point_data_Mk , y_point_data_Mk];
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; ( ( (1 * x_point_data_Mk) + (3 * x_point_data_Mkplus1) ) / 4 ) , ( ( (1 * y_point_data_Mk) + (3 * y_point_data_Mkplus1) ) / 4 ) ];
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; ( ( (2 * x_point_data_Mk) + (2 * x_point_data_Mkplus1) ) / 4 ) , ( ( (2 * y_point_data_Mk) + (2 * y_point_data_Mkplus1) ) / 4 ) ];
                troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot; ( ( (3 * x_point_data_Mk) + (1 * x_point_data_Mkplus1) ) / 4 ) , ( ( (3 * y_point_data_Mk) + (1 * y_point_data_Mkplus1) ) / 4 ) ];
            end
         end
         troncon_data_prim_pseudo_mot = [troncon_data_prim_pseudo_mot;  data( indice_dernier_point_du_pseudo_mot_en_cours , :)];

         data_prim = [data_prim; 0, 0; troncon_data_prim_pseudo_mot; 0, 0];

      else
         
         data_prim = [data_prim; 0, 0; data( indice_premier_point_du_pseudo_mot_en_cours : indice_dernier_point_du_pseudo_mot_en_cours , : ); 0, 0];
         
      end
      

  end
  
%%%%HHHHHHHHHHHHHH%BBBBBBBBBB%%%%%%
%%%%HHHHHHHHHHHHHH%BBBBBBBBBB%%%%%%
%%%%HHHHHHHHHHHHHH%BBBBBBBBBB%%%%%%


