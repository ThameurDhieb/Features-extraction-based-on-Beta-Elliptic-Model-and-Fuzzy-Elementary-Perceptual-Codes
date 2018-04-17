%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -SCRIPT    : extraction des caracteristiques par l'approche beta_eliptique             %
% -SYNTAX    : Lecture_traitement_save_param_phrase_online                               %
% -ENTREE    : chemin d'une suite de fichiers de mots dont l'ordre temporel de leurs     %
%              trajectoire off line est reconstruit                                      %
% -SORTIES   : fichier texte contenant les cractéristiques extraites des graphèmes       %
%              segmentés                                                                 %
% -copyright : REGIM-ENIS                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mmatrice_param_1, mmatrice_param_2] = Lecture_traitement_phrase_online_beta_elliptique_preclass_2(Data, tline, numero_mot, numero_scripteur, chemin_acces_dossier_Base_Param_I_formes_visuelles, chemin_acces_dossier_Base_Param_I_formes_visuelles_FEPC,  chemin_acces_dossier_Base_Param_II_formes_visuelles, chemin_acces_dossier_Base_Param_II_formes_visuelles_FEPC, nom_du_fichier_echantillon_sans_extension, elimin_trait_excedent, Nombre_de_traits_double_d_une_forme_visuelle, Nombre_de_traits_chevauches_d_une_forme_visuelle)



rayon_filtre_V = 2;              %%% Valeur initiale
sigma_p_filtre_V = 0.7;          %%% Valeur initiale

rayon_filtre_traject = 3;        %%% Valeur initiale
sigma_p_filtre_traject = 3.5;    %%% Valeur initiale



%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  
%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  
%%%%%--------------Décomposition d'un paragraphe manuscrit Online en lignes d'écriture manuscrite -------------%%%%%
%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  
%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  

  max_y_Data = max( Data(:,2) );
  figure(1);
  axis equal;
  hold on;
  plot(Data(:,1), Data(:,2),'.b');
  
  
  [Data_repare] = reparation_de_l_echantillonnage_faible_des_pseudo_mots( Data );
  Data = Data_repare;

  
  max_y_Data_brut_init = max_y_Data;

  rap_min_longueur_hauteur_ligne_ecriture_complete = 5;


%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Initialisation des paramètres de la Décomposition en lignes de texte %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  indice_points_de_debut_de_la_Ligne_d_ecriture_manuscrite = 1;
  indice_points_de_fin_de_la_Ligne_d_ecriture_manuscrite = size(Data,1);
  
  indice_points_de_debut_du_pseudo_mot_courant = 1;
  indice_points_de_fin_du_pseudo_mot_courant = size(Data,1);
  
  Vecteur_des_Lignes_d_ecriture_manuscrite = [indice_points_de_debut_de_la_Ligne_d_ecriture_manuscrite , indice_points_de_fin_de_la_Ligne_d_ecriture_manuscrite];
  vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt = [indice_points_de_debut_de_la_Ligne_d_ecriture_manuscrite, indice_points_de_fin_de_la_Ligne_d_ecriture_manuscrite];
  Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite = [];
  Vect_taille_horiz_et_vertec_des_lignes_d_ecrit_manuscrit = [];
  Pure_Vect_Points_ligne_ecriture = [];
  nbr_de_points_pures_inclus = 0;
  nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture = 0;
  parcourt_d_un_pseudo_mot = 0;
  passage_d_un_pseudo_mot_a_un_autre = 1;
  changement_ligne_ecriture_manuscrite = 0;
  
  indice_dernier_point_pseudo_mot_preced = 1;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%



%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% recherche des indices des points limites pour la decomposition %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% de la suite de pseudo_mots saisis en ligne d'écriture manuscrite %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for k = 1 : size(Data,1)
      x_Data = Data(k,1);
      y_Data = Data(k,2);
      
      if ( x_Data ~= 0 ) | ( y_Data ~= 0 )
         x_point_Data = x_Data;
         y_point_Data = y_Data;

         indice_point_courant = k;

         if passage_d_un_pseudo_mot_a_un_autre == 1;
            nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture = nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture + 1;
            vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(end, 2) = indice_dernier_point_pseudo_mot_preced;
            vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt = [vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt; indice_point_courant , size(Data,1)];
         end
         
         changement_ligne_ecriture_manuscrite = 0;
         
         if (passage_d_un_pseudo_mot_a_un_autre == 1) & ( nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture > 3 ) %& (parcourt_d_un_pseudo_mot == 0)
            dist_horiz_entre_deb_pseudo_courant_fin_pseudo_preced = abs( x_point_Data - x_fin_pseudo_mot_preced );
            
            if ( (dist_horiz_entre_deb_pseudo_courant_fin_pseudo_preced / hauteur_ligne_d_ecriture) > rap_min_longueur_hauteur_ligne_ecriture_complete )
               changement_ligne_ecriture_manuscrite = 1;
               Vecteur_des_Lignes_d_ecriture_manuscrite(end, 2) = indice_dernier_point_pseudo_mot_preced;
               Vecteur_des_Lignes_d_ecriture_manuscrite = [Vecteur_des_Lignes_d_ecriture_manuscrite; indice_point_courant , size(Data,1)];

               Vect_taille_horiz_et_vertec_des_lignes_d_ecrit_manuscrit = [Vect_taille_horiz_et_vertec_des_lignes_d_ecrit_manuscrit; hauteur_ligne_d_ecriture , longueur_ligne_d_ecriture];

               %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt = size( vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt , 1);
               vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(1:(size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt - 1) , :);
               size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt = size( vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt , 1);
                           
               if ( (size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt - 1) >= 10 )
                  numero_pseudo_mot_milieu_de_la_ligne_de_texte = fix ( size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt / 2 );
                  
                  indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(1,1);
                  indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(numero_pseudo_mot_milieu_de_la_ligne_de_texte,2);
                  indices_points_deb_de_la_deuxieme_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(numero_pseudo_mot_milieu_de_la_ligne_de_texte + 1,1);
                  indices_points_fin_de_la_deuxieme_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt,2);
                  indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt = [indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt , indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt];
                  indices_points_deb_et_fin_de_la_deuxieme_moitier_dela_ligne_txt = [indices_points_deb_de_la_deuxieme_moitier_de_la_ligne_txt , indices_points_fin_de_la_deuxieme_moitier_de_la_ligne_txt];                
                  Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite = [Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite; indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt; indices_points_deb_et_fin_de_la_deuxieme_moitier_dela_ligne_txt];
                 
               else
                  indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(1,1);
                  indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt( size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt ,2);
                  indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt = [indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt , indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt];
                  Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite = [Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite; indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt];

               end
                       
               vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt = [indice_point_courant , size(Data,1)];
               nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture = 1;
               Pure_Vect_Points_ligne_ecriture = [];
            end
         end

         Pure_Vect_Points_ligne_ecriture = [Pure_Vect_Points_ligne_ecriture; x_point_Data , y_point_Data];
         nbr_de_points_pures_inclus = nbr_de_points_pures_inclus + 1;

         parcourt_d_un_pseudo_mot = 1;
         passage_d_un_pseudo_mot_a_un_autre = 0;
      else
         if ( parcourt_d_un_pseudo_mot == 1 )
            parcourt_d_un_pseudo_mot = 0;
            passage_d_un_pseudo_mot_a_un_autre = 1;
            
            x_fin_pseudo_mot_preced = x_point_Data;
            indice_dernier_point_pseudo_mot_preced = indice_point_courant;
         end
         if ( nbr_de_pseudo_mots_parcourus_dans_la_ligne_d_ecriture >= 1 ) %>= 2 )
            max_y_pure_points = max(Pure_Vect_Points_ligne_ecriture(:,2));
            max_x_pure_points = max(Pure_Vect_Points_ligne_ecriture(:,1));
            min_y_pure_points = min(Pure_Vect_Points_ligne_ecriture(:,2));
            min_x_pure_points = min(Pure_Vect_Points_ligne_ecriture(:,1));
            hauteur_ligne_d_ecriture = abs( max_y_pure_points - min_y_pure_points );
            longueur_ligne_d_ecriture = abs( max_x_pure_points - min_x_pure_points );
         end
         
      end
  end

  Vect_taille_horiz_et_vertec_des_lignes_d_ecrit_manuscrit = [Vect_taille_horiz_et_vertec_des_lignes_d_ecrit_manuscrit; hauteur_ligne_d_ecriture , longueur_ligne_d_ecriture];

  size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt = size( vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt , 1);
  if ( (size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt - 1) >= 7 )
     numero_pseudo_mot_milieu_de_la_ligne_de_texte = fix ( size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt / 2 );
     indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(1,1);
     indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(numero_pseudo_mot_milieu_de_la_ligne_de_texte,2);
     indices_points_deb_de_la_deuxieme_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(numero_pseudo_mot_milieu_de_la_ligne_de_texte + 1,1);
     indices_points_fin_de_la_deuxieme_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt,2);
     indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt = [indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt , indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt];
     indices_points_deb_et_fin_de_la_deuxieme_moitier_dela_ligne_txt = [indices_points_deb_de_la_deuxieme_moitier_de_la_ligne_txt , indices_points_fin_de_la_deuxieme_moitier_de_la_ligne_txt];
     Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite = [Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite; indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt; indices_points_deb_et_fin_de_la_deuxieme_moitier_dela_ligne_txt];
  else
     indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt(1,1);
     indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt = vect_indices_pnts_deb_et_fin_paws_inclus_dans_ligne_txt( size_vect_indis_pnts_deb_et_fin_paws_inclu_ds_line_txt ,2);
     indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt = [indices_points_deb_de_la_premiere_moitier_de_la_ligne_txt , indices_points_fin_de_la_premiere_moitier_de_la_ligne_txt];
     Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite = [Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite; indices_points_deb_et_fin_de_la_premiere_moitier_dela_ligne_txt];
  end
  
nbr_de_lignes_d_ecriture_manuscrite = size( Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite , 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%



%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  
%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  
%%%%%------------  Début de la boucle de lecture des données d'écriture manuscrite ligne par ligne   ----------%%%%%
%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  
%%%%%---------********---------^^^^^^^^^^----------:::::::::::::---------^^^^^^^^^^---------*********----------%%%%%  

for indice_ligne_ecriture = 1 : nbr_de_lignes_d_ecriture_manuscrite

%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% lecture de la matrice des points successifs de la trajectoire %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%  d'écriture manuscrite d'une ligne de texte au plus, suivant  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%           le mode de fonctionnment interne défini             %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  indice_point_deb_ligne_ecriture_courante = Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite(indice_ligne_ecriture , 1);
  indice_point_fin_ligne_ecriture_courante = Vecteur_des_troncons_de_lignes_d_ecriture_manuscrite(indice_ligne_ecriture , 2);
  data = Data(indice_point_deb_ligne_ecriture_courante : indice_point_fin_ligne_ecriture_courante ,:);
  data = [data; 0 , 0; 0 , 0];
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&%%%%%
  a = data; 
  x = a(:,1);
  y = a(:,2);
  max_y_data = max(data(:,2));
  max_x_data = max(data(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%                    Normalisation en hauteur du mot                %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  points_equidistants = [];
  Pression_levee_point = 0;
  Pression_bessee_point = 1;

  for k = 1 : size(data,1)
      if ( data(k,1) ~= 0 ) | ( data(k,2) ~= 0 )
         points_equidistants = [points_equidistants; data(k,1), data(k,2), Pression_bessee_point];
      else
         points_equidistants = [points_equidistants; data(k,1), data(k,2), Pression_levee_point];
      end
  end

  [points_equidistants]= normalisation_H2_hauteur(points_equidistants, size(points_equidistants,1));
  points_equidistants = points_equidistants(:,1:2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  longueur_nbr_points_min = 1;
  [mat_pseudo_mot,nbr_pseudo_mot,taille_pseudo_mot] = extract_pseudo_mot_p(points_equidistants , longueur_nbr_points_min);


  taille_minimale_pseudo_mot = 7;
  taille_minimale_points_filtres_avec_pression = 4;
  

  mmatrice_param_1 = [];
  mmatrice_param_2 = [];

  points_total = [];
  points_total_filtre = [];
  DAT_par_pseudo_mot_total = [];
  nbr_points_par_pseudo_mot = [];
  Ensemble_max_min_x_y = [];


if nbr_pseudo_mot > 0
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   Filtrage et calcul de l'angle d'inclinaison de la tangent DAT   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  k_p = 0;
  for k = 1 : nbr_pseudo_mot

      data_k = mat_pseudo_mot(:,:,k);
      data_k = data_k(1 : taille_pseudo_mot(k,2) ,:);
      M = size(data_k,1);
      if M > taille_minimale_pseudo_mot
         
         k_p = k_p + 1;
         numero_pseudo_mot = k_p;
         
         [max_min_x_y, DAT_par_pseudo_mot, points_filtre] = Filtrage_conventionnel_plus_Dat_et_Detection_lim_max_min(data_k, rayon_filtre_traject, sigma_p_filtre_traject);
         
         Num_pseudo_mot = [];
         for i = 1 : M
             Num_pseudo_mot = [Num_pseudo_mot; k_p];
         end

         DAT_par_pseudo_mot_total = [DAT_par_pseudo_mot_total; DAT_par_pseudo_mot, Num_pseudo_mot];
         Ensemble_max_min_x_y = [Ensemble_max_min_x_y; max_min_x_y];
         points_total = [points_total; data_k, Num_pseudo_mot];
         points_total_filtre = [points_total_filtre; points_filtre, Num_pseudo_mot];
         
         nbr_points_par_pseudo_mot = [nbr_points_par_pseudo_mot; numero_pseudo_mot M];
      end    

  end

  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  ajout de la pression du stylo et constitution de la matrice des  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%               de points pseudo mots à trois indices               %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%         (numéro pseudomot , x y pression , noméro point)          %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%       avec un indice de longueur variable : numero points         %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  taille_de_tous_les_pseudo_mots_en_nbr_de_points = [];
  pression_baissee = 1;
  pression_levee   = 0;
  k_p = 0;
  for k = 1 : nbr_pseudo_mot
      numero_pseudo_mot = k;
%      pause
      data_k = mat_pseudo_mot(:,:,k);
      data_k = data_k(1 : taille_pseudo_mot(k,2) ,:);
      M = size(data_k,1);
      points_filtres_pression_pseudo_mot = [];
      if M > taille_minimale_pseudo_mot
         k_p = k_p + 1;
      end
      for i = 1 : M
          points_filtres_pression_pseudo_mot = [points_filtres_pression_pseudo_mot; data_k(i,:), pression_baissee];
      end
      points_filtres_pression_pseudo_mot = [points_filtres_pression_pseudo_mot(1,1) points_filtres_pression_pseudo_mot(1,2) pression_levee; points_filtres_pression_pseudo_mot; points_filtres_pression_pseudo_mot(M,1) points_filtres_pression_pseudo_mot(M,2) pression_levee]; 
      data_k_filtree = data_k;        
      for indice_point = 1 : size(points_filtres_pression_pseudo_mot,1)
          mat_pseudo_mot_de_points_filtres_avec_pression(indice_point,:,numero_pseudo_mot) = points_filtres_pression_pseudo_mot(indice_point,:);
      end
      for indice_point = 1 : size(data_k_filtree,1)
          mat_data_k_filtree_pseudo_mot(indice_point,:,numero_pseudo_mot) = data_k_filtree(indice_point,:);
      end
      taille_de_tous_les_pseudo_mots_en_nbr_de_points = [taille_de_tous_les_pseudo_mots_en_nbr_de_points; M , size(points_filtres_pression_pseudo_mot,1)];
  end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  Traitement des pseudo - mots pour leurs segmentations en traits  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      beta - elliptiques et l'extraction de leurs paramètres       %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ representation graphique ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
  nombre_totale_de_graphemes = 0;
  nbr_pseudo_mot_traite = 0;
  k_p = 0;
  k_q = 0; 
  figure (51);
  subplot(2,1,1);
  hold off;
  subplot(2,1,2);
  hold off;
  for k = 1 : nbr_pseudo_mot
      numero_pseudo_mot = k;
      points_filtres_avec_pression = mat_pseudo_mot_de_points_filtres_avec_pression(:,:,k);
      points_filtres_avec_pression = points_filtres_avec_pression(1 : taille_de_tous_les_pseudo_mots_en_nbr_de_points(k,2) , :);
      data_k_filtree = mat_data_k_filtree_pseudo_mot(:,:,k);
      data_k_filtree = data_k_filtree(1 : taille_de_tous_les_pseudo_mots_en_nbr_de_points(k,1) , :);
      M = taille_de_tous_les_pseudo_mots_en_nbr_de_points(k,1);
      
      if M > taille_minimale_pseudo_mot
         k_p = k_p + 1;
         
         figure(51);
         subplot(2,1,1);
         axis equal;
         plot(points_filtres_avec_pression(:,1),points_filtres_avec_pression(:,2),'.b');
         hold on;
         
      end
%       pause
  end

%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%
%^^^ Traitement des pseudo - mots : segmentation en traits Bêta - elliptique + extraction des caractéristiques ^^^%
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%

  nombre_totale_de_graphemes = 0;
  nbr_pseudo_mot_traite = 0;
  k_p = 0;
  k_q = 0;
  
  for k = 1 : nbr_pseudo_mot
      numero_pseudo_mot = k;
%      pause
      points_filtres_avec_pression = mat_pseudo_mot_de_points_filtres_avec_pression(:,:,k);
      points_filtres_avec_pression = points_filtres_avec_pression(1 : taille_de_tous_les_pseudo_mots_en_nbr_de_points(k,2) , :);
      data_k_filtree = mat_data_k_filtree_pseudo_mot(:,:,k);
      data_k_filtree = data_k_filtree(1 : taille_de_tous_les_pseudo_mots_en_nbr_de_points(k,1) , :);
      M = taille_de_tous_les_pseudo_mots_en_nbr_de_points(k,1);
      
      if M > taille_minimale_pseudo_mot
         k_p = k_p + 1;
         
         if ( size(points_filtres_avec_pression,1) > taille_minimale_points_filtres_avec_pression )  
            nbr_pseudo_mot_traite = nbr_pseudo_mot_traite + 1;
            
            [vect_param, matrice_param, param_trajectoire_BC, DAT_BC, indiceTemps_deb, points, nbr_trait, e, Vect_Erreurs, famille_tg_horiz, KMH_deb_fin_traits_chevauches] = new_beta_ellipse_stretegy_Beta_chevauchees_preclass(points_filtres_avec_pression, elimin_trait_excedent, rayon_filtre_V, sigma_p_filtre_V);
            vect_param_1 = vect_param;
            
            mmatrice_param_1 = [mmatrice_param_1, vect_param_1];


            
            [vect_param, matrice_param, param_trajectoire_BC, DAT_BC, indiceTemps_deb, points, nbr_trait, e, Vect_Erreurs_BC, KMH_deb_fin_traits_double] = new_beta_ellipt_stretegy_Beta_compsnt_entrnmnt_preclass(points_filtres_avec_pression, elimin_trait_excedent, rayon_filtre_V, sigma_p_filtre_V);
            vect_param_2 = vect_param;

            mmatrice_param_2 = [mmatrice_param_2, vect_param_2];
            figure (51);
            subplot(2,1,2);
            axis equal;
            plot(points(:,1), points(:,2), 'Color',[.3 .7 .999], 'lineStyle', '-.');%'.b');
            hold on;
            
            %%% pause(0.05);
            
            nbr_traits_methode_simplifiee = size(vect_param_2 , 2);
            nbr_traits_methode_chevauchee = size(vect_param_1 , 2);
            indice_traits_deb = 1;
            indice_traits_arrivee_preced = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            while (indice_traits_arrivee_preced < nbr_traits_methode_chevauchee)
                KMH_deb_troncon_trajectoire_a_enregistrer = KMH_deb_fin_traits_chevauches( indice_traits_deb , 1 );
                KMH_fin_troncon_trajectoire_a_enregistrer = min(   (size(points,1))    ,    ( KMH_deb_fin_traits_chevauches( (min( nbr_traits_methode_chevauchee , (indice_traits_deb + Nombre_de_traits_chevauches_d_une_forme_visuelle - 1) )) , 2 ) )    );
                nombre_de_traits_considres = Nombre_de_traits_chevauches_d_une_forme_visuelle - (    (indice_traits_deb + Nombre_de_traits_chevauches_d_une_forme_visuelle - 1)     -     ( min( nbr_traits_methode_chevauchee , (indice_traits_deb + Nombre_de_traits_chevauches_d_une_forme_visuelle - 1) ) )     );
                nbr_traits_methode_consideree = nbr_traits_methode_chevauchee;
                
                indice_trait_deb_methode_simplifiee = 0;
                iiiiih = 1;
                while (iiiiih <= nbr_traits_methode_simplifiee) & ( indice_trait_deb_methode_simplifiee == 0 )
                    KMH_deb_troncon_traits_simplifie_iiiiih = KMH_deb_fin_traits_double( iiiiih , 1 );
                    if ( KMH_deb_troncon_traits_simplifie_iiiiih == KMH_deb_troncon_trajectoire_a_enregistrer )
                       indice_trait_deb_methode_simplifiee = iiiiih;
                    end
                    iiiiih = iiiiih + 1;
                end
                if ( indice_trait_deb_methode_simplifiee == 0 )
                   indice_trait_deb_methode_simplifiee = 1;
                   iiiiih = 1;
                   while (iiiiih <= nbr_traits_methode_simplifiee) & ( indice_trait_deb_methode_simplifiee == 1 )
                       KMH_deb_troncon_traits_simplifie_iiiiih = KMH_deb_fin_traits_double( iiiiih , 1 );
                       if ( abs( KMH_deb_troncon_traits_simplifie_iiiiih - KMH_deb_troncon_trajectoire_a_enregistrer ) <= 1 )
                          indice_trait_deb_methode_simplifiee = iiiiih;
                       end
                       iiiiih = iiiiih + 1;
                   end
                end
                nombre_de_traits_considres_methode_simplifiee = Nombre_de_traits_double_d_une_forme_visuelle - (    (indice_trait_deb_methode_simplifiee + Nombre_de_traits_double_d_une_forme_visuelle - 1)     -     ( min( nbr_traits_methode_simplifiee , (indice_trait_deb_methode_simplifiee + Nombre_de_traits_double_d_une_forme_visuelle - 1) ) )     );
                indice_traits_arrivee_methode_simplifiee = indice_trait_deb_methode_simplifiee + nombre_de_traits_considres_methode_simplifiee - 1;
                KMH_deb_troncon_traits_simplifies_debut = KMH_deb_fin_traits_double( indice_trait_deb_methode_simplifiee , 1 );
                KMH_fin_troncon_traits_simplifies_arrivee = KMH_deb_fin_traits_double( indice_traits_arrivee_methode_simplifiee , 2 );
                indice_traits_arrivee = indice_traits_deb + nombre_de_traits_considres - 1;               
                Vect_x_y_points_limites_des_traits = [];
                Vect_angle_courbure_variation = [];
                Vect_delta_x_point_deb_point_fin_trait_courant = [];
                Vect_delta_y_point_deb_point_fin_trait_courant = [];
                Vect_delta_x_point_fin_point_deb_trait_courant = [];
                Vect_delta_y_point_fin_point_deb_trait_courant = [];

                for iiiiih = indice_traits_deb : indice_traits_arrivee
                    KMH_deb_trait_chevauche_iiiiih = KMH_deb_fin_traits_chevauches( iiiiih , 1 );
                    KMH_fin_trait_chevauche_iiiiih = KMH_deb_fin_traits_chevauches( iiiiih , 2 );
                    x_deb_trait_chevauche_iiiiih = points(min(KMH_deb_trait_chevauche_iiiiih,size(points,1)),1);
                    y_deb_trait_chevauche_iiiiih = points(min(KMH_deb_trait_chevauche_iiiiih,size(points,1)),2);
                    x_fin_trait_chevauche_iiiiih = points(min(KMH_fin_trait_chevauche_iiiiih,size(points,1)),1);
                    y_fin_trait_chevauche_iiiiih = points(min(KMH_fin_trait_chevauche_iiiiih,size(points,1)),2);
                    
                    delta_x_point_deb_point_fin_trait_courant = x_fin_trait_chevauche_iiiiih - points( KMH_deb_fin_traits_chevauches( indice_traits_deb , 1 ) ,  1);
                    delta_y_point_deb_point_fin_trait_courant = y_fin_trait_chevauche_iiiiih - points( KMH_deb_fin_traits_chevauches( indice_traits_deb , 1 ) ,  2);
                    delta_x_point_fin_point_deb_trait_courant = points( min(KMH_deb_fin_traits_chevauches( indice_traits_arrivee , 2 ),size(points,1)) ,  1)  -  x_deb_trait_chevauche_iiiiih;
                    delta_y_point_fin_point_deb_trait_courant = points( min(KMH_deb_fin_traits_chevauches( indice_traits_arrivee , 2 ),size(points,1)) ,  2)  -  y_deb_trait_chevauche_iiiiih;

                    Vect_x_y_points_limites_des_traits = [Vect_x_y_points_limites_des_traits; x_deb_trait_chevauche_iiiiih, y_deb_trait_chevauche_iiiiih, x_fin_trait_chevauche_iiiiih, y_fin_trait_chevauche_iiiiih];
                    
                    Vect_delta_x_point_deb_point_fin_trait_courant = [Vect_delta_x_point_deb_point_fin_trait_courant; delta_x_point_deb_point_fin_trait_courant];
                    Vect_delta_y_point_deb_point_fin_trait_courant = [Vect_delta_y_point_deb_point_fin_trait_courant; delta_y_point_deb_point_fin_trait_courant];
                    Vect_delta_x_point_fin_point_deb_trait_courant = [Vect_delta_x_point_fin_point_deb_trait_courant; delta_x_point_fin_point_deb_trait_courant];
                    Vect_delta_y_point_fin_point_deb_trait_courant = [Vect_delta_y_point_fin_point_deb_trait_courant; delta_y_point_fin_point_deb_trait_courant];

                    Vect_angle_courbure_variation = [Vect_angle_courbure_variation; vect_param_1( 8 , iiiiih ) , vect_param_1( 9 , iiiiih ) ];
                end
                KMH_fin_trait_mi_chemin = KMH_deb_fin_traits_chevauches( fix( (indice_traits_deb + indice_traits_arrivee) / 2 ) , 2 );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                figure (51);
                subplot(2,1,2);
                axis equal;
                plot( points(KMH_deb_troncon_trajectoire_a_enregistrer:KMH_fin_troncon_trajectoire_a_enregistrer,1) , points(KMH_deb_troncon_trajectoire_a_enregistrer:KMH_fin_troncon_trajectoire_a_enregistrer,2) , 'Color' , [.9 .0 .0] , 'lineStyle' , '-.');%'.b');
                hold on;
                
                %%% pause
                               
                nbr_traits_exacte = size(Vect_x_y_points_limites_des_traits , 1);               
                min_x_ensemble_segment_N_traits = ( min( [Vect_x_y_points_limites_des_traits( : , 1 ) ; Vect_x_y_points_limites_des_traits( : , 3 )] ) );
                min_y_ensemble_segment_N_traits = ( min( [Vect_x_y_points_limites_des_traits( : , 2 ) ; Vect_x_y_points_limites_des_traits( : , 4 )] ) );
                max_x_ensemble_segment_N_traits = ( max( [Vect_x_y_points_limites_des_traits( : , 1 ) ; Vect_x_y_points_limites_des_traits( : , 3 )] ) );
                max_y_ensemble_segment_N_traits = ( max( [Vect_x_y_points_limites_des_traits( : , 2 ) ; Vect_x_y_points_limites_des_traits( : , 4 )] ) );
                hauteur__normalisee_segment_N_traits = abs(   max_y_ensemble_segment_N_traits  -   min_y_ensemble_segment_N_traits   );
                largeur__normalisee_segment_N_traits = abs(   max_x_ensemble_segment_N_traits  -   min_x_ensemble_segment_N_traits   );
                
                rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm = ( ( Vect_x_y_points_limites_des_traits( 1 , 1 ) ) - ( min_x_ensemble_segment_N_traits ) ) / largeur__normalisee_segment_N_traits;
                rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm = ( ( Vect_x_y_points_limites_des_traits( 1 , 2 ) ) - ( min_y_ensemble_segment_N_traits ) ) / hauteur__normalisee_segment_N_traits;
                rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm = ( ( Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 3 ) ) - ( min_x_ensemble_segment_N_traits ) ) / largeur__normalisee_segment_N_traits;
                rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm = ( ( Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 4 ) ) - ( min_y_ensemble_segment_N_traits ) ) / hauteur__normalisee_segment_N_traits;
                rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm = ( ( Vect_x_y_points_limites_des_traits( max( fix( nbr_traits_exacte / 2 ) , 1 ) , 3 ) ) - ( min_x_ensemble_segment_N_traits ) ) / largeur__normalisee_segment_N_traits;
                rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm = ( ( Vect_x_y_points_limites_des_traits( max( fix( nbr_traits_exacte / 2 ) , 1 ) , 4 ) ) - ( min_y_ensemble_segment_N_traits ) ) / hauteur__normalisee_segment_N_traits;
                
                delta_x_deb_mi_chemin = Vect_delta_x_point_deb_point_fin_trait_courant( max( fix( nbr_traits_exacte / 2 ) , 1 ) );
                delta_y_deb_mi_chemin = Vect_delta_y_point_deb_point_fin_trait_courant( max( fix( nbr_traits_exacte / 2 ) , 1 ) );
                delta_x_mi_chemin_fin = Vect_delta_x_point_fin_point_deb_trait_courant( min( fix( nbr_traits_exacte / 2 ) + 1 , nbr_traits_exacte ) );
                delta_y_mi_chemin_fin = Vect_delta_y_point_fin_point_deb_trait_courant( min( fix( nbr_traits_exacte / 2 ) + 1 , nbr_traits_exacte ) );
                
                Teta_G_deb_des_N_traits = atan( tan( Vect_angle_courbure_variation( 1 , 1 ) ) );
                Teta_G_fin_des_N_traits = atan( tan( Vect_angle_courbure_variation( nbr_traits_exacte , 1 ) ) );
                indice_trait_mi_chemin = max( fix( nbr_traits_exacte / 2 ) , 1 );
                min_x_mi_chemin_segment = ( min( [Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 1 ) ; Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 3 )] ) );
                min_y_mi_chemin_segment = ( min( [Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 2 ) ; Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 4 )] ) );
                max_x_mi_chemin_segment = ( max( [Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 1 ) ; Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 3 )] ) );
                max_y_mi_chemin_segment = ( max( [Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 2 ) ; Vect_x_y_points_limites_des_traits( 1:indice_trait_mi_chemin , 4 )] ) );
                hauteur__normalisee_mi_chemin_segment = abs(   max_y_mi_chemin_segment  -   min_y_mi_chemin_segment   );
                largeur__normalisee_mi_chemin_segment = abs(   max_x_mi_chemin_segment  -   min_x_mi_chemin_segment   );
                rap_hauteur_mi_chemin_sur_hauteur_totale = hauteur__normalisee_mi_chemin_segment / ( max(hauteur__normalisee_segment_N_traits , 1.0e-7) );
                rap_largeur_mi_chemin_sur_largeur_totale = largeur__normalisee_mi_chemin_segment / ( max(largeur__normalisee_segment_N_traits , 1.0e-7) );
                rap_largeur_mi_chemin_sur_hauteur_totale = largeur__normalisee_mi_chemin_segment / ( max(hauteur__normalisee_segment_N_traits , 1.0e-7) );
                mean_abs_angle_courbure_variation_of_signicatif_strokes = 0;
                nbr_traits_exacte_a_considerer_dans_le_calcul_de_ce_parametre = 0;
                for iiiiih = 1 : nbr_traits_exacte
                    if ( ( abs(Vect_x_y_points_limites_des_traits( iiiiih , 1 ) - Vect_x_y_points_limites_des_traits( iiiiih , 3 )) / largeur__normalisee_segment_N_traits  >= 0.1 ) | ( abs(Vect_x_y_points_limites_des_traits( iiiiih , 2 ) - Vect_x_y_points_limites_des_traits( iiiiih , 4 )) / hauteur__normalisee_segment_N_traits  >= 0.1 ) )
                       mean_abs_angle_courbure_variation_of_signicatif_strokes = mean_abs_angle_courbure_variation_of_signicatif_strokes + abs( atan(tan(Vect_angle_courbure_variation(iiiiih , 1))) );
                       nbr_traits_exacte_a_considerer_dans_le_calcul_de_ce_parametre = nbr_traits_exacte_a_considerer_dans_le_calcul_de_ce_parametre + 1;
                    end
                end
                mean_abs_angle_courbure_variation_of_signicatif_strokes = mean_abs_angle_courbure_variation_of_signicatif_strokes / nbr_traits_exacte_a_considerer_dans_le_calcul_de_ce_parametre;
               
                if ( nbr_traits_exacte >= 2 )
                   if ( ( abs(Vect_x_y_points_limites_des_traits( 1 , 1 ) - Vect_x_y_points_limites_des_traits( 1 , 3 )) / largeur__normalisee_segment_N_traits  >= 0.1 ) | ( abs(Vect_x_y_points_limites_des_traits( 1 , 2 ) - Vect_x_y_points_limites_des_traits( 1 , 4 )) / hauteur__normalisee_segment_N_traits  >= 0.1 ) )
                      courbure_absolue_init_des_N_traits = abs( Vect_angle_courbure_variation(1 , 1) - Vect_angle_courbure_variation(1 , 2) );
                      courbure_relative_init_des_N_traits = ( Vect_angle_courbure_variation(1 , 1) - Vect_angle_courbure_variation(1 , 2) );
                   else
                      courbure_absolue_init_des_N_traits = 0;
                      courbure_relative_init_des_N_traits = 0;
                   end
                   if ( ( abs(Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 1 ) - Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 3 )) / largeur__normalisee_segment_N_traits  >= 0.1 ) | ( abs(Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 2 ) - Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 4 )) / hauteur__normalisee_segment_N_traits  >= 0.1 ) )
                      courbure_absolue_fin_des_N_traits = abs( Vect_angle_courbure_variation(nbr_traits_exacte , 2) - Vect_angle_courbure_variation(nbr_traits_exacte , 1) );
                      courbure_relative_fin_des_N_traits = ( Vect_angle_courbure_variation(nbr_traits_exacte , 2) - Vect_angle_courbure_variation(nbr_traits_exacte , 1) );
                   else
                      courbure_absolue_fin_des_N_traits = 0;
                      courbure_relative_fin_des_N_traits = 0;
                   end
                   courbure_absolue_des_N_traits = (courbure_absolue_init_des_N_traits + courbure_absolue_fin_des_N_traits) / 2;
                   courbure_relative_des_N_traits = (courbure_relative_init_des_N_traits + courbure_relative_fin_des_N_traits) / 2;
                else
                   courbure_absolue_des_N_traits = 0;
                   courbure_relative_des_N_traits = 0;
                end
                for iiiiih = 2 : nbr_traits_exacte
                   if ( ( abs(Vect_x_y_points_limites_des_traits( iiiiih , 1 ) - Vect_x_y_points_limites_des_traits( iiiiih , 3 )) / largeur__normalisee_segment_N_traits  >= 0.1 ) | ( abs(Vect_x_y_points_limites_des_traits( iiiiih , 2 ) - Vect_x_y_points_limites_des_traits( iiiiih , 4 )) / hauteur__normalisee_segment_N_traits  >= 0.1 ) )
                      courbure_absolue_des_N_traits = courbure_absolue_des_N_traits + ( abs( Vect_angle_courbure_variation(iiiiih , 1) - Vect_angle_courbure_variation(iiiiih -1 , 1) ) );
                      courbure_relative_des_N_traits = courbure_relative_des_N_traits + ( Vect_angle_courbure_variation(iiiiih , 1) - Vect_angle_courbure_variation(iiiiih -1 , 1) );
                   end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                uuh = 5;
                Teta_voisinage_avant_point_mi_chemin = zeros(5,1);
                Teta_voisinage_apres_point_mi_chemin = zeros(5,1);
                for uuh = 1 : 5
                    if (KMH_fin_trait_mi_chemin > size(points,1))
                       KMH_fin_trait_mi_chemin = size(points,1);
                    end
                      
                    denum = (points(max((KMH_fin_trait_mi_chemin - uuh),1) , 1) - points(KMH_fin_trait_mi_chemin , 1));
                    if denum == 0
                       denum= 1.0e-7;
                    end
                    cotion = (points(max((KMH_fin_trait_mi_chemin - uuh),1) , 2) - points(KMH_fin_trait_mi_chemin , 2)) / denum;
                    Teta_voisinage_avant_point_mi_chemin(uuh) = atan( cotion );

                    denum = (points(min((KMH_fin_trait_mi_chemin + uuh),size(points,1)) , 1) - points(KMH_fin_trait_mi_chemin , 1));
                    if denum == 0
                       denum= 1.0e-7;
                    end
                    cotion = (points(min((KMH_fin_trait_mi_chemin + uuh),size(points,1)) , 2) - points(KMH_fin_trait_mi_chemin , 2)) / denum;
                    Teta_voisinage_apres_point_mi_chemin(uuh) = atan( cotion );
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                size_Vect_ang_curv_var = size(Vect_angle_courbure_variation , 1);                  
                   le_segment_N_traits_est_il_classe = 0;

                   if (indice_traits_deb <= 2) & (indice_traits_arrivee < nbr_traits_methode_consideree)
                      pre_pre_classe = 'Segments in the beginning\';
                      %% choix_menu = menu('Choisir la classe du segment considéré','Hampe lam Debut................','Nabra Debut....................','Courbe Ouverte à Droite Debut.','Occlusion Debut...............','Large Occlusion Debut.........','7a Debut ..............');
                      %% choix_menu = choix_menu
                      
                      %% if (choix_menu == 1)
                         if ( (hauteur__normalisee_segment_N_traits >= 2) & ((largeur__normalisee_segment_N_traits/hauteur__normalisee_segment_N_traits)<=0.75) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.8) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.30) & ( (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.70)  | ((rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.80) & (courbure_absolue_des_N_traits >= 2*pi/4)) ) & (courbure_absolue_des_N_traits <= 3.5*pi/4)  )

                         pre_classe = 'beginning ascending shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 2)
                      elseif ( (hauteur__normalisee_segment_N_traits <= 40) & (min_y_ensemble_segment_N_traits >= 45) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.6) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.8) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.45 ) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <=0.45) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <=0.4) & (courbure_absolue_des_N_traits <= 1.5*pi) & (abs(Vect_angle_courbure_variation(1 , 1)) >= 1.2*pi/4 ) )
  
                      pre_classe = 'i in the beginning\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 3)
                       elseif ( (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.45) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.45) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.75) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.45) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.2) & (courbure_absolue_des_N_traits >= 3*pi/4) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.5)  )
                       pre_classe = 'opened right curve\';
                        
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)
                      elseif ( (hauteur__normalisee_segment_N_traits <= 70) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.7 ) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.4) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm > (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm + (1/15))) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.2) & (rap_largeur_mi_chemin_sur_hauteur_totale > 0.15) & (courbure_absolue_des_N_traits >= 3*pi/4) & (abs(courbure_relative_des_N_traits) >= 3*pi/4) & (abs(Vect_angle_courbure_variation(min(3,size_Vect_ang_curv_var) , 1) - Vect_angle_courbure_variation(2 , 1))<=1.4*pi/2) ) & (size_Vect_ang_curv_var >= 3)
                         pre_classe = 'broad occlusion in the beginning\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 5)
                      elseif ( (hauteur__normalisee_segment_N_traits <= 70) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.25) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.3) & (courbure_absolue_des_N_traits >= 3.5*pi/4) & (abs(courbure_relative_des_N_traits) >= 3.5*pi/4) )
                         pre_classe = 'occlusion e beginning\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 6)
                        elseif ( (hauteur__normalisee_segment_N_traits >= 1.50) & ((rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.45)|((rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.45)&(abs(Vect_angle_courbure_variation(1 , 1))>=0.7*pi/2))) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.65) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 0.75) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes >= (1.2*pi/4))  )
                         pre_classe = 'beginning descending  shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                        elseif ( (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.70) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.75) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.3) & (courbure_absolue_des_N_traits >= 3*pi/4) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.6)  )                         

                            pre_classe = 'opened left curve\';
                         le_segment_N_traits_est_il_classe = 1;
                      end
                      
                      if (le_segment_N_traits_est_il_classe == 0)
                           if ((hauteur__normalisee_segment_N_traits <= 0.90) & (min_y_ensemble_segment_N_traits >= 45))

                            pre_classe = 'residual segments in the median\';
                            le_segment_N_traits_est_il_classe = 1;
                           elseif ((hauteur__normalisee_segment_N_traits >= 0.90) & (min_y_ensemble_segment_N_traits >= 45))                              

                            pre_classe = 'residual segments in the top\';
                            le_segment_N_traits_est_il_classe = 1;
                         else
                            pre_classe = 'residual segments in the low\';
                            le_segment_N_traits_est_il_classe = 1;
                         end                          
                      end

                      
                   elseif (indice_traits_deb > 2) & (indice_traits_arrivee < nbr_traits_methode_consideree)
                      pre_pre_classe = 'Segments in the middle\';
                      %% choix_menu = choix_menu

                      %% if (choix_menu == 1)
                          if (  (hauteur__normalisee_segment_N_traits >= 1) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.4) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.4) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.70) & (courbure_absolue_des_N_traits <= 2.2*pi) & (courbure_absolue_des_N_traits >= 0.8*pi) & (abs(Teta_voisinage_apres_point_mi_chemin(2)) >= 1.35*pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(2)) >= 1.35*pi/4) )                             
                          
                          pre_classe = 'medium ascending shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 2)
                      elseif (  (hauteur__normalisee_segment_N_traits <= 40) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.50) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.35) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.95) & (courbure_absolue_des_N_traits <= 2.2*pi) & (courbure_absolue_des_N_traits >= 0.8*pi/2) & (((abs(Teta_voisinage_apres_point_mi_chemin(5)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(5)) >= pi/4))|((abs(Teta_voisinage_apres_point_mi_chemin(2)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(2)) >= pi/4))|((abs(Teta_voisinage_apres_point_mi_chemin(3)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(3)) >= pi/4))|((abs(Teta_voisinage_apres_point_mi_chemin(4)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(4)) >= pi/4))) & (sign(Teta_voisinage_apres_point_mi_chemin(3)) ~= sign(Teta_voisinage_avant_point_mi_chemin(3))) )
                         pre_classe = 'e in the medium\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 3)

                      elseif (  (hauteur__normalisee_segment_N_traits <= 45) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.60) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.60) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.15) & (courbure_absolue_des_N_traits <= 2.2*pi) & (courbure_absolue_des_N_traits >= 0.8*pi/2) & (((abs(Teta_voisinage_apres_point_mi_chemin(5)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(5)) >= pi/4))|((abs(Teta_voisinage_apres_point_mi_chemin(2)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(2)) >= pi/4))|((abs(Teta_voisinage_apres_point_mi_chemin(3)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(3)) >= pi/4))|((abs(Teta_voisinage_apres_point_mi_chemin(4)) >= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(4)) >= pi/4))) & (sign(Teta_voisinage_apres_point_mi_chemin(3)) ~= sign(Teta_voisinage_avant_point_mi_chemin(3))) )

                         pre_classe = 'n in the medium\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 3)
                      
                      elseif (  (hauteur__normalisee_segment_N_traits <= 60) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.2) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.2) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.90) & (courbure_absolue_des_N_traits <= 2.2*pi) & (courbure_absolue_des_N_traits >= 1.1*pi) & (abs(Teta_voisinage_apres_point_mi_chemin(2)) <= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(2)) <= pi/4) )
                         pre_classe = 'medium occlusion\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)
                      elseif (  (hauteur__normalisee_segment_N_traits <= 70) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.5) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.5) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.40) & (courbure_absolue_des_N_traits <= 2.2*pi) & (courbure_absolue_des_N_traits >= 0.8*pi) & (abs(Teta_voisinage_apres_point_mi_chemin(2)) <= 1.5*pi/4) & (abs(Teta_voisinage_apres_point_mi_chemin(2) - Teta_voisinage_avant_point_mi_chemin(2)) <= 0.4*pi/4) & (abs(Vect_angle_courbure_variation(min(3,size_Vect_ang_curv_var) , 1) - Vect_angle_courbure_variation(2 , 1)) <= pi) ) & (size_Vect_ang_curv_var >= 3)
                         pre_classe = 'broad occlusion in the medium\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 5)
                      elseif ( (hauteur__normalisee_segment_N_traits >= 1) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.7) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.30) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 0.65) & (courbure_absolue_des_N_traits <= 5.5*pi/4)  & (mean_abs_angle_courbure_variation_of_signicatif_strokes >= (1.2*pi/4))  )

                         pre_classe = 'medium descending  shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                      elseif ( (hauteur__normalisee_segment_N_traits <= 80) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.6 ) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.3) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm > (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm + (1/15))) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.3) & (rap_largeur_mi_chemin_sur_hauteur_totale > 0.15) & (courbure_absolue_des_N_traits >= 3*pi/4) & (abs(courbure_relative_des_N_traits) >= 3*pi/4) & (abs(Vect_angle_courbure_variation(min(3,size_Vect_ang_curv_var) , 1) - Vect_angle_courbure_variation(2 , 1))<=1.4*pi/2) ) & (size_Vect_ang_curv_var >= 3)

                         pre_classe = 'begin of occlusion in the medium\';
                         le_segment_N_traits_est_il_classe = 1;
                      elseif ( (hauteur__normalisee_segment_N_traits <= 40) & (min_y_ensemble_segment_N_traits >= 45) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.6) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.8) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.45 ) & (courbure_absolue_des_N_traits <= 1.5*pi) & (abs(Vect_angle_courbure_variation(1 , 1)) >= 1.2*pi/4 ) )
                         pre_classe = 'medium valley\';
                         le_segment_N_traits_est_il_classe = 1;
                         
                      elseif ( (hauteur__normalisee_segment_N_traits <= 70) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.6) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.3) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.75) & (courbure_absolue_des_N_traits >= 2.5*pi/4) & (courbure_absolue_des_N_traits >= 2.2*pi)  )
                         pre_classe = 'r in the medium\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)

                      end
                      
                      if (le_segment_N_traits_est_il_classe == 0)
                           if ((hauteur__normalisee_segment_N_traits <= 0.90) & (min_y_ensemble_segment_N_traits >= 50))

                            pre_classe = 'residual segments in the median\';
                            le_segment_N_traits_est_il_classe = 1;
                           elseif ((hauteur__normalisee_segment_N_traits >= 0.90) & (min_y_ensemble_segment_N_traits >= 50))

                            pre_classe = 'residual segments in the top\';
                            le_segment_N_traits_est_il_classe = 1;
                         else
                            pre_classe = 'residual segments in the low\';
                            le_segment_N_traits_est_il_classe = 1;
                         end                          
                      end
                      
                      
                      
                   elseif (indice_traits_deb > 1) & (indice_traits_arrivee == nbr_traits_methode_consideree)
                      pre_pre_classe = 'Segments in the end\';
                      %% choix_menu = choix_menu
                      
                      hbhbhb = 3;
                      
                      %% if (choix_menu == 1)
                        if ( (hauteur__normalisee_segment_N_traits >= 1) & (min_y_ensemble_segment_N_traits >= 34) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.52) & (hbhbhb == 3) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.35) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.85) & (courbure_absolue_des_N_traits <= 2*pi) & (courbure_absolue_des_N_traits >= 0.6*pi/2)  )
                         pre_classe = 'end ascending shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 2)
                      elseif ( (hauteur__normalisee_segment_N_traits <= 40) & (min_y_ensemble_segment_N_traits >= 47) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.5) & (hbhbhb == 3) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.25) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.85) & (courbure_absolue_des_N_traits <= 2*pi) & (courbure_absolue_des_N_traits >= 0.6*pi/2)  )
                         pre_classe = 'i in the end\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 3)                     
                      elseif ( (hauteur__normalisee_segment_N_traits <= 90) & (rap_delt_x_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.3) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.6) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.6) & (courbure_absolue_des_N_traits >= pi/4) & (courbure_absolue_des_N_traits >= pi)  )
                         pre_classe = 'end occlusion\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)
                      elseif (  (hauteur__normalisee_segment_N_traits <= 65) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.25) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.25) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.80) & (courbure_absolue_des_N_traits <= 2*pi) & (courbure_absolue_des_N_traits >= 1.1*pi) & (abs(Teta_voisinage_apres_point_mi_chemin(2)) <= pi/4) & (abs(Teta_voisinage_avant_point_mi_chemin(2)) <= pi/4) )
                         pre_classe = 'e in the end\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)
                        elseif ( (hauteur__normalisee_segment_N_traits <= 90) & (min_y_ensemble_segment_N_traits <= 75) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.3) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.65) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.5) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.45) & (courbure_absolue_des_N_traits <= 2.5*pi) & (courbure_absolue_des_N_traits >= 0.5*pi/2) )

                         pre_classe = 'curvy leg end\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 5)
                        elseif ( (hauteur__normalisee_segment_N_traits <= 95) & (min_y_ensemble_segment_N_traits <= 75) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.12) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.4) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.45) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.26) & (courbure_absolue_des_N_traits <= 2*pi) & (courbure_absolue_des_N_traits >= pi/2) )
                         pre_classe = 'pocket leg end\';
                         le_segment_N_traits_est_il_classe = 1;
                        elseif ( (hauteur__normalisee_segment_N_traits >= 1) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.5) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.5) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 0.8) & (courbure_absolue_des_N_traits <= 5.5*pi/4)  )

                         pre_classe = 'end descending shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                      end
                      
                      if (le_segment_N_traits_est_il_classe == 0)
                              if ((hauteur__normalisee_segment_N_traits <= 0.90) & (min_y_ensemble_segment_N_traits >= 50))
                            pre_classe = 'residual segments in the median\';
                            le_segment_N_traits_est_il_classe = 1;
                           elseif ((hauteur__normalisee_segment_N_traits >= 0.90) & (min_y_ensemble_segment_N_traits >= 50))
                            pre_classe = 'residual segments in the top\';
                            le_segment_N_traits_est_il_classe = 1;
                         else
                            pre_classe = 'residual segments in the low\';
                            le_segment_N_traits_est_il_classe = 1;
                         end                          
                      end                   
                      
                   elseif (indice_traits_deb == 1) & (indice_traits_arrivee == nbr_traits_methode_consideree)
                      pre_pre_classe = 'isolated segments\';
                      %% choix_menu = choix_menu
                      
                      %% if (choix_menu == 1)
                      if ( (hauteur__normalisee_segment_N_traits >= 1) & ((rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.40)|((rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.40)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 0.6) & (courbure_absolue_des_N_traits <= 5.5*pi/4) ))
                         pre_classe = 'isolated shaft\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 2)
                      elseif ( (hauteur__normalisee_segment_N_traits <= 92) & (min_y_ensemble_segment_N_traits <= 90) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.2) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.3) & (rap_delta_x_pnt_deb_lim_inf_rectangle_capable_sur_H_norm <= 0.70) & (rap_delta_x_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.3) & (courbure_absolue_des_N_traits <= 2.2*pi) & (courbure_absolue_des_N_traits >= pi/2) )

                         pre_classe = 'isolated pocket\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)
                        elseif ( (hauteur__normalisee_segment_N_traits <= 85) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm >= 0.4 ) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.1) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.6) & (courbure_absolue_des_N_traits >= 3*pi/4) & (abs(courbure_relative_des_N_traits) >= 3*pi/4) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.2) )
                         pre_classe = 'isolated begin of occlusion\';
                         le_segment_N_traits_est_il_classe = 1;
                        elseif ( (hauteur__normalisee_segment_N_traits >= 0.5) & (rap_delt_y_pnt_mi_chm_lim_inf_rectangl_capable_sur_H_norm <= 0.8 ) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.1) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm >= 0.1) & (courbure_absolue_des_N_traits >= pi/4) & (abs(courbure_relative_des_N_traits) >= pi/4) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.2) & (sign(atan(tan(Vect_x_y_points_limites_des_traits( 1 , 1 ))))~=sign(atan(tan(Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 1 ))))) &         ((abs(atan(tan(Vect_x_y_points_limites_des_traits( 1 , 1 )))) <= 0.65*pi/2)&(abs(atan(tan(Vect_x_y_points_limites_des_traits( nbr_traits_exacte , 1 )))) <= 0.65*pi/2))     )
                         pre_classe = 'isolated O occlusion\';
                         le_segment_N_traits_est_il_classe = 1;
                      %% elseif (choix_menu == 4)

                      
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 36) & (max_y_ensemble_segment_N_traits >= 95) & (min_y_ensemble_segment_N_traits <= 63) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.83) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.35) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.3) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 1.9) )
%%%%%%%%%%%%%%%thameur
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 3) & (max_y_ensemble_segment_N_traits >= 80) & (min_y_ensemble_segment_N_traits <= 50) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.63) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.45) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.3) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 1.9) )
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 30) & (max_y_ensemble_segment_N_traits >= 70) & (min_y_ensemble_segment_N_traits <= 70) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.63) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.45) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.3) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 1.9) )
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 10) & (max_y_ensemble_segment_N_traits >= 50) & (min_y_ensemble_segment_N_traits <= 75) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.63) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.45) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.3) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 1.9) )
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 1) & (max_y_ensemble_segment_N_traits >= 50) & (min_y_ensemble_segment_N_traits <= 75) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.63) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.45) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.3) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 1.9) )
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 1) & (max_y_ensemble_segment_N_traits >= 50) & (min_y_ensemble_segment_N_traits <= 85) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.43) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.65) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.3) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 1.9) )
%                       elseif ( (hauteur__normalisee_segment_N_traits >= 0.5) & (max_y_ensemble_segment_N_traits >= 50) & (min_y_ensemble_segment_N_traits <= 75) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.63) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.5) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (1.4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.2) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) <= 2.9) )
                      elseif ( (hauteur__normalisee_segment_N_traits >= 0.5) & (max_y_ensemble_segment_N_traits >= 50) & (min_y_ensemble_segment_N_traits <= 90) & (rap_delta_y_pnt_deb_lim_inf_rectangle_capable_sur_H_norm >= 0.43) & (rap_delta_y_pnt_fin_lim_inf_rectangle_capable_sur_H_norm <= 0.65) & (courbure_absolue_des_N_traits <= 5.5*pi/4) & (mean_abs_angle_courbure_variation_of_signicatif_strokes <= (4*pi/4)) & ((largeur__normalisee_segment_N_traits / hauteur__normalisee_segment_N_traits) >= 0.15) )

                         pre_classe = 'isolated j\';
                         le_segment_N_traits_est_il_classe = 1;                   
                      end
                      
                      if (le_segment_N_traits_est_il_classe == 0)
                        if ((hauteur__normalisee_segment_N_traits <= 0.5) & (min_y_ensemble_segment_N_traits >= 50))

                            pre_classe = 'residual segments in the median\';
                            le_segment_N_traits_est_il_classe = 1;
                           
                         elseif ((hauteur__normalisee_segment_N_traits >= 0.5) & (min_y_ensemble_segment_N_traits >= 50))

                            pre_classe = 'residual segments in the top\';
                            le_segment_N_traits_est_il_classe = 1;
                         else
                            pre_classe = 'residual segments in the low\';
                            le_segment_N_traits_est_il_classe = 1;
                         end                          
                      end

                      
                   end
                   
                if (le_segment_N_traits_est_il_classe == 1)
                      
                   pre_pre_classe;
                   pre_classe;             
                   Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe = [chemin_acces_dossier_Base_Param_I_formes_visuelles, pre_pre_classe];
                   existence_repertoire = exist( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe );
                   end
                   Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe_FEPC = [chemin_acces_dossier_Base_Param_I_formes_visuelles_FEPC, pre_pre_classe];
                   existence_repertoire = exist( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe_FEPC , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe_FEPC );
                   end
                   Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe = [chemin_acces_dossier_Base_Param_II_formes_visuelles, pre_pre_classe];
                   existence_repertoire = exist( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe );
                   end
                   
                   Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe_FEPC = [chemin_acces_dossier_Base_Param_II_formes_visuelles_FEPC, pre_pre_classe];
                   existence_repertoire = exist( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe_FEPC , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe_FEPC );
                   end
                   
                   Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe = [Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe, pre_classe];
                   Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe_fig = [Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe, pre_classe,'fig\'];
                   existence_repertoire = exist( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe );
                   end
                    Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe_FEPC = [Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_pre_classe_FEPC, pre_classe];
                   existence_repertoire = exist( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe_FEPC , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe_FEPC );
                   end
                  Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe = [Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe, pre_classe];
                   Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe_fig = [Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe, pre_classe,'fig\'];
                   existence_repertoire = exist( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe );
                   end
                   Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe_FEPC = [Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_pre_classe_FEPC, pre_classe];
                   existence_repertoire = exist( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe_FEPC , 'dir' );
                   if (existence_repertoire == 0)
                      mkdir( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe_FEPC );
                   end
                   liste_files_I = dir( fullfile( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe , '*.txt') );
                   nombre_des_fichiers_echantillons_type_I = size(liste_files_I , 1);
                   
                   liste_files_I_FEPC = dir( fullfile( Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe_FEPC , '*.txt') );
                   nombre_des_fichiers_echantillons_type_I_FEPC = size(liste_files_I_FEPC , 1);
                   liste_files_II = dir( fullfile( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe , '*.txt') );
                   nombre_des_fichiers_echantillons_type_II = size(liste_files_II , 1);
                   liste_files_II_FEPC = dir( fullfile( Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe_FEPC , '*.txt') );
                   nombre_des_fichiers_echantillons_type_II_FEPC = size(liste_files_II_FEPC , 1);
                   numero_echantillon_type_I_a_enregistrer = nombre_des_fichiers_echantillons_type_I + 1;
                   numero_echantillon_type_I_a_enregistrer_en_caracter = num2str( numero_echantillon_type_I_a_enregistrer );
                   numero_echantillon_type_I_a_enregistrer_FEPC = nombre_des_fichiers_echantillons_type_I_FEPC + 1;
                   numero_echantillon_type_I_a_enregistrer_en_caracter_FEPC = num2str( numero_echantillon_type_I_a_enregistrer_FEPC );
                   numero_echantillon_type_II_a_enregistrer = nombre_des_fichiers_echantillons_type_II + 1;
                   numero_echantillon_type_II_a_enregistrer_en_caracter = num2str( numero_echantillon_type_II_a_enregistrer );
                   numero_echantillon_type_II_a_enregistrer_FEPC = nombre_des_fichiers_echantillons_type_II_FEPC + 1;
                   numero_echantillon_type_II_a_enregistrer_en_caracter_FEPC = num2str( numero_echantillon_type_II_a_enregistrer_FEPC );
                   pre_classe_sans_slash = pre_classe(1 : end - 1);
                   nom_du_fichier_Type_I_d_enregistrement_sans_extension = [ pre_classe_sans_slash , '_AOBS_' , numero_echantillon_type_I_a_enregistrer_en_caracter ];
                   nom_du_fichier_Type_I_d_enregistrement_sans_extension_FEPC = [ pre_classe_sans_slash , '_AOBSFEPC' , numero_echantillon_type_I_a_enregistrer_en_caracter_FEPC ];
                   nom_du_fichier_Type_II_d_enregistrement_sans_extension = [ pre_classe_sans_slash , '_SBS_' , numero_echantillon_type_II_a_enregistrer_en_caracter ];
                   nom_du_fichier_Type_II_d_enregistrement_sans_extension_FEPC = [ pre_classe_sans_slash , '_SBSFEPC' , numero_echantillon_type_II_a_enregistrer_en_caracter_FEPC ];
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   vect_param_1_p = vect_param_1( : , indice_traits_deb : indice_traits_arrivee );
                   for iiiiih = 1 : ( Nombre_de_traits_chevauches_d_une_forme_visuelle - (indice_traits_arrivee - indice_traits_deb + 1) )
                       vect_param_1_p = [ vect_param_1_p , zeros( size( vect_param_1 , 1) , 1 ) ];
                   end
                   vect_vect_param_1 = vect_param_1_p(:);
                   
                   vect_param_2_p = vect_param_2( : , indice_trait_deb_methode_simplifiee : indice_traits_arrivee_methode_simplifiee );
                   for iiiiih = 1 : ( Nombre_de_traits_double_d_une_forme_visuelle  - (indice_traits_arrivee_methode_simplifiee - indice_trait_deb_methode_simplifiee + 1) )
                       vect_param_2_p = [ vect_param_2_p , zeros( size( vect_param_2 , 1) , 1 ) ];
                   end
                   vect_vect_param_2 = vect_param_2_p(:);
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                   
                   size_vect_vect_param_1 = size(vect_vect_param_1)
                   size_vect_vect_param_2 = size(vect_vect_param_2)
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
                   %%%%%%%********%%%%%%%   Enregistrement des données paramétriques extraites   %%%%%%%********%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   %  pause
                   csvwrite([Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe  , nom_du_fichier_Type_I_d_enregistrement_sans_extension , '.txt'], vect_vect_param_1 );
                   x=vect_vect_param_1([8,18,28,38],:);
                   x=x(:);
                   vect_vect_param_1_FEPC = ExtractCPEF1(x,vect_vect_param_1)
                   csvwrite([Chemin_acces_Base_I_Forme_Visuelle_dossier_pre_classe_FEPC  , nom_du_fichier_Type_I_d_enregistrement_sans_extension_FEPC , '.txt'], vect_vect_param_1_FEPC );
                   csvwrite([Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe , nom_du_fichier_Type_II_d_enregistrement_sans_extension , '.txt'], vect_vect_param_2 );
      
                   x=vect_vect_param_2([12,26],:);
                   x=x(:);
                   vect_vect_param_2_FEPC = ExtractCPEF2(x,vect_vect_param_2)
                   csvwrite([Chemin_acces_Base_II_Forme_Visuelle_dossier_pre_classe_FEPC  , nom_du_fichier_Type_II_d_enregistrement_sans_extension_FEPC , '.txt'], vect_vect_param_2_FEPC );


                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                 
                   
                   
                end

                   %%% pause
                   
%                    pause(0.05);
                   
                   figure (51);
                   subplot(2,1,2);
                   axis equal;
                   plot(points(KMH_deb_troncon_trajectoire_a_enregistrer:KMH_fin_troncon_trajectoire_a_enregistrer,1), points(KMH_deb_troncon_trajectoire_a_enregistrer:KMH_fin_troncon_trajectoire_a_enregistrer,2), 'Color',[.3 .7 .999], 'lineStyle', '-.');%'.b');
                   hold on;
                   
                   
                    %%%%%%%%%%%%%%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%*************%%%%%%%%%   Enregistrement des figures   %%%%%%%%%*************%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
                   
                   indice_traits_deb = indice_traits_deb + 1;
                   indice_traits_arrivee_preced = indice_traits_arrivee;

            end
            
    
         else
            vect_param_1 = []; nbr_trait = 0; Vect_Erreurs = [];
            vect_param_2 = []; 
         end
      else
         vect_param_1 = []; nbr_trait = 0; Vect_Erreurs = [];
         vect_param_2 = []; 
      end
      
%       pause
      
  end
  
 
else
 nbr_pseudo_mot = 1;
 nbr_pseudo_mot_traite = 1;
 vect_param_1 = ones(10,1); %[1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
 mmatrice_param_1 = vect_param_1;
 
 vect_param_2 = ones(14,1); %[1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
 mmatrice_param_2 = vect_param_2;
 
 figure(47);
%  pause;

end

if size(mmatrice_param_1,1) == 0
 nbr_pseudo_mot = 1;
 nbr_pseudo_mot_traite = 1;
 vect_param_1 = ones(10,1); %[1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
 mmatrice_param_1 = vect_param_1;
 
 vect_param_2 = ones(14,1); %[1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1; 1];
 mmatrice_param_2 = vect_param_2;
 
 figure(48);
%  pause;
end
   


pause(0.5);

end
