function [vecteur_param, matrice_param, param_trajectoire_BC, DAT_BC, indiceTemps_deb, points, nbr_trait, e, Vect_Erreurs_BC, famille_tg_horiz] = new_beta_ellipse_stretegy_Beta_chevauchees_s(data_k, elimin_trait_excedent, rayon_filtre_V, sigma_p_filtre_V)
%function [matrice_param, nbr_trait,e,Vect_Erreurs_BC] = new_beta_ellipse_stretegy_Beta_chevauchees(data_k, elimin_trait_excedent, rayon_filtre_V, sigma_p_filtre_V)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction new_beta_ellipse_stretegy_Beta_chevauchees permet de décomposer le signal        %%%%%
%%%%% enligne en segements séparés par une levée suivie d'une posée du                             %%%%%
%%%%% stylo. Chaque segment sera ensuite traité dans la phase de                                   %%%%%
%%%%% modélisation dans ses aspects : vitesse (dynamique)et trajectoire (statique)                 %%%%%
%%%%% Cette fonction permet d'adapter le niveau de sensibilité de détection des                    %%%%%
%%%%% points de double inflexion de la vitesse suivant la précision de la                          %%%%%
%%%%% reconstrction elliptique. Elle permet aussi suivant le choix de l'utilisateur                %%%%%
%%%%% d'eliminer les traits excédentaires.                                                         %%%%%
%%%%% Entrées : signal en ligne, choix de l'utilisateur                                            %%%%%
%%%%% Sorties : matrice de paramètres de modélisation, nombre de traits, Erreurs de reconstruction %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Vect_Erreurs_BC = [ 0 0 0 0; 0 0 0 0; 0 0 0 0];
M = size(data_k,1);

%OF-line -> ON-line% [data_norm] = normalisation_H(data_k,M);
data_norm = data_k;

tf = ((size(data_norm,1))-1)/100;
pas = 0.01;
t = 0 : pas : tf ; 
T = t';
data_norm = [data_norm T];

[N_seg,data_mot] = copie_segmentation_tablet_2(data_norm);

  
for j = 1 : N_seg

    nbre_beta = [];
    X = [];
    Y = [];
    T = [];
    ii = data_mot(:,j);                % ce vecteur marque par ii(i) ~= 0 les points d'indice i appartenant au segment j
                                       % et par ii(i) = 0 ceux qui ne lui appartiennent pas
    XXX = data_mot(:,j+N_seg);
    YYY = data_mot(:,j+2*N_seg);
    PPP = data_mot(:,j+3*N_seg);
    TTT = data_mot(:,j+4*N_seg);

%%%%%%%%%%%%%%%%%%%% elimination des zeros :reconstruction de X,Y et T
%%%%%%%%%%%%%%%%%%%% les points d'indice ii = 0 à la colonne j des vecteurs X_seg et Y_seg
%%%%%%%%%%%%%%%%%%%% n'appartiennent pas au segment j
   for i = 1:M
       if ii(i) == 0          
           X = X;
           Y = Y;
           T = T;        
           
       else
           X = [X;XXX(i)];
           Y = [Y;YYY(i)];
           T = [T;TTT(i)]; 
       end 
      
   end
   hh = size(T); 

%    figure (1);
%    plot(X,Y,'g-');
%    hold on;
  


   [points,V,DPV,T,t, X_init, Y_init, T_init, pas_t] = copie_pre_traitement_seg_H_p(X, Y, T, rayon_filtre_V, sigma_p_filtre_V);
 

   Rap_lim1 = 0.285;%0.5;%0.7;             % seuil de detection grossiere des points d'inflexion en montée 
   Rap_lim2 = 0.2;%0.45;%0.7;              % seuil de detection grossiere des points d'inflexion en descente
   Rap_lim3 = 1;%0.8;%0.5;                 % seuil de detection fine des points d'inflexion en montée 
   Rap_lim4 = 0.7;%0.6;%0.45;              % seuil de detection fine des points d'inflexion en descente 
   RAP_lim = [Rap_lim1; Rap_lim2; Rap_lim3; Rap_lim4];

   Sensibilite_plus = [];
   Sensibilite_plus_BC = [];

   [fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1_i] = copie_extract_betas_chevauchees_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus_BC, pas_t);
    
   num_fig = 4;
   [param_trajectoire_BC, param_grandeur_BC, e_indice_BC, Vect_Erreurs_BC, Sensibilite_plus_BC, DAT_BC] = construction_par_arc_ellipse_adapte_beta_chevauchees(fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1_i, points, T, num_fig);

    
   Sensibilite_plus_memoire_BC = Sensibilite_plus_BC;
   if size(Sensibilite_plus_BC, 1) > 0   %requete d'augmenter la sensibilite de detection des points de double inflexion
      [fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1_i] = copie_extract_betas_chevauchees_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus_BC, pas_t);
      num_fig = 32;
      [param_trajectoire_BC, param_grandeur_BC, e_indice_BC, Vect_Erreurs_BC, Sensibilite_plus_BC, DAT_BC] = construction_par_arc_ellipse_adapte_beta_chevauchees(fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1_i, points, T, num_fig);
   end

   [reponse, X_rest, Y_rest, T_rest] = detection_segments_excedentaires_beta_chevauchees(param_grandeur_BC, tableau_affectation_intervalles, points, X_init, Y_init, T_init);

   nbr_repassage = 0;
   while (reponse == 0)&( elimin_trait_excedent == 'oui' )&(nbr_repassage < 4)%3)%4)

%       figure (11);
%       hold on;
%       plot(X_rest,Y_rest,'g-');
       
      Sensibilite_plus_BC = Sensibilite_plus_memoire_BC;
      [points,V,DPV,T,t, X_init, Y_init, T_init, pas_t] = copie_pre_traitement_seg_H_Re_p(X_rest,Y_rest,T_rest, rayon_filtre_V, sigma_p_filtre_V);
      [fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1_i] = copie_extract_betas_chevauchees_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus_BC, pas_t);
      num_fig = 31;
      [param_trajectoire_BC, param_grandeur_BC, e_indice_BC, Vect_Erreurs_BC, Sensibilite_plus_BC, DAT_BC] = construction_par_arc_ellipse_adapte_beta_chevauchees(fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1_i, points, T, num_fig);       
      [reponse, X_rest, Y_rest, T_rest] = detection_segments_excedentaires_beta_chevauchees(param_grandeur_BC, tableau_affectation_intervalles, points, X_init, Y_init, T_init);
      nbr_repassage = nbr_repassage + 1;
   end
   
   indiceTemps_deb = tableau_affectation_intervalles(1,7);          %= 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   famille_tg_horiz = [];
%    [famille_tg_horiz] = Detection_tg_horiz(DAT_BC, points, tableau_affectation_intervalles);
%    [familles_ligne_de_base] = Detection_ligne_de_base(DAT_BC, points, tableau_affectation_intervalles);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [modele_sortie] = collecte_donnees_Beta_chevauchees_etendue(fonctions_beta, param_trajectoire_BC);

    matrice_param = modele_sortie;
    vecteur_param = modele_sortie;
    nbr_trait = size(matrice_param, 1);
    e = e_indice_BC;

    mat = vecteur_param;
    mat1 = mat';
    %%% mat2 = mat1(:);
    %%% vecteur_param = mat2;
    vecteur_param = mat1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end

        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%




