function [vecteur_param, matrice_param, param_trajectoire_BC, DAT_BC, indiceTemps_deb, points, nbr_trait, e, Vect_Erreurs_BC, KMH_deb_fin_traits] = new_beta_ellipt_stretegy_Beta_compsnt_entrnmnt_preclass(data_k, elimin_trait_excedent, rayon_filtre_V, sigma_p_filtre_V)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction new_beta_ellipse_stretegy_Beta_composante_entrainement permet décomposer le signal %%%%%
%%%%% enligne en segements séparés par une levée suivie d'une posée du                               %%%%%
%%%%% stylo. Chaque segment sera ensuite traité dans la phase de                                     %%%%%
%%%%% modélisation dans ses aspects : vitesse (dynamique)et trajectoire (statique)                   %%%%%
%%%%% Cette fonction permet d'adapter le niveau de sensibilité de détection des                      %%%%%
%%%%% points de double inflexion de la vitesse suivant la précision de la                            %%%%%
%%%%% reconstrction elliptique. Elle permet aussi suivant le choix de l'utilisateur                  %%%%%
%%%%% d'eliminer les traits excédentaires.                                                           %%%%%
%%%%% Entrées : signal en ligne, choix de l'utilisateur                                              %%%%%
%%%%% Sorties : matrice de paramètres de modélisation, nombre de traits, Erreurs de reconstruction   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Vect_Erreurs_BC = [ 0 0 0 0; 0 0 0 0; 0 0 0 0];

M = size(data_k,1);

%OF-line -> ON-line% [data_norm] = normalisation_H(data_k,M);
data_norm = data_k;

tf= ((size(data_norm,1))-1)/100;
pas = 0.01;
t = 0 : pas : tf ; 
T = t';
data_norm = [data_norm T];

%[N_seg,data_mot] = copie_segmentation_tablet(data_norm);
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
    for i = 1 : M
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

    %figure (1);
    %hold on;
    %plot(X,Y,'g-');
 

%%%%%%%%% traitement de la liste X Y  du segment N° j   %%%%%%%%%%
 
    [points,V,DPV,T,t, X_init, Y_init, T_init, pas_t] = copie_pre_traitement_seg_H_p(X, Y, T, rayon_filtre_V, sigma_p_filtre_V);

    Rap_lim1 = 0.285;%0.5;%0.7;             % seuil de detection grossiere des points d'inflexion en montée 
    Rap_lim2 = 0.2;%0.45;%0.7;              % seuil de detection grossiere des points d'inflexion en descente
    Rap_lim3 = 1;%0.8;%0.5;                 % seuil de detection fine des points d'inflexion en montée 
    Rap_lim4 = 0.7;%0.6;%0.45;              % seuil de detection fine des points d'inflexion en descente 
    RAP_lim = [Rap_lim1; Rap_lim2; Rap_lim3; Rap_lim4];

    Sensibilite_plus = [];
    Sensibilite_plus_BC = [];

    [IMPULSION_BETA, ENTRAINEMENT, vect_X0XCX1] = copie_extract_beta_mot_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus, pas_t);
   
    num_fig = 4;
    [param_trajectoire_BC, param_grandeur, e_indice, Vect_Erreurs_BC, Sensibilite_plus, DAT_BC] = construction_par_arc_ellipse_adapte(IMPULSION_BETA,points,T, num_fig);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Sensibilite_plus_memoire = Sensibilite_plus;
    if size(Sensibilite_plus, 1) > 0   %requete d'augmenter la sensibilite de detection des points de double inflexion
       [IMPULSION_BETA, ENTRAINEMENT, vect_X0XCX1] = copie_extract_beta_mot_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus, pas_t);
       num_fig = 32;
       [param_trajectoire_BC, param_grandeur, e_indice, Vect_Erreurs_BC, Sensibilite_plus, DAT_BC] = construction_par_arc_ellipse_adapte(IMPULSION_BETA,points,T, num_fig);
    end

    [reponse, X_rest, Y_rest, T_rest] = detection_segments_excedentaires(param_grandeur, IMPULSION_BETA, points, X_init, Y_init, T_init);

    nbr_repassage = 0;
    while (reponse == 0)&( elimin_trait_excedent == 'oui' )&(nbr_repassage < 4)%3)%4)

       % figure (11);
       % hold on;
       % plot(X_rest,Y_rest,'g-');

       Sensibilite_plus = Sensibilite_plus_memoire;
       [points,V,DPV,T,t, X_init, Y_init, T_init, pas_t] = copie_pre_traitement_seg_H_Re_p(X_rest, Y_rest, T_rest, rayon_filtre_V, sigma_p_filtre_V);
       %[points,V,DPV,T,t, X_init, Y_init, T_init] = copie_pre_traitement_seg_H_filt_lineair(X,Y,T);
       [IMPULSION_BETA, ENTRAINEMENT, vect_X0XCX1] = copie_extract_beta_mot_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus, pas_t);
       num_fig = 31;
       [param_trajectoire_BC, param_grandeur, e_indice, Vect_Erreurs_BC, Sensibilite_plus, DAT_BC] = construction_par_arc_ellipse_adapte(IMPULSION_BETA,points,T, num_fig);
       
       [reponse, X_rest, Y_rest, T_rest] = detection_segments_excedentaires(param_grandeur, IMPULSION_BETA, points, X_init, Y_init, T_init);
       nbr_repassage = nbr_repassage + 1;
    end


    indiceTemps_deb = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [modele_sortie, KMH_deb_fin_traits] = collecte_donnees_etendue_preclass(IMPULSION_BETA, ENTRAINEMENT, param_trajectoire_BC);

    matrice_param = modele_sortie;
    vecteur_param = modele_sortie;
    nbr_trait = size(matrice_param, 1);
    e = e_indice;

    mat = matrice_param;
    mat1 = mat';
    %% mat2 = mat1(:);
    %% vecteur_param = mat2;
    vecteur_param = mat1;
    
end


        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%


