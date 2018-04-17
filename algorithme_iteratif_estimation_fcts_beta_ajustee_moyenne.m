function [fonctions_beta_ajustees_moyennes_m] = algorithme_iteratif_estimation_fcts_beta_ajustee_moyenne( fonctions_beta_ajustee_a_gauche , fonctions_beta_ajustee_a_droite );

fonctions_beta_ajustees_moyennes_m = [];

nbr_beta_fctn = size(fonctions_beta_ajustee_a_gauche , 1);
for indice_beta_fctn = 1 : nbr_beta_fctn
    
    t0 = ( fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 1 ) + fonctions_beta_ajustee_a_droite( indice_beta_fctn , 1 ) ) / 2;
    tc = ( fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 2 ) + fonctions_beta_ajustee_a_droite( indice_beta_fctn , 2 ) ) / 2;
    t1 = ( fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 3 ) + fonctions_beta_ajustee_a_droite( indice_beta_fctn , 3 ) ) / 2;

    p_ajust_gauche = fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 4 );
    q_ajust_gauche = fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 5 );
    
    p_ajust_droite = fonctions_beta_ajustee_a_droite( indice_beta_fctn , 4 );
    q_ajust_droite = fonctions_beta_ajustee_a_droite( indice_beta_fctn , 5 );
    
    %%%%% Initialisation
    p_ajust_droite_valeur_iteration_preced = p_ajust_droite;
    q_ajust_droite_valeur_iteration_preced = q_ajust_droite;
    p_ajust_gauche_valeur_iteration_preced = p_ajust_gauche;
    q_ajust_gauche_valeur_iteration_preced = q_ajust_gauche;
    
    p_valeur_iteration_preced = p_ajust_droite_valeur_iteration_preced;
    q_valeur_iteration_preced = q_ajust_gauche_valeur_iteration_preced;
    
    p_valeur_iteration_courante = (p_ajust_gauche_valeur_iteration_preced + p_ajust_droite_valeur_iteration_preced) / 2;
    %q_valeur_iteration_courante = q_ajust_gauche_valeur_iteration_preced;
    q_valeur_iteration_courante = (p_valeur_iteration_courante * (((t1 - tc)/(tc - t0))));
    
    i = 1;
    delta_p = 4; delta_q = 4;
    
    while (i <= 50) & ( ( delta_p > 0.05 ) | ( delta_q > 0.05 ) )
          p_valeur_iteration_courante = (p_valeur_iteration_courante + p_ajust_droite_valeur_iteration_preced) / 2;
          q_valeur_iteration_courante = (p_valeur_iteration_courante * (((t1 - tc)/(tc - t0))));

          p_ajust_droite_valeur_iteration_courante = p_valeur_iteration_courante;
          %q_ajust_droite_valeur_iteration_courante = q_valeur_iteration_courante;

          q_valeur_iteration_courante = (q_valeur_iteration_courante + q_ajust_gauche_valeur_iteration_preced) / 2;
          p_valeur_iteration_courante = (q_valeur_iteration_courante * (((tc - t0)/(t1 - tc))));
          
          q_ajust_gauche_valeur_iteration_courante = q_valeur_iteration_courante;
          %p_ajust_gauche_valeur_iteration_courante = p_valeur_iteration_courante;
          
          p_valeur_iteration_courante = (p_valeur_iteration_courante + p_ajust_gauche_valeur_iteration_preced) / 2;
          q_valeur_iteration_courante = (p_valeur_iteration_courante * (((t1 - tc)/(tc - t0))));
          
          p_ajust_gauche_valeur_iteration_courante = p_valeur_iteration_courante;
          %q_ajust_gauche_valeur_iteration_courante = q_valeur_iteration_courante;
          
          q_valeur_iteration_courante = (q_valeur_iteration_courante + q_ajust_gauche_valeur_iteration_preced) / 2;
          p_valeur_iteration_courante = (q_valeur_iteration_courante * (((tc - t0)/(t1 - tc))));

          q_ajust_droite_valeur_iteration_courante = q_valeur_iteration_courante;
          %p_ajust_droite_valeur_iteration_courante = p_valeur_iteration_courante;
          
          delta_p = abs( p_valeur_iteration_courante - p_valeur_iteration_preced );
          delta_q = abs( q_valeur_iteration_courante - q_valeur_iteration_preced );
          
          
          p_valeur_iteration_preced = p_valeur_iteration_courante;
          q_valeur_iteration_preced = q_valeur_iteration_courante;

          p_ajust_droite_valeur_iteration_preced = p_ajust_droite_valeur_iteration_courante;
          q_ajust_droite_valeur_iteration_preced = q_ajust_droite_valeur_iteration_courante;
          p_ajust_gauche_valeur_iteration_preced = p_ajust_gauche_valeur_iteration_courante;
          q_ajust_gauche_valeur_iteration_preced = q_ajust_gauche_valeur_iteration_courante;

          i = i + 1;
          
    end
    
    %% i_houcine = i
    %% delta_p = delta_p
    %% delta_q = delta_q
    %% pause;
    
    K_hauteur = ( fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 6 ) + fonctions_beta_ajustee_a_droite( indice_beta_fctn , 6 ) ) / 2;
    rong_extremum = ( fonctions_beta_ajustee_a_gauche( indice_beta_fctn , 7 ) + fonctions_beta_ajustee_a_droite( indice_beta_fctn , 7 ) ) / 2;
    fonction_beta_ajustee_moyenne_resultante = [t0, tc, t1, p_valeur_iteration_courante, q_valeur_iteration_courante, K_hauteur, rong_extremum];
    
    fonctions_beta_ajustees_moyennes_m = [fonctions_beta_ajustees_moyennes_m; fonction_beta_ajustee_moyenne_resultante];
end






%%%       fonctions_beta_ajustee_a_gauche = [fonctions_beta_ajustee_a_gauche; t0, tc, t1, p, q, K_hauteur_actuel, rong_extremum ];


%%%       q = (p * (((t1 - tc)/(tc - t0))));
%%%       p = (q * (((tc - t0)/(t1 - tc))));



