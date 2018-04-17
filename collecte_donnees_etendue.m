function [modele_sortie] = collecte_donnees_etendue(IMPULSION_BETA, ENTRAINEMENT, param_trajectoire)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction collecte_donnees permet de collecter les paramètre des  %%%%%
%%%%% impulsions Bêta, de la composante d'entrainement et des arcs        %%%%%
%%%%% d'ellipse modélisant le signal en ligne dans une seule matrice qui  %%%%%
%%%%% représente le modèle de sortie                                      %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nbr_seg = size(IMPULSION_BETA, 1);

T0 = IMPULSION_BETA(:,6);
Tc = IMPULSION_BETA(:,1);
T1 = IMPULSION_BETA(:,7);
P = IMPULSION_BETA(:,3);
K = IMPULSION_BETA(:,5);
Vini = ENTRAINEMENT(:,3);
Vfin = ENTRAINEMENT(:,4);

DELTA_T = [];
RAP_Tc = [];
P0SITION_TRAIT = [];
Rapport_Amplitude_Impulsions_Entrainement = [];

for  i = 1 : nbr_seg
     delta_T = T1(i) - T0(i);
     DELTA_T = [DELTA_T; delta_T];
     rap_Tc = (Tc(i) - T0(i)) / delta_T;
     RAP_Tc = [RAP_Tc; rap_Tc];
     
     position_trait = 1;
     if (i == 1)&(nbr_seg >= 2)
        position_trait = 1;                 %%%%% premier trait au début de la trajectoire
     elseif (i == 2)&(nbr_seg >= 3)
        position_trait = 2;                 %%%%% 2 ème trait au début de la trajectoire
     elseif (i > 2) & (i < (nbr_seg-1)) & (nbr_seg >= 5)
        position_trait = 3;                 %%%%% trait au milieu de la trajectoire
     elseif (i > 2) & (i == (nbr_seg-1))
        position_trait = 4;                 %%%%% avant dernier trait à la fin de la trajectoire
     elseif (i > 1) & (i == nbr_seg)
        position_trait = 5;                 %%%%% dernier trait à la fin de la trajectoire
     elseif (i == 1) & (i == nbr_seg)
        position_trait = 0.1;               %%%%% trait isolé
     end
     P0SITION_TRAIT = [P0SITION_TRAIT; position_trait];

     K_i = IMPULSION_BETA(i,5);
     Vini_i = ENTRAINEMENT(i,3);
     Vfin_i = ENTRAINEMENT(i,4);
     
     V_entrainement_moyenne = ( Vini_i + Vfin_i ) / 2;
     
     if ( V_entrainement_moyenne == 0 )
        V_entrainement_moyenne = 1.0e-10;
     end

     %% rapport_amplitude_impulsion_vitesse_entrainement_minimale = K_i / ( V_entrainement_moyenne );
     %% rapport_amplitude_impulsion_vitesse_entrainement_minimale = (K_i + 1) / (V_entrainement_moyenne + 1);
     rapport_amplitude_impulsion_vitesse_entrainement_minimale = (K_i + 0.1) / (V_entrainement_moyenne + 0.1);
     Rapport_Amplitude_Impulsions_Entrainement = [Rapport_Amplitude_Impulsions_Entrainement; rapport_amplitude_impulsion_vitesse_entrainement_minimale];

end

RAP_Ki_SUR_Entrainement = Rapport_Amplitude_Impulsions_Entrainement;

% % % LAMDA = [];
% % % for  i = 1 : nbr_seg
% % %      K_i = IMPULSION_BETA(i,5);
% % %      Vini_i = ENTRAINEMENT(i,3);
% % %      Vfin_i = ENTRAINEMENT(i,4);
% % %      V_max_i = K_i + ( ( Vini_i + Vfin_i ) / 2 );
% % %      
% % %      rapport_amplitude_impulsions_successive_vitesse_1 = max(V_max_i , Vini_i) / min(V_max_i , Vini_i);
% % %      grand_axe_a_1 = param_trajectoire(i,1);
% % %      petit_axe_b_1 = param_trajectoire(i,2);
% % %      rap_a_1_sur_b_1 = grand_axe_a_1 / petit_axe_b_1;
% % %      angle_delta_1 = abs( param_trajectoire(i,4) - param_trajectoire(i,6) );
% % %      coef_angle_delta_1 = sin( angle_delta_1 ) / ( sin( atan( rap_a_1_sur_b_1 * ( tan( angle_delta_1 ) ) ) ) );
% % %      lamda_1 = (1/3) * ( abs( ( log( rapport_amplitude_impulsions_successive_vitesse_1 ) ) / ( log( rap_a_1_sur_b_1 * coef_angle_delta_1 )  ) ) );
% % %
% % %      rapport_amplitude_impulsions_successive_vitesse_2 = max(V_max_i , Vfin_i) / min(V_max_i , Vfin_i);
% % %      grand_axe_a_2 = param_trajectoire(i,5);
% % %      petit_axe_b_2 = param_trajectoire(i,3);
% % %      rap_a_2_sur_b_2 = grand_axe_a_2 / petit_axe_b_2;
% % %      angle_delta_2 = abs( param_trajectoire(i,7) - param_trajectoire(i,4) );
% % %      coef_angle_delta_2 = sin( angle_delta_2 ) / ( sin( atan( rap_a_2_sur_b_2 * ( tan( angle_delta_2 ) ) ) ) );
% % %      lamda_2 = (1/3) * ( abs( ( log( rapport_amplitude_impulsions_successive_vitesse_2 ) ) / ( log( rap_a_2_sur_b_2 * coef_angle_delta_2 )  ) ) );
% % %      
% % %      lamda = (lamda_1 + lamda_2) / 2;
% % %      LAMDA = [LAMDA; lamda];
% % % end


%% A1 = param_trajectoire(:,1);
A1 = param_trajectoire(:,8);
%% B1 = param_trajectoire(:,2);
B1 = param_trajectoire(:,11);
%% B2 = param_trajectoire(:,3);
B2 = param_trajectoire(:,17);


%% TETA = param_trajectoire(:,4);
%% TETA = atan( tan( TETA ) );
%% TETA_P1 = param_trajectoire(:,6);
%% TETA_P1 = atan( tan( TETA_P1 ) );
%% TETA_P2 = param_trajectoire(:,7);
%% TETA_P2 = atan( tan( TETA_P2 ) );


teta_p1_deb_DAT = param_trajectoire(1,6);
teta_p1_deb = atan( tan( teta_p1_deb_DAT ) );

TETA = param_trajectoire(:,4) - (teta_p1_deb_DAT - teta_p1_deb);
TETA_P1 = param_trajectoire(:,6) - (teta_p1_deb_DAT - teta_p1_deb);
TETA_P2 = param_trajectoire(:,7) - (teta_p1_deb_DAT - teta_p1_deb);



%%modele_sortie = [T0, Tc, T1, P, K, Vini, Vfin, A1, B1, B2, TETA];
%modele_sortie = [DELTA_T, RAP_Tc, P, K, Vini, Vfin, A1, B1, B2, TETA];
modele_sortie = [DELTA_T, RAP_Tc, P, K, Vini, Vfin, RAP_Ki_SUR_Entrainement, A1, B1, B2, TETA_P1, TETA, TETA_P2, P0SITION_TRAIT];


modele_sortie = real( modele_sortie );


%save collecte_donnees modele_sortie;

        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%


        
        
        
        
        





                       
%%%%%                     [                   1         ,2         ,3         ,4   ,5         ,6      ,7      ,8      ,11      ,14    ,17     ,20       ,23 ,24 ,25 ,26 ]

%%%%% param_trajectoire = [param_trajectoire; grand_axe1,petit_axe1,petit_axe2,teta,grand_axe2,teta_p1,teta_p2,A1a1aa1,B1b1bb1,A2a2aa2,B2b2bb2,M_3KMH_HB,X10,Y10,X20,Y20];
