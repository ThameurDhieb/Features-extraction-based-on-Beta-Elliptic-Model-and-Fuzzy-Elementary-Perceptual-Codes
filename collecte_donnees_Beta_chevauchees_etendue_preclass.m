function [modele_sortie, KMH_deb_fin_traits] = collecte_donnees_Beta_chevauchees_etendue_preclass(fonctions_beta, param_trajectoire)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction collecte_donnees_Beta_chevauchees permet de collecter   %%%%%
%%%%% les paramètre des impulsions Bêta chevauchées ainsi que ceux et des %%%%%
%%%%% arcs d'ellipse modélisant le signal en ligne dans une seule matrice %%%%%
%%%%% qui représente le modèle de sortie                                  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nbr_seg = size(fonctions_beta, 1);

T0 = fonctions_beta(:,1);
Tc = fonctions_beta(:,2);
T1 = fonctions_beta(:,3);
P = fonctions_beta(:,4);
K = fonctions_beta(:,6);

DELTA_T = [];
RAP_Tc = [];
P0SITION_TRAIT = [];
Rapport_Amplitude_Impulsions_Successives_Vitesse = [];
for  i = 1 : nbr_seg
     delta_T = T1(i) - T0(i);
     DELTA_T = [DELTA_T; delta_T];
     rap_Tc = (Tc(i) - T0(i)) / delta_T;
     RAP_Tc = [RAP_Tc; rap_Tc];
     
     
     position_trait = 1;
     if (i == 1)&(nbr_seg >= 2)
        position_trait = 1;                 %%%%% premier trait au début de la trajectoire
     elseif (i == 2)&(nbr_seg >= 3)
        position_trait = 2;                 %%%%% 2 ème trait depuis le début de la trajectoire
     elseif (i == 3)&(nbr_seg >= 6)
        position_trait = 3;                 %%%%% 3 ème trait depuis le début de la trajectoire
     elseif (i == 4)&(nbr_seg >= 8)
        position_trait = 4;                 %%%%% 4 ème trait depuis le début de la trajectoire
     elseif (i > 4) & (i < (nbr_seg-3)) & (nbr_seg >= 9)
        position_trait = 5;                 %%%%% trait au milieu de la trajectoire
     elseif (i > 4) & (i == (nbr_seg-3))
        position_trait = 6;                 %%%%% avant avant avant dernier trait à la fin de la trajectoire
     elseif (i > 3) & (i == (nbr_seg-2))
        position_trait = 7;                 %%%%% avant avant dernier trait à la fin de la trajectoire
     elseif (i > 2) & (i == (nbr_seg-1))
        position_trait = 8;                 %%%%% avant dernier trait à la fin de la trajectoire
     elseif (i > 1) & (i == nbr_seg)
        position_trait = 9;                 %%%%% dernier trait à la fin de la trajectoire
     elseif (i == 1) & (i == nbr_seg)
        position_trait = 0.1;               %%%%% trait isolé
     end
     P0SITION_TRAIT = [P0SITION_TRAIT; position_trait];


     K_i = fonctions_beta(i,6);
     %% if (i < nbr_seg)
     %%    K_i_plus_1 = fonctions_beta(i+1,6);
     %%    %rapport_amplitude_impulsions_successives_vitesse = K_i / K_i_plus_1;
     %%    rapport_amplitude_impulsions_successives_vitesse = max(K_i , K_i_plus_1) / min(K_i , K_i_plus_1);
     %% else
     %%    rapport_amplitude_impulsions_successives_vitesse = 1;
     %% end
     if (i > 1)
        K_i_moins_1 = fonctions_beta(i-1,6);
        if (K_i_moins_1 == 0)
           K_i_moins_1 = 1.0e-10;
        end
        %% rapport_amplitude_impulsions_successives_vitesse = K_i / K_i_moins_1;
        %% rapport_amplitude_impulsions_successives_vitesse = (K_i + 1) / (K_i_moins_1 + 1);
        rapport_amplitude_impulsions_successives_vitesse = (K_i + 0.1) / (K_i_moins_1 + 0.1);
        %% rapport_amplitude_impulsions_successives_vitesse = max(K_i , K_i_moins_1) / min(K_i , K_i_moins_1);
     else
        K_i_moins_1 = 1.0e-10;
        %% rapport_amplitude_impulsions_successives_vitesse = (K_i + 1) / (K_i_moins_1 + 1);
        rapport_amplitude_impulsions_successives_vitesse = (K_i + 0.1) / (K_i_moins_1 + 0.1);
        %% rapport_amplitude_impulsions_successives_vitesse = 2.5;
     end
     Rapport_Amplitude_Impulsions_Successives_Vitesse = [Rapport_Amplitude_Impulsions_Successives_Vitesse; rapport_amplitude_impulsions_successives_vitesse];
end

RAP_Ki_SUR_Kimoins1 = Rapport_Amplitude_Impulsions_Successives_Vitesse;



% % % LAMDA = [];
% % % for  i = 1 : nbr_seg
% % %      K_i = fonctions_beta(i,6);
% % %      if (i < nbr_seg)
% % %         K_i_plus_1 = fonctions_beta(i+1,6);
% % %         rapport_amplitude_impulsions_successives_vitesse = max(K_i , K_i_plus_1) / min(K_i , K_i_plus_1);
% % %         grand_axe_a = param_trajectoire(i,1);
% % %         petit_axe_b = param_trajectoire(i,2);
% % %         rap_a_sur_b = grand_axe_a / petit_axe_b;
% % %         angle_delta = abs( param_trajectoire(i,4) - param_trajectoire(i,3) );
% % %         coef_angle_delta = sin( angle_delta ) / ( sin( atan( rap_a_sur_b * ( tan( angle_delta ) ) ) ) );
% % %         lamda = (1/3) * ( abs( ( log( rapport_amplitude_impulsions_successives_vitesse ) ) / ( log( rap_a_sur_b * coef_angle_delta )  ) ) );
% % %         
% % %      else
% % %         lamda = 1/3;
% % %      end
% % %      LAMDA = [LAMDA; lamda];
% % % end


%% A = param_trajectoire(:,1);
A = param_trajectoire(:,5);
%% B = param_trajectoire(:,2);
B = param_trajectoire(:,8);




%% TETA = param_trajectoire(:,3);
%% TETA = atan( tan( TETA ) );
%% 
%% TETA_P = param_trajectoire(:,4);
%% TETA_P = atan( tan( TETA_P ) );

teta_deb_DAT = param_trajectoire(1,3);
teta_deb = atan( tan( teta_deb_DAT ) );

TETA = param_trajectoire(:,3) - (teta_deb_DAT - teta_deb);
TETA_P = param_trajectoire(:,4) - (teta_deb_DAT - teta_deb);




A1 = A(1 : length(A)-1); B1 = B(1 : length(B)-1); TETA1 = TETA(1 : length(TETA)-1);


%%% modele_sortie = [DELTA_T, RAP_Tc, P, K, A1, B1, TETA1];

%modele_sortie = [DELTA_T, RAP_Tc, P, K, A, B, TETA];

%modele_sortie = [DELTA_T, RAP_Tc, P, K, RAP_Ki_SUR_Kipls1, A, B, TETA, TETA_P, P0SITION_TRAIT];
modele_sortie = [DELTA_T, RAP_Tc, P, K, RAP_Ki_SUR_Kimoins1, A, B, TETA, TETA_P, P0SITION_TRAIT];

modele_sortie = real( modele_sortie );


%save collecte_donnees modele_sortie;

KMH_deb_fin_traits = [];

for  i = 1 : nbr_seg
     KMH_M1 = param_trajectoire(i,11);
     KMH_M2 = param_trajectoire(i,12);
     KMH_deb_trait_courant = min([KMH_M1 , KMH_M2]);
     KMH_fin_trait_courant = max([KMH_M1 , KMH_M2]);
     KMH_deb_fin_traits = [KMH_deb_fin_traits; KMH_deb_trait_courant KMH_fin_trait_courant];
end



        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%




%%%%%                     [                   1         , 2         , 3   , 4      , 5      , 8      , 11       , 13 , 14 ];

        
%%%%% param_trajectoire = [param_trajectoire; grand_axe1, petit_axe1, teta, teta_p1, A1a1aa1, B1b1bb1, M_2KMH_HB, X10, Y10];


