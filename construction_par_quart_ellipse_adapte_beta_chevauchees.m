function  [param_trajectoire, param_grandeur, e_indice, Vect_Erreurs, Sensibilite_plus, DAT] = construction_par_quart_ellipse_adapte_beta_chevauchees(fonctions_beta, tableau_affectation_intervalles, vect_X0XCX1, points, T, num_fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction construction_par_quart_ellipse_adapte permet de modéliser la trajectoire           %%%%%
%%%%% par des arcs d'ellipse. trois méthodes sont utilisés pour générer les arcs :                   %%%%%
%%%%% quart d'ellipse, cercle en projection oblique et deux points deux tangentes.                   %%%%%
%%%%% La fonction permet de calculer une variation continue de l'angle d'inclinaison                 %%%%%
%%%%% de la tangente au tracé permettant de mieux caractérisée la géométrie de la trajectoire.       %%%%%
%%%%% Entrées: instants de segmentation temporelle (approche impuls. Bêta chevauchées), trajectoire  %%%%%
%%%%% Sorties: paramètres elliptiques de modélisation de la trajectoire, Erreurs de reconstruction   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


e_indice = 1;
Sensibilite_plus = [];

maxi_x = max(points(:,1));% maxi_x est le maximum des abscisses
maxi_y = max(points(:,2));% maxi_y est le maximum des ordonnées
mini_x = min(points(:,1));% mini_x est le minimum des abscisses
mini_y = min(points(:,2));% mini_y est le minimum des ordonnées
dimension_Max = max(maxi_x - mini_x , maxi_y  - mini_y);


T0 = vect_X0XCX1(:,2);
XYm1 = [];
XYm2 = [];
param_ellipse = [];
param_trajectoire = [];
param_grandeur = [];
pas = 0.01;
nbr_intervalles = size(tableau_affectation_intervalles, 1);
indiceTemps_deb = tableau_affectation_intervalles(1,7);          %= 1;
indiceTemps_fin = tableau_affectation_intervalles(nbr_intervalles,8);

k1 = indiceTemps_deb;
k2 = indiceTemps_fin;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Angle d'inclinaison de la tangente au tracé, DAT. %%%%%%%%%%%%

DAT = [];
type_point = [];
Ttemps = [];
u = 2; %2;
dk = (2*u) + 1;
Lp = length(points);
pas = 0.01;

for k = k1 : k2,
    if (k - u < 1)
        cotion = (points(k + u,2) - points(1,2)) / (points(k + u,1) - points(1,1));
        dat_k = atan( cotion );
    elseif (k + u > Lp)
%        cotion = (points(Lp,2) - points(k - u,2)) / (points(Lp,1) - points(k - u,1));
        cotion = (points(Lp,2) - points(k-1 ,2)) / (points(Lp,1) - points(k-1 ,1));
        
        dat_k = atan( cotion );        
    else
        num = (points(k + u,2) - points(k - u,2));
        denum = (points(k + u,1) - points(k - u,1));
        cotion = num / denum;
        
        dat_k = atan( cotion );        
    end

    DAT = [DAT; dat_k];
    type_point = [type_point, 0];
    Ttemps = [Ttemps,k*pas];
    
    
end
memDAT = DAT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lhb = length(T0);
for i = 2 : Lhb %Lhb-1
    tempsT0 = T0(i);
    indice_i = round(tempsT0/pas) - k1 + 1;
    type_point(indice_i) = 1;                          % marquer les minimums de vitesse ou de rayon de courbure
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% éviter les discontinuités dans la variation de l'angle d'inclinaison de la tangente DAT

tour_ajout = 0;
L_DAT = length(DAT);
sens_preced = 0;

for i = 2 : L_DAT,
    tour_ajout = fix((DAT(i) - DAT(i-1))/(pi))*pi;
    DAT(i) = DAT(i) - tour_ajout;

    diff_angle_1 = abs( (DAT(i)) - DAT(i-1) );
    diff_angle_2 = abs( (DAT(i) + pi) - DAT(i-1) );
    diff_angle_3 = abs( (DAT(i) - pi) - DAT(i-1) );
    
    x = [diff_angle_1, diff_angle_2, diff_angle_3];
    diff_angle_min = min(x);
    if diff_angle_min == diff_angle_1
        DAT(i) = DAT(i);
        hu = 0;
    elseif diff_angle_min == diff_angle_2
        DAT(i) = DAT(i) + pi;
        hu = pi;
    elseif diff_angle_min == diff_angle_3
        DAT(i) = DAT(i) - pi;
        hu = -pi;
    end

    sens_actuel = sign( DAT(i) - DAT(i-1) );
    if (sens_actuel ~= sens_preced)&(sens_preced ~= 0)           %sign( DAT(i-1) - DAT(i-2) )
       if (type_point(i) == 1)|(type_point(i-1) == 1)            % minimum de vitesse ou de rayon de courbure
          indice_k = i + k1 - 1;
          [continuite,inversion] = verification_changement_sens_parcour(indice_k,points,T);
          if (continuite == 0)&(inversion == 1)
             DAT(i) = DAT(i) + (sens_preced*(pi));
             sens_actuel = sign( DAT(i) - DAT(i-1) );
          else
             DAT(i) = DAT(i);
          end          
       end
    end
    sens_preced = sens_actuel;
            
end
    
memDAT1 = DAT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ******************************************************* %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constitution des vecteurs caractéristiques des points extrémums XYm1 et XYm2

temps_T1_T2 = [];
uu = 0;%1;%2; %1;
for i = 1 : nbr_intervalles
    V_M_dep = tableau_affectation_intervalles(i,5);
    V_M_arr = tableau_affectation_intervalles(i,6);

    if V_M_dep <= V_M_arr 
       indiceM1 = tableau_affectation_intervalles(i,7);
       indiceM2 = tableau_affectation_intervalles(i,8);
    else
       indiceM1 = tableau_affectation_intervalles(i,8);
       indiceM2 = tableau_affectation_intervalles(i,7);
    end

    teta_M1 = DAT(indiceM1 - indiceTemps_deb + 1);
    
    if ( (indiceM2 - indiceTemps_deb + 1 - uu) > 0 ) & ( (indiceM2 - indiceTemps_deb + 1 + uu) < L_DAT )
       teta_M2 = ( DAT(indiceM2 - indiceTemps_deb + 1 + uu) + DAT(indiceM2 - indiceTemps_deb + 1 - uu) ) / 2;
    else
       teta_M2 = DAT(indiceM2 - indiceT0_deb + 1);
    end

    XYm1 = [XYm1; indiceM1 points(indiceM1,1) points(indiceM1,2) teta_M1];
    XYm2 = [XYm2; indiceM2 points(indiceM2,1) points(indiceM2,2) teta_M2];
    temps_T1_T2 = [temps_T1_T2; tableau_affectation_intervalles(i,1), tableau_affectation_intervalles(i,2)];

end


%%%%%%%%%%%%% initialisation des paramètres de mesures des erreurs de la reconstruction elliptique %%%%%%%%%%%%%
                                       ERR_Moy_D_q = 0;      %Erreur  moyenne en distance pour la méthode quart d'ellipse
                                       ERR_Max_D_q = 0;      %Erreur maximale en distance pour la méthode quart d'ellipse
                                       ERR_Moy_P_q = 0;      %Erreur  moyenne en pourcent pour la méthode quart d'ellipse
                                       ERR_Max_P_q = 0;      %Erreur maximale en pourcent pour la méthode quart d'ellipse
       
                                       ERR_Moy_D_a = 0;      %Erreur  moyenne en distance pour la méthode projection oblique des arcs d'ellipse
                                       ERR_Max_D_a = 0;      %Erreur maximale en distance pour la méthode projection oblique des arcs d'ellipse
                                       ERR_Moy_P_a = 0;      %Erreur  moyenne en pourcent pour la méthode projection oblique des arcs d'ellipse
                                       ERR_Max_P_a = 0;      %Erreur maximale en pourcent pour la méthode projection oblique des arcs d'ellipse

                                       ERR_Moy_D_t = 0;      %Erreur  moyenne en distance pour la méthode deux points deux tangentes
                                       ERR_Max_D_t = 0;      %Erreur maximale en distance pour la méthode deux points deux tangentes
                                       ERR_Moy_P_t = 0;      %Erreur  moyenne en pourcent pour la méthode deux points deux tangentes
                                       ERR_Max_P_t = 0;      %Erreur maximale en pourcent pour la méthode deux points deux tangentes
       
                                       nbr_points_total = 0; %longeur totale des arcs d'ellipse en points
                                       nbr_arc = 0;          %nombre total des arcs d'ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

KMH_deb = indiceTemps_deb;

M_2KMH = [];
M_2KMH_HB = [];
for i = 1 : nbr_intervalles
  KMH1 = XYm1(i,1);
  KMH2 = XYm2(i,1);

  M_2KMH_HB = [KMH1, KMH2];
  
  M_2Donnees = [XYm1(i,:), XYm2(i,:)];
    

%%%%%%%%%%% Détermination des parmètres du modèle elliptique (de repérage) a, b, teta %%%%%%%%%%%
    XM1 = XYm1(i,2); 
    YM1 = XYm1(i,3);
    XM2 = XYm2(i,2);
    YM2 = XYm2(i,3);

    if abs(KMH2 - KMH1) >= 1
       M_2points = [XM2, YM2, XM1, YM1];
       M_2KMH = [KMH2, KMH1];
       teta = XYm2(i,4);

%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%
%%% cas particulier de points ou la direction de la tangente au point M2 = M_G est confondue avec la direction (M1 M2)
       if atan( tan(teta) ) == atan( (YM2 - YM1) / (XM2 - XM1) )
          teta = teta + (pi / 1000);
       end
%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%

       [A1,B1,X10,Y10,Teta1,Drap1,Erreur_MAX1,Ordre_ok1,HA1] = parametres_quart_ellipse(M_2points,M_2KMH,teta,num_fig,1);
       Teta_G = teta;
       Teta_P = XYm1(i,4);
       teta_p1 = Teta_P;

%        if teta == Teta1;
%           Teta_G = XYm2(i,4);
%           Teta_P = XYm1(i,4);
%        else
%           Teta_G = XYm1(i,4);
%           Teta_P = XYm2(i,4);
%        end

       [a1,b1,x10,y10,teta1,drap1,erreur_MAX1,ordre_ok1,hA1] = parametres_quart_ellipse_oblique(A1,B1,X10,Y10,Teta_G,Teta_P,M_2points,M_2KMH,KMH_deb, DAT,num_fig,1);
       [aa1,bb1,xx10,yy10,tteta1,ddrap1,eerreur_MAX1,oordre_ok1,hhA1] = parametres_arc_ellipse_2pts2tg(A1,B1,X10,Y10,Teta_G,Teta_P,M_2points,M_2KMH,KMH_deb, DAT,num_fig,1);


                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%       
       liste_points = [];
       for indice_point = min(KMH1,KMH2) : max(KMH1,KMH2)
           liste_points = [liste_points; points(indice_point,1) points(indice_point,2)];
       end
       ortho_dist = [A1, B1];
       
       parametre_ellipse = [A1, B1, X10, Y10, Teta1];
       [Err_moy_d_q1, Err_Max_d_q1, Err_moy_p_q1, Err_Max_p_q1] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       parametre_ellipse = [a1, b1, x10, y10, teta1];
       [Err_moy_d_a1, Err_Max_d_a1, Err_moy_p_a1, Err_Max_p_a1] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       parametre_ellipse = [aa1, bb1, xx10, yy10, tteta1];
       [Err_moy_d_t1, Err_Max_d_t1, Err_moy_p_t1, Err_Max_p_t1] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);

       longueur_segment =  sqrt( ( (XM2-XM1)^2 ) + ( (YM2-YM1)^2 ) );
    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

    else
       A1 = 0; B1 = 0; X10 = XYm2(i,2); Y10 = XYm2(i,3); Teta1 = XYm2(i,4); teta = XYm2(i,4);
       a1 = 0; b1 = 0; x10 = XYm2(i,2); y10 = XYm2(i,3); teta1 = XYm2(i,4); teta = XYm2(i,4);
       aa1 = 0; bb1 = 0; xx10 = XYm2(i,2); yy10 = XYm2(i,3); tteta1 = XYm2(i,4); teta = XYm2(i,4);

       Err_moy_d_q1 = 0; Err_Max_d_q1 = 0; Err_moy_p_q1 = 0; Err_Max_p_q1 = 0;
       Err_moy_d_a1 = 0; Err_Max_d_a1 = 0; Err_moy_p_a1 = 0; Err_Max_p_a1 = 0;
       Err_moy_d_t1 = 0; Err_Max_d_t1 = 0; Err_moy_p_t1 = 0; Err_Max_p_t1 = 0;
       teta_p1 = XYm1(i,4);
   
       longueur_segment =  0.5; 
       teta_p1 = XYm1(i,4);

    end


                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Enregistrement des donnees geometriques \\%

  if  Err_moy_p_t1 < Err_moy_p_q1
      grand_axe1 = aa1;
      petit_axe1 = bb1;
  else
      grand_axe1 = A1;
      petit_axe1 = B1;
  end

%% param_trajectoire = [param_trajectoire; grand_axe1, petit_axe1, teta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% param_trajectoire complet pour la du Verification_rapport_vitesses_axes %%%%%%

A1a1aa1 = [A1, a1, aa1];
B1b1bb1 = [B1, b1, bb1];

%param_trajectoire = [param_trajectoire; grand_axe1, petit_axe1, teta, teta_p1, A1a1aa1, B1b1bb1, M_2KMH_HB ];
param_trajectoire = [param_trajectoire; grand_axe1,petit_axe1,teta,teta_p1,A1a1aa1,B1b1bb1,M_2KMH_HB];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

%  teta_limitee = atan( tan( teta ) );
%  param_trajectoire = [param_trajectoire; grand_axe1, petit_axe1, teta_limitee];
 
  
param_grandeur = [param_grandeur; A1, B1, longueur_segment];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Detection et enregistrement des erreurs de reconstruction geometriques intolerables
  if ( min( Err_moy_d_t1 , Err_moy_d_q1 ) / dimension_Max ) * 100 > 4.85%4.9%4.8%4.5%4%2%5%6%7%8
     Sensibilite_plus = [Sensibilite_plus; 1, temps_T1_T2(i,:)];    % requete d'augmenter la sensibilite
  end                                                               % de detection des points de double inflexion
                                                                    % dans l'intervalle de temps [T1(i) T2(i)]

% %// Detection et enregistrement des erreurs de reconstruction geometriques intolerables
%   if ( min( Err_Max_d_t1 , Err_Max_d_q1 ) / dimension_Max ) * 100 > 8%10%7%15
%      Sensibilite_plus = [Sensibilite_plus; 1, temps_T1_T2(i,:)];    % requete d'augmenter la sensibilite
%   end                                                               % de detection des points de double inflexion
%                                                                     % dans l'intervalle de temps [T0(i) Tc(i)]


% % %   if ( Err_moy_p_t1 < Err_moy_p_q1 )
% % %      param_trajectoire = [param_trajectoire; aa1, bb1, teta];
% % %      Err_moy_p_segment = Err_moy_p_t1;
% % %      Err_moy_Max_p_segment = Err_Max_p_t1;
% % %   else
% % %      param_trajectoire = [param_trajectoire; A1, B1, teta];
% % %      Err_moy_p_segment = Err_moy_p_q1;
% % %      Err_moy_Max_p_segment = Err_Max_p_q1;
% % %   end
% % % 
% % %   if (Err_moy_Max_p_segment > 15)%15) %(Err_moy_p_segment > 15)%20)
% % %      Sensibilite_plus = 1;          %requete d'augmenter la sensibilite de detection des points de double inflexion
% % %   end

%%  param_trajectoire = [param_trajectoire; A1, B1, teta, Teta_P];

%  if (A1 < B1)
%     e_indice = 0; % cas aberrant
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%
       nbr_points_M1M2 = abs(KMH2 - KMH1) + 1;

       if (abs(KMH2 - KMH1) >= 1)
          Err_moy_d_q = (Err_moy_d_q1*nbr_points_M1M2);
          Err_Max_d_q = (Err_Max_d_q1*nbr_points_M1M2);
          Err_moy_p_q = (Err_moy_p_q1*nbr_points_M1M2);
          Err_Max_p_q = (Err_Max_p_q1*nbr_points_M1M2);

          if (Err_moy_d_a1 < 1e+15)
             Err_moy_d_a = (Err_moy_d_a1*nbr_points_M1M2);
             Err_Max_d_a = (Err_Max_d_a1*nbr_points_M1M2);
             Err_moy_p_a = (Err_moy_p_a1*nbr_points_M1M2);
             Err_Max_p_a = (Err_Max_p_a1*nbr_points_M1M2);
          else
             Err_moy_d_a = Err_moy_d_q;
             Err_Max_d_a = Err_Max_d_q;
             Err_moy_p_a = Err_moy_p_q;
             Err_Max_p_a = Err_Max_p_q;
          end

          Err_moy_d_t = (Err_moy_d_t1*nbr_points_M1M2);
          Err_Max_d_t = (Err_Max_d_t1*nbr_points_M1M2);
          Err_moy_p_t = (Err_moy_p_t1*nbr_points_M1M2);
          Err_Max_p_t = (Err_Max_p_t1*nbr_points_M1M2);

       else
          Err_moy_d_q = 0;
          Err_Max_d_q = 0;
          Err_moy_p_q = 0;
          Err_Max_p_q = 0;

          Err_moy_d_a = 0;
          Err_Max_d_a = 0;
          Err_moy_p_a = 0;
          Err_Max_p_a = 0;

          Err_moy_d_t = 0;
          Err_Max_d_t = 0;
          Err_moy_p_t = 0;
          Err_Max_p_t = 0;

       end

       ERR_Moy_D_q = ERR_Moy_D_q + Err_moy_d_q;
       ERR_Max_D_q = ERR_Max_D_q + Err_Max_d_q;
       ERR_Moy_P_q = ERR_Moy_P_q + Err_moy_p_q;
       ERR_Max_P_q = ERR_Max_P_q + Err_Max_p_q;
       
       ERR_Moy_D_a = ERR_Moy_D_a + Err_moy_d_a;
       ERR_Max_D_a = ERR_Max_D_a + Err_Max_d_a;
       ERR_Moy_P_a = ERR_Moy_P_a + Err_moy_p_a;
       ERR_Max_P_a = ERR_Max_P_a + Err_Max_p_a;

       ERR_Moy_D_t = ERR_Moy_D_t + Err_moy_d_t;
       ERR_Max_D_t = ERR_Max_D_t + Err_Max_d_t;
       ERR_Moy_P_t = ERR_Moy_P_t + Err_moy_p_t;
       ERR_Max_P_t = ERR_Max_P_t + Err_Max_p_t;

       
       nbr_points_total = nbr_points_total + nbr_points_M1M2 - 1;
       nbr_arc = nbr_arc + 1;
       
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

end

                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%       

ERR_Moy_D_q = ERR_Moy_D_q  / nbr_points_total;
ERR_Max_D_q = ERR_Max_D_q  / nbr_points_total;
ERR_Moy_P_q = ERR_Moy_P_q  / nbr_points_total;
ERR_Max_P_q = ERR_Max_P_q  / nbr_points_total;
       
ERR_Moy_D_a = ERR_Moy_D_a  / nbr_points_total;
ERR_Max_D_a = ERR_Max_D_a  / nbr_points_total;
ERR_Moy_P_a = ERR_Moy_P_a  / nbr_points_total;
ERR_Max_P_a = ERR_Max_P_a  / nbr_points_total;

ERR_Moy_D_t = ERR_Moy_D_t  / nbr_points_total;
ERR_Max_D_t = ERR_Max_D_t  / nbr_points_total;
ERR_Moy_P_t = ERR_Moy_P_t  / nbr_points_total;
ERR_Max_P_t = ERR_Max_P_t  / nbr_points_total;


Vect_Erreurs = [ERR_Moy_D_q ERR_Max_D_q ERR_Moy_P_q ERR_Max_P_q;
                ERR_Moy_D_a ERR_Max_D_a ERR_Moy_P_a ERR_Max_P_a;
                ERR_Moy_D_t ERR_Max_D_t ERR_Moy_P_t ERR_Max_P_t];


                   %HBH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HBH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %HBH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HBH

