function  [param_trajectoire, param_grandeur, e_indice, Vect_Erreurs, Sensibilite_plus, DAT] = construction_par_arc_ellipse_adapte(IMPULSION_BETA,points,T, num_fig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction construction_par_quart_ellipse_adapte permet de modéliser la trajectoire           %%%%%
%%%%% par des arcs d'ellipse. trois méthodes sont utilisés pour générer les arcs :                   %%%%%
%%%%% quart d'ellipse, cercle en projection oblique et deux points deux tangentes.                   %%%%%
%%%%% La fonction permet de calculer une variation continue de l'angle d'inclinaison                 %%%%%
%%%%% de la tangente au tracé permettant de mieux caractérisée la géométrie de la trajectoire.       %%%%%
%%%%% Entrées: instants de segmentation temporelle (approche Bêta + comp. entrainement), trajectoire %%%%%
%%%%% Sorties: paramètres elliptiques de modélisation de la trajectoire, Erreurs de reconstruction   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


e_indice = 1;  % cas normal
Sensibilite_plus = [];

maxi_x = max(points(:,1));% maxi_x est le maximum des abscisses
maxi_y = max(points(:,2));% maxi_y est le maximum des ordonnées
mini_x = min(points(:,1));% mini_x est le minimum des abscisses
mini_y = min(points(:,2));% mini_y est le minimum des ordonnées
dimension_Max = max(maxi_x - mini_x , maxi_y  - mini_y);

T0 = IMPULSION_BETA(:,6);
T1 = IMPULSION_BETA(:,7);
Tc = IMPULSION_BETA(:,1);

XYm1 = [];
XYm3 = [];
XYm2 = [];

param_ellipse = [];
param_trajectoire = [];
param_grandeur = [];

%%% pas = 0.01;
pas = T(2) - T(1);

L = length(T0);

indiceT0_deb = round(T0(1)/pas);
if indiceT0_deb == 0
   indiceT0_deb = 1;
end
indiceT1_deb = round(T1(1)/pas) + 1;
indiceTc_deb = round(Tc(1)/pas) + 1;
indiceT1_fin = round(T1(L)/pas) + 1;

if indiceT1_fin > size(points,1)
   indiceT1_fin = size(points,1);
end


k1 = indiceT0_deb;
k2 = indiceT1_fin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Angle d'inclinaison de la tangente au tracé, DAT. %%%%%%%%%%%%

DAT = [];
type_point = [];
Ttemps = [];
u = 2; %2;
dk = (2*u) + 1;
Lp = length(points);
%pas = 0.01;
pas = T(2) - T(1);;

cotion_preced = 1;

for k = k1 : k2,
    
    if (k - u < 1)
        num = (points(k + u,2) - points(1,2));
        denum = (points(k + u,1) - points(1,1));
        if denum == 0
           denum = 1.0e-10;
        end
        %cotion = (points(k + u,2) - points(1,2)) / (points(k + u,1) - points(1,1));
        cotion = num / denum;
        dat_k = atan( cotion );
    elseif (k + u > Lp)
% %        cotion = (points(Lp,2) - points(k - u,2)) / (points(Lp,1) - points(k - u,1));
%         cotion = (points(Lp,2) - points(k-1 ,2)) / (points(Lp,1) - points(k-1 ,1));
        num = (points(Lp,2) - points(k-1 ,2));
        denum = (points(Lp,1) - points(k-1 ,1));
        if denum == 0
           denum = 1.0e-10;
        end
        cotion = num / denum;
        
        dat_k = atan( cotion );        
    else
        num = (points(k + u,2) - points(k - u,2));
        denum = (points(k + u,1) - points(k - u,1));
        
        if denum == 0
           denum = 1.0e-10;
           cotion = abs( num / denum ) * sign(cotion_preced);
        else
           cotion = num / denum;
        end
        cotion_preced = cotion;
                    
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

    sens_actuel = sign( DAT(i) - DAT(i-1) );  %1) );
    if (type_point(i) == 1)|(type_point(i-1) == 1)                 % minimum de vitesse ou de rayon de courbure
       indice_k = i + k1 - 1;
       [continuite,inversion] = verification_changement_sens_parcour(indice_k,points,T);        
       if (sens_actuel ~= sens_preced)&(sens_preced ~= 0)          %sign( DAT(i-1) - DAT(i-2) )
          if (continuite == 0)&(inversion == 1)
             DAT(i) = DAT(i) + (sens_preced*(pi));
             sens_actuel = sign( DAT(i) - DAT(i-1) );
          else
             DAT(i) = DAT(i);
          end

       else
          if (continuite == 0)&(inversion == 1)&( abs(DAT(i) - DAT(i-1)) < (pi/8) )&(type_point(i) == 1)
             DAT(i) = DAT(i) + (sens_preced*(pi/3));
             sens_actuel = sign( DAT(i) - DAT(i-1) );
          else
             DAT(i) = DAT(i);
          end
       end
       
    end

    sens_preced = sens_actuel;
            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(points,1) > size(DAT,1)
   indice_point_fin = size(points,1);
   for indice_point = size(DAT,1) + 1 : indice_point_fin
       DAT = [DAT; 0];
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
memDAT1 = DAT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ******************************************************* %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialisation de XYm1, XYm2, XYm3 les vecteurs caractéristiques des points extrémums

teta_M1 = DAT(indiceT0_deb - indiceT0_deb + 1);
teta_M3 = DAT(indiceT1_deb - indiceT0_deb + 1);
teta_M2 = DAT(indiceTc_deb - indiceT0_deb + 1);

XYm1 = [indiceT0_deb points(indiceT0_deb,1) points(indiceT0_deb,2) teta_M1];
XYm3 = [indiceT1_deb points(indiceT1_deb,1) points(indiceT1_deb,2) teta_M3];
XYm2 = [indiceTc_deb points(indiceTc_deb,1) points(indiceTc_deb,2) teta_M2];

% suite de la constitution des vecteurs caractéristiques des points extrémums
uu = 0;%1;%2; %1;
for i = 2 : L

    indiceT0 = round(T0(i)/pas) + 1;
    indiceT1 = round(T1(i)/pas) + 1;
    indiceTc = round(Tc(i)/pas) + 1;

    if (indiceT0 - indiceT0_deb + 1 - uu) > 0
       teta_M1 = ( DAT(indiceT0 - indiceT0_deb + 1) + DAT(indiceT0 - indiceT0_deb + 1 - uu) ) / 2;
    else
       teta_M1 = DAT(indiceT0 - indiceT0_deb + 1);
    end
    if (indiceT1 - indiceT0_deb + 1 + 1) < L_DAT
       teta_M3 = ( DAT(indiceT1 - indiceT0_deb + 1) + DAT(indiceT1 - indiceT0_deb + 1 + uu) ) / 2;
    else
       teta_M3 = DAT(indiceT1 - indiceT0_deb + 1);
    end

    teta_M2 = DAT(indiceTc - indiceT0_deb + 1);

    XYm1 = [XYm1; indiceT0 points(indiceT0,1) points(indiceT0,2) teta_M1];
    XYm3 = [XYm3; indiceT1 points(indiceT1,1) points(indiceT1,2) teta_M3];
    XYm2 = [XYm2; indiceTc points(indiceTc,1) points(indiceTc,2) teta_M2];          

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

KMH_deb = indiceT0_deb;

M_3KMH = [];
M_3KMH_HB = [];
for i = 1 : L
  KMH1 = XYm1(i,1);
  KMH3 = XYm3(i,1);
  KMH2 = XYm2(i,1);

  M_3KMH_HB = [KMH1, KMH2, KMH3];
  
  if  KMH3 - KMH1 >= 2

    M_3points(1) = XYm1(i,2);
    M_3points(2) = XYm1(i,3);
    M_3points(3) = XYm3(i,2);
    M_3points(4) = XYm3(i,3);
    M_3points(5) = XYm2(i,2);
    M_3points(6) = XYm2(i,3);
    M_3KMH = [XYm1(i,1),XYm3(i,1),XYm2(i,1)];
    M_3Donnees = [XYm1(i,:),XYm3(i,:),XYm2(i,:)];

%%%%%%%%%%% Détermination des parmètres du modèle elliptique (de repérage) a, b, teta %%%%%%%%%%%
    XM1 = XYm1(i,2); 
    YM1 = XYm1(i,3);
    XM3 = XYm3(i,2);
    YM3 = XYm3(i,3);
    XM2 = XYm2(i,2);
    YM2 = XYm2(i,3);
    
    if KMH2 - KMH1 > 0
       M_2points = [XM2, YM2, XM1, YM1];
       M_2KMH = [KMH2, KMH1];
       teta = XYm2(i,4);
       
%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%
%%% cas particulier de points ou la direction de la tangente au point M2 = M_G est confondue avec la direction (M1 M2)
       
       num_tan_teta_M2_M1 = (YM2 - YM1);
       denum_tan_teta_M2_M1 = (XM2 - XM1);
       if denum_tan_teta_M2_M1 == 0
          denum_tan_teta_M2_M1 = 1.0e-10;
       end
       %teta_M2_M1 = atan( (YM2 - YM1) / (XM2 - XM1) )
       teta_M2_M1 = atan( num_tan_teta_M2_M1 / denum_tan_teta_M2_M1 );
       
       %if atan( tan(teta) ) == atan( (YM2 - YM1) / (XM2 - XM1) )
       if atan( tan(teta) ) == teta_M2_M1
          teta = teta + (pi / 1000);
       end
%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%
       
       autorisation_de_representation_graphique_de_l_arc = 0;
       
       [A1,B1,X10,Y10,Teta1,Drap1,Erreur_MAX1,Ordre_ok1,HA1] = parametres_quart_ellipse(M_2points, M_2KMH, teta, num_fig,autorisation_de_representation_graphique_de_l_arc);       
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

       [a1,b1,x10,y10,teta1,drap1,erreur_MAX1,ordre_ok1,hA1] = parametres_quart_ellipse_oblique(A1, B1, X10, Y10, Teta_G, Teta_P, M_2points, M_2KMH, KMH_deb, DAT, num_fig, autorisation_de_representation_graphique_de_l_arc);
       [aa1,bb1,xx10,yy10,tteta1,ddrap1,eerreur_MAX1,oordre_ok1,hhA1] = parametres_arc_ellipse_2pts2tg(A1, B1, X10, Y10, Teta_G, Teta_P, M_2points, M_2KMH,KMH_deb, DAT, num_fig, autorisation_de_representation_graphique_de_l_arc);


                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%       
       liste_points = [];
       for indice_point = KMH1 : KMH2
           liste_points = [liste_points; points(indice_point,1) points(indice_point,2)];
       end
       ortho_dist = [A1, B1];
       
       parametre_ellipse = [A1, B1, X10, Y10, Teta1];
       [Err_moy_d_q1, Err_Max_d_q1, Err_moy_p_q1, Err_Max_p_q1] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       parametre_ellipse = [a1, b1, x10, y10, teta1];
       [Err_moy_d_a1, Err_Max_d_a1, Err_moy_p_a1, Err_Max_p_a1] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       parametre_ellipse = [aa1, bb1, xx10, yy10, tteta1];
       [Err_moy_d_t1, Err_Max_d_t1, Err_moy_p_t1, Err_Max_p_t1] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

    else
       A1 = 0; B1 = 0; X10 = XYm2(i,2); Y10 = XYm2(i,3); Teta1 = XYm2(i,4); teta = XYm2(i,4);
       a1 = 0; b1 = 0; x10 = XYm2(i,2); y10 = XYm2(i,3); teta1 = XYm2(i,4); teta = XYm2(i,4);
       aa1 = 0; bb1 = 0; xx10 = XYm2(i,2); yy10 = XYm2(i,3); tteta1 = XYm2(i,4); teta = XYm2(i,4);

       Err_moy_d_q1 = 0; Err_Max_d_q1 = 0; Err_moy_p_q1 = 0; Err_Max_p_q1 = 0;
       Err_moy_d_a1 = 0; Err_Max_d_a1 = 0; Err_moy_p_a1 = 0; Err_Max_p_a1 = 0;
       Err_moy_d_t1 = 0; Err_Max_d_t1 = 0; Err_moy_p_t1 = 0; Err_Max_p_t1 = 0;
       teta_p1 = XYm1(i,4);
    end

    if KMH3 - KMH2 > 0
       M_2points = [XM2, YM2, XM3, YM3];
       M_2KMH = [KMH2, KMH3];
       teta = XYm2(i,4);
       
%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%
%%% cas particulier de points ou la direction de la tangente au point M3 = M_G est confondue avec la direction (M2 M3)
       
       num_tan_teta_M3_M2 = (YM3 - YM2);
       denum_tan_teta_M3_M2 = (XM3 - XM2);
       if denum_tan_teta_M3_M2 == 0
          denum_tan_teta_M3_M2 = 1.0e-10;
       end
       %teta_M3_M2 = atan( (YM3 - YM2) / (XM3 - XM2) )
       teta_M3_M2 = atan( num_tan_teta_M3_M2 / denum_tan_teta_M3_M2 );
       
       %if atan( tan(teta) ) == atan( (YM3 - YM2) / (XM3 - XM2) )
       if atan( tan(teta) ) == teta_M3_M2
          teta = teta + (pi / 1000);
       end
%%%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^%%%
       
       
       autorisation_de_representation_graphique_de_l_arc = 0;

       [A2,B2,X20,Y20,Teta2,Drap2,Erreur_MAX2,Ordre_ok2,HA2] = parametres_quart_ellipse(M_2points, M_2KMH, teta, num_fig, autorisation_de_representation_graphique_de_l_arc);

       Teta_G = teta;
       Teta_P = XYm3(i,4);
       teta_p2 = Teta_P;

%        if teta == Teta2;
%           Teta_G = XYm2(i,4);
%           Teta_P = XYm3(i,4);
%        else
%           Teta_G = XYm3(i,4);
%           Teta_P = XYm2(i,4);
%        end

       [a2,b2,x20,y20,teta2,drap2,erreur_MAX2,ordre_ok2,hA2] = parametres_quart_ellipse_oblique(A2, B2, X20, Y20, Teta_G, Teta_P, M_2points, M_2KMH, KMH_deb, DAT, num_fig, autorisation_de_representation_graphique_de_l_arc);
       [aa2,bb2,xx20,yy20,tteta2,ddrap2,eerreur_MAX2,oordre_ok2,hhA2] = parametres_arc_ellipse_2pts2tg(A2, B2, X20, Y20, Teta_G, Teta_P, M_2points, M_2KMH,KMH_deb, DAT, num_fig, autorisation_de_representation_graphique_de_l_arc);

                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%
       liste_points = [];
       for indice_point = KMH2 : KMH3
           liste_points = [liste_points; points(indice_point,1) points(indice_point,2)];
       end
       ortho_dist = [A2, B2];
       
       parametre_ellipse = [A2, B2, X20, Y20, Teta2];
       [Err_moy_d_q2, Err_Max_d_q2, Err_moy_p_q2, Err_Max_p_q2] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       parametre_ellipse = [a2, b2, x20, y20, teta2];
       [Err_moy_d_a2, Err_Max_d_a2, Err_moy_p_a2, Err_Max_p_a2] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);
       parametre_ellipse = [aa2, bb2, xx20, yy20, tteta2];
       [Err_moy_d_t2, Err_Max_d_t2, Err_moy_p_t2, Err_Max_p_t2] = mesure_erreur_reconstruction_arc(liste_points,parametre_ellipse,ortho_dist);

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else
       A2 = 0; B2 = 0; X20 = XYm2(i,2); Y20 = XYm2(i,3); Teta2 = XYm2(i,4); teta = XYm2(i,4);
       a2 = 0; b2 = 0; x20 = XYm2(i,2); y20 = XYm2(i,3); teta2 = XYm2(i,4); teta = XYm2(i,4);
       aa2 = 0; bb2 = 0; xx20 = XYm2(i,2); yy20 = XYm2(i,3); tteta2 = XYm2(i,4); teta = XYm2(i,4);

       Err_moy_d_q2 = 0; Err_Max_d_q2 = 0; Err_moy_p_q2 = 0; Err_Max_p_q2 = 0;
       Err_moy_d_a2 = 0; Err_Max_d_a2 = 0; Err_moy_p_a2 = 0; Err_Max_p_a2 = 0;
       Err_moy_d_t2 = 0; Err_Max_d_t2 = 0; Err_moy_p_t2 = 0; Err_Max_p_t2 = 0;
       teta_p2 = XYm3(i,4);
    end
  
    longueur_segment =  (sqrt( ( (XM2-XM1)^2 ) + ( (YM2-YM1)^2 ) )) + (sqrt( ( (XM3-XM2)^2 ) + ( (YM3-YM2)^2 ) ));
    
  else
    A1 = 0; B1 = 0; X10 = XYm1(i,2); Y10 = XYm1(i,3); Teta1 = XYm1(i,4); teta = XYm1(i,4);
    A2 = 0; B2 = 0; X20 = XYm3(i,2); Y20 = XYm3(i,3); Teta2 = XYm3(i,4); teta = XYm1(i,4);

    a1 = 0; b1 = 0; x10 = XYm2(i,2); y10 = XYm2(i,3); teta1 = XYm2(i,4); teta = XYm2(i,4);
    a2 = 0; b2 = 0; x20 = XYm2(i,2); y20 = XYm2(i,3); teta2 = XYm2(i,4); teta = XYm2(i,4);

    aa1 = 0; bb1 = 0; xx10 = XYm2(i,2); yy10 = XYm2(i,3); tteta1 = XYm2(i,4); teta = XYm2(i,4);
    aa2 = 0; bb2 = 0; xx20 = XYm2(i,2); yy20 = XYm2(i,3); tteta2 = XYm2(i,4); teta = XYm2(i,4);

    Err_moy_d_q1 = 0; Err_Max_d_q1 = 0; Err_moy_p_q1 = 0; Err_Max_p_q1 = 0;
    Err_moy_d_a1 = 0; Err_Max_d_a1 = 0; Err_moy_p_a1 = 0; Err_Max_p_a1 = 0;
    Err_moy_d_t1 = 0; Err_Max_d_t1 = 0; Err_moy_p_t1 = 0; Err_Max_p_t1 = 0;
    
    Err_moy_d_q2 = 0; Err_Max_d_q2 = 0; Err_moy_p_q2 = 0; Err_Max_p_q2 = 0;
    Err_moy_d_a2 = 0; Err_Max_d_a2 = 0; Err_moy_p_a2 = 0; Err_Max_p_a2 = 0;
    Err_moy_d_t2 = 0; Err_Max_d_t2 = 0; Err_moy_p_t2 = 0; Err_Max_p_t2 = 0;
    
    longueur_segment =  0.5;     % minimale ( par défaut )
    
    teta_p1 = XYm1(i,4);
    teta_p2 = XYm3(i,4);

  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Enregistrement des donnees geometriques \\%

  if  Err_moy_p_t1 < Err_moy_p_q1
      grand_axe1 = aa1;
      petit_axe1 = bb1;
  else
      grand_axe1 = A1;
      petit_axe1 = B1;
  end
  if  Err_moy_p_t2 < Err_moy_p_q2
      grand_axe2 = aa2;
      petit_axe2 = bb2;
  else
      grand_axe2 = A2;
      petit_axe2 = B2;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% param_trajectoire complet pour la du Verification_rapport_vitesses_axes %%%%%%
A1a1aa1 = [A1, a1, aa1];
B1b1bb1 = [B1, b1, bb1];
A2a2aa2 = [A2, a2, aa2];
B2b2bb2 = [B2, b2, bb2];
%param_trajectoire = [param_trajectoire; grand_axe1,petit_axe1,petit_axe2,teta,grand_axe2,teta_p1,teta_p2,A1a1aa1,B1b1bb1,A2a2aa2,B2b2bb2,M_3KMH_HB];
param_trajectoire = [param_trajectoire; grand_axe1,petit_axe1,petit_axe2,teta,grand_axe2,teta_p1,teta_p2,A1a1aa1,B1b1bb1,A2a2aa2,B2b2bb2,M_3KMH_HB,X10,Y10,X20,Y20];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
param_grandeur = [param_grandeur; A1, B1, A2, B2, longueur_segment];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%// Detection et enregistrement des erreurs de reconstruction geometriques intolerables

  if ( min( Err_moy_d_t1 , Err_moy_d_q1 ) / dimension_Max ) * 100 > 4.85%4.9%4.8%4.5%4%2%5%6%7%8
     Sensibilite_plus = [Sensibilite_plus; 1 T0(i) Tc(i)];          % requete d'augmenter la sensibilite
  end                                                               % de detection des points de double inflexion
                                                                    % dans l'intervalle de temps [T0(i) Tc(i)]
  if ( min( Err_moy_d_t2 , Err_moy_d_q2 ) / dimension_Max ) * 100 > 4.85%4.9%4.8%4.5%4%2%5%6%7%8
     Sensibilite_plus = [Sensibilite_plus; 1 Tc(i) T1(i)];          % requete d'augmenter la sensibilite
  end                                                               % de detection des points de double inflexion
                                                                    % dans l'intervalle de temps [Tc(i) T1(i)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%
                    
       nbr_points_M1M2 = KMH2 - KMH1 + 1;
       nbr_points_M2M3 = KMH3 - KMH2 + 1;

       if (KMH3 - KMH1 >= 2)
          Err_moy_d_q = (Err_moy_d_q1*nbr_points_M1M2) + (Err_moy_d_q2*nbr_points_M2M3);
          Err_Max_d_q = (Err_Max_d_q1*nbr_points_M1M2) + (Err_Max_d_q2*nbr_points_M2M3);
          Err_moy_p_q = (Err_moy_p_q1*nbr_points_M1M2) + (Err_moy_p_q2*nbr_points_M2M3);
          Err_Max_p_q = (Err_Max_p_q1*nbr_points_M1M2) + (Err_Max_p_q2*nbr_points_M2M3);
          
          
          if (Err_moy_d_a1 < 1e+15)&(Err_moy_d_a2 < 1e+15)
             Err_moy_d_a = (Err_moy_d_a1*nbr_points_M1M2) + (Err_moy_d_a2*nbr_points_M2M3);
             Err_Max_d_a = (Err_Max_d_a1*nbr_points_M1M2) + (Err_Max_d_a2*nbr_points_M2M3);
             Err_moy_p_a = (Err_moy_p_a1*nbr_points_M1M2) + (Err_moy_p_a2*nbr_points_M2M3);
             Err_Max_p_a = (Err_Max_p_a1*nbr_points_M1M2) + (Err_Max_p_a2*nbr_points_M2M3);
          else
             Err_moy_d_a = Err_moy_d_q;
             Err_Max_d_a = Err_Max_d_q;
             Err_moy_p_a = Err_moy_p_q;
             Err_Max_p_a = Err_Max_p_q;
          end
             

          Err_moy_d_t = (Err_moy_d_t1*nbr_points_M1M2) + (Err_moy_d_t2*nbr_points_M2M3);
          Err_Max_d_t = (Err_Max_d_t1*nbr_points_M1M2) + (Err_Max_d_t2*nbr_points_M2M3);
          Err_moy_p_t = (Err_moy_p_t1*nbr_points_M1M2) + (Err_moy_p_t2*nbr_points_M2M3);
          Err_Max_p_t = (Err_Max_p_t1*nbr_points_M1M2) + (Err_Max_p_t2*nbr_points_M2M3);

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

       
       nbr_points_total = nbr_points_total + nbr_points_M1M2 + nbr_points_M2M3 - 1;
       nbr_arc = nbr_arc + 1;
       
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

end

                    %%%%%%%%%%%%%%%%% calcul des erreurs de reconstruction %%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%% moyennes enregistrées pour le caractère entier %%%%%%%%%%%%

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

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                   %HB% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %HB%
           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   %HB% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %HB%

    
    
    