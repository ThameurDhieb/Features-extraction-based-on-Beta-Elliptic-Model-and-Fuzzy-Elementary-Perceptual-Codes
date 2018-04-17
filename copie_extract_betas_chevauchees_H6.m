function [fonctions_beta_sortie, tableau_affectation_intervalles_sortie, vect_X0XCX1]=copie_extract_betas_chevauchees_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus, pas)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction extract_beta calcule les  5 parametres caracteristiques des differents beta chevauchées  %%%%%
%%%%% contenues dans un signal de vitesse défini                                                           %%%%%
%%%%% les parametres d'entree etant le signal de vitesse curviligne  filtré, lui meme                      %%%%%
%%%%% les paramètres de sortie sont les caractéristiques des Betas et de la composante d'entrainement      %%%%%
%%%%% constituant le signal de vitesse curviligne à savoir le t0, tc, t1, K_hauteur, p et q                %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ajout_impulsion_finale_nulle = 'oui';
%ajout_impulsion_finale_nulle = 'non';


res = [];
res1 = [];

resVL = [];
resVPorteusePM = [];
resVimpulsive_betaPx = [];
resVreconstruite =[];

result1 = [];
Nbre_inff = 0;
Nbre_inf = 0;
extremum = [];
XC = [];
X0 = [];

X1 = [];
P = [];
Q = [];
DELTA = [];
HAUTEUR = [];
ENTRAINEMENT = [];

vxmax = [];
vymax = [];
ind_vmin = [];
ind_vmax = [];
v_min = [];
v_max = [];
vxmin = [];
vymin = [];
Nbre_max = [];
extrema = [];

YY = [];
XX = [];

%%%%%%%%
temp_depart = T(1);
for i = 1 : length(T)
    T(i) = T(i) - temp_depart;
end
for kor = 1 : size( Sensibilite_plus, 1 )
    Sensibilite_plus( kor, 2 ) = Sensibilite_plus( kor, 2 ) - temp_depart;
    Sensibilite_plus( kor, 3 ) = Sensibilite_plus( kor, 3 ) - temp_depart;
end
%%%%%%%%

XX = T;
YY = V;
ZZ = DPV;

%%% pas = 0.01;

%-------------------------------------derivee seconde de la vitesse
delta_t = diff(t);
delta_DPV = diff(DPV);
DSV = (delta_DPV)./(delta_t);
DSV = [0; DSV];
%------------------------------------------------------------------
%-------------------------------------longueur du vecteur vitesse
l = length(YY);
% ----------------------------------- initialisation des paramètres
vymin = [vymin; YY(1)];
vxmin = [vxmin; XX(1)];
ind_vmin = [ind_vmin; (XX(1)/pas+1)];
v_min = [v_min; ( max(round(XX(1)/pas),1) ) XX(1) YY(1)] ;
extrema = [extrema; ( max(round(XX(1)/pas),1) ) XX(1) YY(1) 0];
   
%----------------------------------- detection des points d'inflexion des minima et des maxima
               
%----------------------------------- detection des maxima

i = 3;
while i <= l

  if YY(i-1) >= YY(i-2)
     vxmax = vxmax;
     vymax = vymax;
     if YY(i-1) >= YY(i)
        vymax = [vymax; YY(i-1)];
        vxmax = [vxmax; XX(i-1)];
        ind_vmax = [ind_vmax; (round(XX(i-1)/pas) + 1)];
        v_max = [ v_max; (round(XX(i-1)/pas) + 1) XX(i-1) YY(i-1)];
        extrema = [extrema; (round(XX(i-1)/pas) + 1) XX(i-1) YY(i-1) 2];
      end;
  end;

%---------------------------------- detection des minima

  if YY(i-1) <= YY(i-2)
     vxmin = vxmin;
     vymin = vymin;
     if YY(i-1) <= YY(i)
        vymin = [vymin; YY(i-1)];
        vxmin = [vxmin; XX(i-1)];
        ind_vmin = [ind_vmin; (round(XX(i-1)/pas) + 1)];
        v_min = [v_min; (round(XX(i-1)/pas) + 1) XX(i-1) YY(i-1)];
         
        extrema = [extrema; (round(XX(i-1)/pas) + 1) XX(i-1) YY(i-1) 0];
     end;
  end;

  i=i+1;   
 
end;

%-------------------------------------- le dernier points 

v_min = [v_min; (round(XX(i-1)/pas) + 1) XX(l) YY(l)];
extrema = [extrema; (round(XX(i-1)/pas) + 1) XX(l) YY(l) 0];
                         
%--------------------------------------------------------------------------
 
p = 2; 
i = 1;
 
 %--------------------- Elimination des points voisins ( bruit )
extremum = [];
for i = 1 : size(extrema,1)
    extremum = [extremum; extrema(i,:)];
    hh = extremum;
end   
     
j=1;
  
%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%

while (j <= size(extremum,1) - 2)

  if (extremum(j,4) == 0)
      uv = 1;
      while (extremum(j+uv,4) ~= 2)&(j+uv <= size(extremum,1) - 1)
          uv = uv + 1;
      end    
      if (extremum(j+uv,4) == 2)&(j+uv <= size(extremum,1) - 1)
         uvp = 1;
         while (extremum(j+uv+uvp,4) ~= 0)&(j+uv+uvp <= size(extremum,1))
            uvp = uvp + 1;
         end    
         if (extremum(j+uv+uvp,4) == 0)&(j+uv+uvp <= size(extremum,1))
            X0 = [X0; extremum(j,:)];
            XC = [XC;extremum(j+uv,:)];
            X1 = [X1;extremum(j+uv+uvp,:)];
         end
      end
     j = j+uv+uvp - 1;
  end    

    j = j + 1;
end

%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vect_X0XCX1 = [X0, XC, X1];

%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
% mettre dans l'ordre les indices des nouveaux points détectés %
Suite_Extremums_Ordonnee = [];
oui = 1;
non = 0;
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%


                %HB% %%%%%%%%%%%%%% Detection des points de distortion %%%%%%%%%%%%%%% %HB%
           %HB% %%%%%%%%%%%%%% Detection des points d'inflexion particuliers %%%%%%%%%%%%%%% %HB%

nbr_intervalles_Sensibilite_plus = size( Sensibilite_plus, 1);

Rap_lim1 = RAP_lim(1);      % seuil de detection grossiere des points de double inflexion en montée 
Rap_lim2 = RAP_lim(2);      % seuil de detection grossiere des points de double inflexion en descente 

           
j = 1;
nbr_int_init = length(XC(:,1));
nbr_int_fin = nbr_int_init;
l=length(YY);
l_DPV = length(DPV);
delta_i = 1;%2;
delta_i2 = 2; %4; %1;

while j <= nbr_int_fin,
   panier_X5 = [];
   panier_X6 = [];
 
   DPV_MMAX_LOCAL = 0;
   DPV_MMIN_LOCAL = 0;
   marque_fin_DPV = 0;

%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
   Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; X0(j,:)];
   maximum_inscrit = non;
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

   for ih = X0(j,1)-X0(1,1)+1 : X1(j,1)-X0(1,1)

%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
       if ( ih >= XC(j,1)-X0(1,1)+1 )&( maximum_inscrit == non )
          Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; XC(j,:)];
          maximum_inscrit = oui;
       end
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%    **************************************    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
       if nbr_intervalles_Sensibilite_plus > 0
          Rap_lim1 = RAP_lim(1);          % seuil de detection grossiere des points de double inflexion en montée 
          Rap_lim2 = RAP_lim(2);          % seuil de detection grossiere des points de double inflexion en descente 

          for kor = 1 : nbr_intervalles_Sensibilite_plus
              temps_lim1 = Sensibilite_plus( kor, 2 );
              temps_lim2 = Sensibilite_plus( kor, 3 );
              if (XX(ih) > temps_lim1)&(XX(ih) < temps_lim2)
                 Rap_lim1 = RAP_lim(3);   % seuil de detection fine des points de double inflexion en montée 
                 Rap_lim2 = RAP_lim(4);   % seuil de detection fine des points de double inflexion en descente 
              end
          end
       end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%    **************************************    %%%%%%%%%%%%%%%%%%%%%%%%%%%%          
       
       if (ih-delta_i > 0)&(ih+delta_i2 < l_DPV)

%HB% Detection des points (d'inflexion excédentaires+) de double inflexion , codés [5] %%%%%%%%%%%%
         if (DPV(ih) >= 0)
            if (DPV(ih) > DPV(ih-delta_i))&(DPV(ih) > DPV(ih+delta_i))
               dernier_DPV_max_local = DPV(ih);
            end
            if (DPV(ih) < DPV(ih-delta_i))&(DPV(ih) < DPV(ih+delta_i))
               marque_ihp = 0;
               premier_suivant_DPV_max_local = dernier_DPV_max_local;
               ihp = ih + 1;
               while (marque_ihp == 0)&(ihp < XC(j,1)-X0(1,1))
                   if (DPV(ihp) > DPV(ihp-delta_i))&(DPV(ihp) > DPV(ihp+delta_i))
                      premier_suivant_DPV_max_local = DPV(ihp);
                      marque_ihp = 1;
                   end
                   ihp = ihp + 1;
               end
                   
%               if ((DPV(ih) / dernier_DPV_max_local) < Rap_lim1)
               if ((DPV(ih) / dernier_DPV_max_local) < Rap_lim1) & ( (ih+X0(1,1)-1 > X0(j,1))&(ih+X0(1,1)-1 < X1(j,1))&(ih+X0(1,1)-1 ~= XC(j,1)) )
                  panier_X5 = [panier_X5 ; ih+X0(1,1)  XX(ih) YY(ih) 5];
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
                  Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; ih+X0(1,1)  XX(ih) YY(ih) 5];
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

               end
            end
         end



         if (DPV(ih) <= 0)
            if (DPV(ih) <= DPV(ih-delta_i))&(DPV(ih) <= DPV(ih+delta_i))
               dernier_DPV_min_local = DPV(ih);
            end
            if (DPV(ih) >= DPV(ih-delta_i))&(DPV(ih) >= DPV(ih+delta_i))
%               if (DPV(ih) / dernier_DPV_min_local) < Rap_lim2
               if ( (DPV(ih) / dernier_DPV_min_local) < Rap_lim2 ) & ( (ih+X0(1,1)-1 > X0(j,1))&(ih+X0(1,1)-1 < X1(j,1))&(ih+X0(1,1)-1 ~= XC(j,1)) )
                  panier_X5 = [panier_X5 ; ih+X0(1,1)  XX(ih) YY(ih) 5];
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
                  Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; ih+X0(1,1)  XX(ih) YY(ih) 5];
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

               end
            end
         end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %HB% Detection des points de dissymétrie  ,   codés [6]%
%          if [(DPV(ih) >= 0)&(DPV(ih)> DPV(ih-delta_i))&(DPV(ih) > DPV(ih+delta_i))]
%             DPV_max_local = DPV(ih);
% %             plot(XX(ih),YY(ih),'*G'); %H%
% 
%             if DPV_max_local > DPV_MMAX_LOCAL
%                DPV_MMAX_LOCAL = DPV_max_local;
%                i_DPV_MMAX_LOCAL = ih;
%             end
%          elseif [(DPV(ih) <= 0)&(DPV(ih)<= DPV(ih-delta_i))&(DPV(ih) <= DPV(i+delta_i))]|[(ih ==  l_DPV - 1 )&(DPV(ih) <= DPV(ih-delta_i))]
%             DPV_min_local = DPV(ih);
%             if DPV_min_local < DPV_MMIN_LOCAL
% 
% 
%                DPV_MMIN_LOCAL = DPV_min_local;
%                i_DPV_MMIN_LOCAL = ih;
%                if (ih ==  l_DPV - 1 )&(DPV(l_DPV) <= DPV(l_DPV - 1))
%                   marque_fin_DPV = 1;
%                   i_DPV_MMIN_LOCAL = l_DPV;
%                end
%             end             
%          end

       end
     end

%      if (DPV_MMAX_LOCAL ~= 0)  %&(ih+delta_i2+2 < l_DPV)
%         RAP_symetrie = (((i_DPV_MMAX_LOCAL - (X0(j,1)-X0(1,1)+1))*pas)/(((XC(j,1)-X0(1,1)+1) - i_DPV_MMAX_LOCAL)*pas));
%         RAP_DPV1 = DPV(i_DPV_MMAX_LOCAL) / min(DPV(round((i_DPV_MMAX_LOCAL + (XC(j,1)-X0(1,1)+1))/2)) , DPV(i_DPV_MMAX_LOCAL + delta_i2));
% %        if ( (RAP_symetrie < 0.4)&(RAP_DPV1 > 3) )
%         if ( (RAP_symetrie < 0.4)&(RAP_DPV1 > 3) ) & ( (ih+X0(1,1)-1 > X0(j,1))&(ih+X0(1,1)-1 < X1(j,1))&(ih+X0(1,1)-1 ~= XC(j,1)) )
%            long_demi_intervalle_impulsion = i_DPV_MMAX_LOCAL - (X0(j,1) - X0(1,1)+1);
%            i_fin_impulsion = i_DPV_MMAX_LOCAL + round((long_demi_intervalle_impulsion));
%            panier_X6 = [i_fin_impulsion+X0(1,1)  XX(i_fin_impulsion) YY(i_fin_impulsion) 6];
% %B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
%            Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; panier_X6];
% %B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

%         end
% 
%      end
%      if (DPV_MMIN_LOCAL ~= 0)%&(i-delta_i2 > 0)
% 
%          if (((X1(j,1)-X0(1,1)+1) - i_DPV_MMIN_LOCAL)*pas) == 0
%             RAP_symetrie = (((i_DPV_MMIN_LOCAL - (XC(j,1)-X0(1,1)+1))*pas))/pas;
%          else
%             RAP_symetrie = (((i_DPV_MMIN_LOCAL - (XC(j,1)-X0(1,1)+1))*pas)/(((X1(j,1)-X0(1,1)+1) - i_DPV_MMIN_LOCAL)*pas));
%          end
%          if (i_DPV_MMIN_LOCAL - delta_i2-1) > 0
%             RAP_DPV2 = DPV(i_DPV_MMIN_LOCAL) / DPV(i_DPV_MMIN_LOCAL - delta_i2-1);
%          else
%             RAP_DPV2 = 1;
%          end
% %         if ( (RAP_symetrie > 10)&(RAP_DPV2 > 6) )
%          if ( (RAP_symetrie > 10)&(RAP_DPV2 > 6) ) & ( (ih+X0(1,1)-1 > X0(j,1))&(ih+X0(1,1)-1 < X1(j,1))&(ih+X0(1,1)-1 ~= XC(j,1)) )
%             if (((X1(j,1)-X0(1,1)+1) - i_DPV_MMIN_LOCAL)*pas) == 0
%                long_demi_intervalle_impulsion = 3;
%             else
%                long_demi_intervalle_impulsion = (X1(j,1)-X0(1,1)+1) - i_DPV_MMIN_LOCAL;
%             end
%             i_debut_impulsion = i_DPV_MMIN_LOCAL - round((long_demi_intervalle_impulsion));
%             panier_X6 = [i_debut_impulsion+X0(1,1)  XX(i_debut_impulsion) YY(i_debut_impulsion) 6];
% %B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%
%             Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; panier_X6];
% %B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

%          end
%      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
   j = j + 1;
 end
Suite_Extremums_Ordonnee = [Suite_Extremums_Ordonnee; X1(size(X1,1),:)];
        %HB% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %HB%
  

%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

%%%%%%%%%% ellimination des extremums confondus
L_SEO = size( Suite_Extremums_Ordonnee, 1);
indice_extremum_a_ellimier = [];
Suite_Extremums_Ordonnee_reduite = [Suite_Extremums_Ordonnee(1,:)];
for j = 2 : L_SEO
    indice_t_extremum_j_moins_1 = Suite_Extremums_Ordonnee(j-1,1);
    indice_t_extremum_j = Suite_Extremums_Ordonnee(j,1);
    t_extremum_j_moins_1 = Suite_Extremums_Ordonnee(j-1,2);
    t_extremum_j = Suite_Extremums_Ordonnee(j,2);
    if ( t_extremum_j == t_extremum_j_moins_1 )|( indice_t_extremum_j == indice_t_extremum_j_moins_1 )
       indice_extremum_a_ellimier = [indice_extremum_a_ellimier; j];
    else
       Suite_Extremums_Ordonnee_reduite = [Suite_Extremums_Ordonnee_reduite; Suite_Extremums_Ordonnee(j,:)];
    end
    
end
Suite_Extremums_Ordonnee = Suite_Extremums_Ordonnee_reduite;

%%%%%%%%%% Définition des instants t0, Tc, et t1 pour chaque intervalle Bêta

L_SEO = size( Suite_Extremums_Ordonnee, 1);

for j = 2 : L_SEO - 1
    debut_intervalle = Suite_Extremums_Ordonnee(j-1,:);
    maximum_intervalle = Suite_Extremums_Ordonnee(j,:);
    fin_intervalle = Suite_Extremums_Ordonnee(j+1,:);
    
    Suite_Intervalles_Beta(j-1,1,:) = debut_intervalle;
    Suite_Intervalles_Beta(j-1,2,:) = maximum_intervalle;
    Suite_Intervalles_Beta(j-1,3,:) = fin_intervalle;
end
    
%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%B%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Définir les paramettres de l'impulsion beta declanchee %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                 dans chaque intervalle                 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% pas = 0.01;

%=================  1er Méthode : estimation de deux impulsions beta une ajutée à gauche  =================%
%=================       et l'autre ajustée à droite. On retient la moyenne des deux      =================%
%%%%%%%%%% initialisation des paramètres

 K_hauteur_preced = 0.00000000001;  % prèsque zéro
 q_preced = 2;
 p_preced = 2;
 t0_preced = 1;
 tc_preced = 10;
 t1_preced = 20;

 t_ajustement1 = Suite_Intervalles_Beta(1,1,2) + (1 * ( Suite_Intervalles_Beta(1,2,2) - Suite_Intervalles_Beta(1,1,2) ) / 2 );%2;
 indice_t_ajustement_sup1 = find(T >= t_ajustement1);
 ind_t1 = 1;
 t_ajustement1 = XX(indice_t_ajustement_sup1(ind_t1));
 Vcurv_t_ajustement1_1 = YY(indice_t_ajustement_sup1(ind_t1));

 t_ajustement2 = Suite_Intervalles_Beta(1,1,2) + (1 * ( Suite_Intervalles_Beta(1,2,2) - Suite_Intervalles_Beta(1,1,2) ) / 2 );%2;
 indice_t_ajustement_sup2 = find(T >= t_ajustement2);
 ind_t2 = 1;
 t_ajustement2 = XX(indice_t_ajustement_sup2(ind_t2));
 Vcurv_t_ajustement1_2 = YY(indice_t_ajustement_sup2(ind_t2));

 
 fonctions_beta_ajustee_a_gauche = [];
 fonctions_beta_ajustee_a_droite = [];
 memoire_points_passage_beta_methode1 = [];
 
 nombre_total_impulsions_beta = size(Suite_Intervalles_Beta,1);
 n_hb = 1;%1%2;
 fact_min = 1;%1.07;

 
 for j = 1 : nombre_total_impulsions_beta,

%%%%%%%%%%   identification des parametres t0 , tc, t1, K_hauteur       
     t0 = Suite_Intervalles_Beta(j,1,2);
     tc = Suite_Intervalles_Beta(j,2,2);
     t1 = Suite_Intervalles_Beta(j,3,2);
     K_hauteur_actuel = Suite_Intervalles_Beta(j,2,3);
     
%%%%%%%%%%   identification des parametres p et q
     K_hauteur_max = max( K_hauteur_actuel, K_hauteur_preced );
     K_hauteur_min = min( K_hauteur_actuel, K_hauteur_preced );
     a_hb_1 = ( Vcurv_t_ajustement1_1 ^ n_hb ) / ( ((K_hauteur_min * fact_min) ^ n_hb) + (K_hauteur_max ^ n_hb) );
     tpc_beta_actuel_1 = a_hb_1 * ( (K_hauteur_actuel / Vcurv_t_ajustement1_1)^(n_hb - 1) );
     tpc_beta_preced_1 = a_hb_1 * ( (K_hauteur_preced / Vcurv_t_ajustement1_1)^(n_hb - 1) );
     valeur_beta_preced_t_ajustement1_1 = K_hauteur_preced * tpc_beta_preced_1;
     valeur_beta_actuel_t_ajustement1_1 = K_hauteur_actuel * tpc_beta_actuel_1;

     K_hauteur_max = max( K_hauteur_actuel, K_hauteur_preced );
     K_hauteur_min = min( K_hauteur_actuel, K_hauteur_preced );
     a_hb_2 = ( Vcurv_t_ajustement1_2 ^ n_hb ) / ( ((K_hauteur_min * fact_min) ^ n_hb) + (K_hauteur_max ^ n_hb) );
     tpc_beta_actuel_2 = a_hb_2 * ( (K_hauteur_actuel / Vcurv_t_ajustement1_2)^(n_hb - 1) );
     tpc_beta_preced_2 = a_hb_2 * ( (K_hauteur_preced / Vcurv_t_ajustement1_2)^(n_hb - 1) );
     valeur_beta_preced_t_ajustement1_2 = K_hauteur_preced * tpc_beta_preced_2;
     valeur_beta_actuel_t_ajustement1_2 = K_hauteur_actuel * tpc_beta_actuel_2;

     p_1 = abs((log(valeur_beta_actuel_t_ajustement1_1/K_hauteur_actuel))/((log((t_ajustement1 - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement1)/(t1 - tc)))));
     p_2 = abs((log(valeur_beta_actuel_t_ajustement1_2/K_hauteur_actuel))/((log((t_ajustement2 - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement2)/(t1 - tc)))));
     p = (p_1 + p_2) / 2;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%% observation des valeurs limites de p et q %%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
     if (p > 7)                                    % %
        p = 7;                                     % %
     elseif (p < 1)                                % %
        p = 1;                                     % %
     end                                           % %
     q = (p * (((t1 - tc)/(tc - t0))));            % %
     if (q > 8)                                    % %
        q = 8;                                     % %
        p = (q * (((tc - t0)/(t1 - tc))));         % %
     elseif (q < 0.8)                              % %
        q = 0.8;                                   % %
        p = (q * (((tc - t0)/(t1 - tc))));         % %
     end                                           % %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     p_preced_droite_1 = abs((log(valeur_beta_preced_t_ajustement1_1/K_hauteur_preced))/((log((t_ajustement1 - t0_preced)/(tc_preced - t0_preced)))+(((t1_preced - tc_preced)/(tc_preced - t0_preced))*log((t1_preced - t_ajustement1)/(t1_preced - tc_preced)))));
     p_preced_droite_2 = abs((log(valeur_beta_preced_t_ajustement1_2/K_hauteur_preced))/((log((t_ajustement2 - t0_preced)/(tc_preced - t0_preced)))+(((t1_preced - tc_preced)/(tc_preced - t0_preced))*log((t1_preced - t_ajustement2)/(t1_preced - tc_preced)))));
     p_preced_droite = (p_preced_droite_1 + p_preced_droite_2) / 2;
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%% observation des valeurs limites de p et q %%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
     if (p_preced_droite > 7)                      % %
        p_preced_droite = 7;                       % %
     elseif (p_preced_droite < 1)                  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        p_preced_droite = 1;                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
     end                                                                                           % %
     q_preced_droite = (p_preced_droite * (((t1_preced - tc_preced)/(tc_preced - t0_preced))));    % %
     if (q_preced_droite > 8)                                                                      % %
        q_preced_droite = 8;                                                                       % %
        p_preced_droite = (q_preced_droite * (((tc_preced - t0_preced)/(t1_preced - tc_preced)))); % %
     elseif (q_preced_droite < 0.8)                                                                % %
        q_preced_droite = 0.8;                                                                     % %
        p_preced_droite = (q_preced_droite * (((tc_preced - t0_preced)/(t1_preced - tc_preced)))); % %
     end                                                                                           % %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     
%%%%%%%%%%   enregistrement des fonctions beta et de leurs parametres            
     rong_extremum = j + 1;

     if j < nombre_total_impulsions_beta + 1
        fonctions_beta_ajustee_a_gauche = [fonctions_beta_ajustee_a_gauche; t0, tc, t1, p, q, K_hauteur_actuel, rong_extremum ];
     end
     if j > 1
        fonctions_beta_ajustee_a_droite = [fonctions_beta_ajustee_a_droite; t0_preced, tc_preced, t1_preced, p_preced_droite, q_preced_droite, K_hauteur_preced, rong_extremum_preced ];
     end
          
     memoire_points_passage_beta_methode1 = [memoire_points_passage_beta_methode1; t_ajustement1, Vcurv_t_ajustement1_1, valeur_beta_preced_t_ajustement1_1, valeur_beta_actuel_t_ajustement1_1];

%%%%%%%%%%   definition de l'instant t_ajustement suivant
     t_ajustement1 = ( tc + (1*( t1 - tc )/ 2) ) ;
     indice_t_ajustement_sup1 = find(XX >= t_ajustement1);
     indice_t_ajustement_inf1 = find(XX <= t_ajustement1);
     ind_t1 = 1;
     t_ajustement1_s = XX(indice_t_ajustement_sup1(ind_t1));
     t_ajustement1_i = XX(indice_t_ajustement_inf1(ind_t1));
     Vcurv_t_ajustement1_s_1 = YY(indice_t_ajustement_sup1(ind_t1));
     Vcurv_t_ajustement1_i_1 = YY(indice_t_ajustement_inf1(ind_t1));
     if (t_ajustement1_s - t_ajustement1_i) ~= 0
        Vcurv_t_ajustement1_1 = Vcurv_t_ajustement1_i_1 + ( (Vcurv_t_ajustement1_s_1  -  Vcurv_t_ajustement1_i_1) * ((t_ajustement1 - t_ajustement1_i) / (t_ajustement1_s - t_ajustement1_i)) );
     else
        Vcurv_t_ajustement1_1 = Vcurv_t_ajustement1_i_1;
     end

%      if (t1 - tc == 0.01)
%         Vcurv_t_ajustement1_1 = ( K_hauteur_actuel + Suite_Intervalles_Beta(j+1,2,3) )/2 ;
%      end

     t_ajustement2 = ( tc + (1*( t1 - tc )/ 2) ) ;
     indice_t_ajustement_sup2 = find(XX >= t_ajustement2);
     indice_t_ajustement_inf2 = find(XX <= t_ajustement2);
     ind_t2 = 1;
     t_ajustement2_s = XX(indice_t_ajustement_sup2(ind_t2));
     t_ajustement2_i = XX(indice_t_ajustement_inf2(ind_t2));
     Vcurv_t_ajustement1_s_2 = YY(indice_t_ajustement_sup2(ind_t2));
     Vcurv_t_ajustement1_i_2 = YY(indice_t_ajustement_inf2(ind_t2));
     if (t_ajustement2_s - t_ajustement2_i) ~= 0
        Vcurv_t_ajustement1_2 = Vcurv_t_ajustement1_i_2 + ( (Vcurv_t_ajustement1_s_2  -  Vcurv_t_ajustement1_i_2) * ((t_ajustement2 - t_ajustement2_i) / (t_ajustement2_s - t_ajustement2_i)) );
     else
        Vcurv_t_ajustement1_2 = Vcurv_t_ajustement1_i_2;
     end

%      if (t1 - tc == 0.01)&(j+1 <= nombre_total_impulsions_beta )
%         Vcurv_t_ajustement1_1 = ( K_hauteur_actuel + Suite_Intervalles_Beta(j+1,2,3) )/2 ;
%         Vcurv_t_ajustement1_2 = ( K_hauteur_actuel + Suite_Intervalles_Beta(j+1,2,3) )/2 ;
%      end


%%%%%%%%%% mémorisation des paramètres du béta précédente

     K_hauteur_preced = K_hauteur_actuel;
     q_preced = q;
     p_preced = p;
     t0_preced = t0;
     tc_preced = tc;
     t1_preced = t1;     
     rong_extremum_preced = rong_extremum;
     
 end     
    
%%%%%%%%%% definition de la derniere fonctions_beta_ajustee_a_droite %%%%%%
 t0 = Suite_Intervalles_Beta(nombre_total_impulsions_beta,2,2);
 tc = Suite_Intervalles_Beta(nombre_total_impulsions_beta,3,2);
 t1 = tc + (tc - t0);
 K_hauteur_actuel = 0.0000000001; % prèsque zéro

%%%%%%%%%%   identification des parametres p et q
 K_hauteur_max = max( K_hauteur_actuel, K_hauteur_preced );
 K_hauteur_min = min( K_hauteur_actuel, K_hauteur_preced );
 a_hb_1 = ( Vcurv_t_ajustement1_1 ^ n_hb ) / ( ((K_hauteur_min * fact_min) ^ n_hb) + (K_hauteur_max ^ n_hb) );
 tpc_beta_actuel_1 = a_hb_1 * ( (K_hauteur_actuel / Vcurv_t_ajustement1_1)^(n_hb - 1) );
 tpc_beta_preced_1 = a_hb_1 * ( (K_hauteur_preced / Vcurv_t_ajustement1_1)^(n_hb - 1) );
 valeur_beta_preced_t_ajustement1_1 = K_hauteur_preced * tpc_beta_preced_1;
 valeur_beta_actuel_t_ajustement1_1 = K_hauteur_actuel * tpc_beta_actuel_1;

 K_hauteur_max = max( K_hauteur_actuel, K_hauteur_preced );
 K_hauteur_min = min( K_hauteur_actuel, K_hauteur_preced );
 a_hb_2 = ( Vcurv_t_ajustement1_2 ^ n_hb ) / ( ((K_hauteur_min * fact_min) ^ n_hb) + (K_hauteur_max ^ n_hb) );
 tpc_beta_actuel_2 = a_hb_2 * ( (K_hauteur_actuel / Vcurv_t_ajustement1_2)^(n_hb - 1) );
 tpc_beta_preced_2 = a_hb_2 * ( (K_hauteur_preced / Vcurv_t_ajustement1_2)^(n_hb - 1) );
 valeur_beta_preced_t_ajustement1_2 = K_hauteur_preced * tpc_beta_preced_2;
 valeur_beta_actuel_t_ajustement1_2 = K_hauteur_actuel * tpc_beta_actuel_2;

 p_1 = abs((log(valeur_beta_actuel_t_ajustement1_1/K_hauteur_actuel))/((log((t_ajustement1 - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement1)/(t1 - tc)))));
 p_2 = abs((log(valeur_beta_actuel_t_ajustement1_2/K_hauteur_actuel))/((log((t_ajustement2 - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement2)/(t1 - tc)))));
 p = (p_1 + p_2) / 2;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% observation des valeurs limites de p et q %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 if (p > 7)                                    % %
    p = 7;                                     % %
 elseif (p < 1)                                % %
    p = 1;                                     % %
 end                                           % %
 q = (p * (((t1 - tc)/(tc - t0))));            % %
 if (q > 8)                                    % %
    q = 8;                                     % %
    p = (q * (((tc - t0)/(t1 - tc))));         % %
 elseif (q < 0.8)                              % %
    q = 0.8;                                   % %
    p = (q * (((tc - t0)/(t1 - tc))));         % %
 end                                           % %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 p_preced_droite_1 = abs((log(valeur_beta_preced_t_ajustement1_1/K_hauteur_preced))/((log((t_ajustement1 - t0_preced)/(tc_preced - t0_preced)))+(((t1_preced - tc_preced)/(tc_preced - t0_preced))*log((t1_preced - t_ajustement1)/(t1_preced - tc_preced)))));
 p_preced_droite_2 = abs((log(valeur_beta_preced_t_ajustement1_2/K_hauteur_preced))/((log((t_ajustement2 - t0_preced)/(tc_preced - t0_preced)))+(((t1_preced - tc_preced)/(tc_preced - t0_preced))*log((t1_preced - t_ajustement2)/(t1_preced - tc_preced)))));
 p_preced_droite = (p_preced_droite_1 + p_preced_droite_2) / 2;
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%% observation des valeurs limites de p et q %%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 if (p_preced_droite > 7)                      % %
    p_preced_droite = 7;                       % %
 elseif (p_preced_droite < 1)                  % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p_preced_droite = 1;                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 end                                                                                           % %
 q_preced_droite = (p_preced_droite * (((t1_preced - tc_preced)/(tc_preced - t0_preced))));    % %
 if (q_preced_droite > 8)                                                                      % %
    q_preced_droite = 8;                                                                       % %
    p_preced_droite = (q_preced_droite * (((tc_preced - t0_preced)/(t1_preced - tc_preced)))); % %
 elseif (q_preced_droite < 0.8)                                                                % %
    q_preced_droite = 0.8;                                                                     % %
    p_preced_droite = (q_preced_droite * (((tc_preced - t0_preced)/(t1_preced - tc_preced)))); % %
 end                                                                                           % %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%   enregistrement des fonctions beta et de leurs parametres            
 if j < nombre_total_impulsions_beta
    fonctions_beta_ajustee_a_gauche = [fonctions_beta_ajustee_a_gauche; t0, tc, t1, p, q, K_hauteur_actuel, rong_extremum ];
 end
 if (j > 1)|(nombre_total_impulsions_beta == 1)
    fonctions_beta_ajustee_a_droite = [fonctions_beta_ajustee_a_droite; t0_preced, tc_preced, t1_preced, p_preced_droite, q_preced_droite, K_hauteur_preced, rong_extremum_preced ];
 end

memoire_points_passage_beta_methode1 = [memoire_points_passage_beta_methode1; t_ajustement1, Vcurv_t_ajustement1_1, valeur_beta_preced_t_ajustement1_1, valeur_beta_actuel_t_ajustement1_1];

% fonctions_beta_ajustee_moyenne = ( (1 * fonctions_beta_ajustee_a_gauche) + (1 * fonctions_beta_ajustee_a_droite) ) / 2;
fonctions_beta_ajustee_moyenne = fonctions_beta_ajustee_a_gauche;
fonctions_beta_ajustee_moyenne_k = [fonctions_beta_ajustee_a_gauche fonctions_beta_ajustee_a_droite];

[fonctions_beta_ajustee_moyenne_m] = algorithme_iteratif_estimation_fcts_beta_ajustee_moyenne( fonctions_beta_ajustee_a_gauche , fonctions_beta_ajustee_a_droite );
%fonctions_beta_ajustee_moyenne_m = ( (1 * fonctions_beta_ajustee_a_gauche) + (1 * fonctions_beta_ajustee_a_droite) ) / 2; 


%============        2ème Méthode : estimation des paramètres des impulsions beta par ajustement        ============%
%============ successif (par intervalle) de la vitesse reconstruite par rapport à la vitesse curviligne ============%

 nombre_total_impulsions_beta = size(Suite_Intervalles_Beta,1);
 Suite_Intervalles_Beta_corrige = Suite_Intervalles_Beta;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (facultatif) %%%%
%%%%%%%%%% correction des estimations excentriques de l'instant tc %%%%
 K_rap_min = 0.6;%0.4;
 K_rap_max = 1.8;%2.5;
 for j = 1 +1 : nombre_total_impulsions_beta -1,
     t0 = Suite_Intervalles_Beta_corrige(j,1,2);
     tc = Suite_Intervalles_Beta_corrige(j,2,2);
     t1 = Suite_Intervalles_Beta_corrige(j,3,2);
     K_rap_centricite = (tc - t0) / (t1 - tc);
     if  ( K_rap_centricite > K_rap_min ) & ( K_rap_centricite < K_rap_max )
         Suite_Intervalles_Beta_corrige(j,2,2) = tc;
         Suite_Intervalles_Beta_corrige(j+1,1,2) = tc;
         Suite_Intervalles_Beta_corrige(j-1,3,2) = tc;
     elseif ( K_rap_centricite > K_rap_max )
         K_rap_centricite = K_rap_max;
         tc_corrige = round( (( t0 + (K_rap_centricite * t1) ) / (1 + K_rap_centricite)) / pas ) * pas;
         Suite_Intervalles_Beta_corrige(j,2,2) = tc_corrige;
         Suite_Intervalles_Beta_corrige(j+1,1,2) = tc_corrige;
         Suite_Intervalles_Beta_corrige(j-1,3,2) = tc_corrige;
     else %( K_rap_centricite < K_rap_min )
         K_rap_centricite = K_rap_min;
         tc_corrige = round( (( t0 + (K_rap_centricite * t1) ) / (1 + K_rap_centricite)) / pas ) * pas;
         Suite_Intervalles_Beta_corrige(j,2,2) = tc_corrige;
         Suite_Intervalles_Beta_corrige(j+1,1,2) = tc_corrige;
         Suite_Intervalles_Beta_corrige(j-1,3,2) = tc_corrige;
     end
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% initialisation des paramètres
 
 %%% pas = 0.01;
 
 K_hauteur_preced = 50;     % facultatif
 q_preced = 20;
 p_preced = 20;
 t0_preced = 1;
 tc_preced = 10;
 t1_preced = 20;

 beta_preced_t_ajustement1 = 0;
 t_ajustement = Suite_Intervalles_Beta_corrige(1,1,2) + (1 * ( Suite_Intervalles_Beta_corrige(1,2,2) - Suite_Intervalles_Beta_corrige(1,1,2) ) / 2 );%2;
 indice_t_ajustement_sup = find(XX >= t_ajustement);
 indice_t_ajustement_inf = find(XX <= t_ajustement);
 ind_t = 1;
 t_ajustement_s = XX(indice_t_ajustement_sup(ind_t));
 t_ajustement_i = XX(indice_t_ajustement_inf(ind_t));
 Vcurv_t_ajustement1_s = YY(indice_t_ajustement_sup(ind_t));
 Vcurv_t_ajustement1_i = YY(indice_t_ajustement_inf(ind_t));
 if (t_ajustement_s - t_ajustement_i) ~= 0
    Vcurv_t_ajustement1 = Vcurv_t_ajustement1_i + ( (Vcurv_t_ajustement1_s  -  Vcurv_t_ajustement1_i) * ((t_ajustement - t_ajustement_i) / (t_ajustement_s - t_ajustement_i)) );
 else
    Vcurv_t_ajustement1 = Vcurv_t_ajustement1_i;
 end
   
 fonctions_beta = [];
 fonctions_beta_p_constante = [];
 memoire_points_passage_beta_methode2 = [];
 nombre_total_impulsions_beta = size(Suite_Intervalles_Beta_corrige,1);

 for j = 1 : nombre_total_impulsions_beta,

%%%%%%%%%%  identification des parametres t0 , tc, t1, K_hauteur       
     t0 = Suite_Intervalles_Beta_corrige(j,1,2);
     tc = Suite_Intervalles_Beta_corrige(j,2,2);
     t1 = Suite_Intervalles_Beta_corrige(j,3,2);
     K_hauteur = Suite_Intervalles_Beta_corrige(j,2,3);

%%%%%%%%%%   identification des parametres p et q
     beta_actuel_t_ajustement1 = (Vcurv_t_ajustement1) - (beta_preced_t_ajustement1);
     if beta_actuel_t_ajustement1 > 0
        beta_actuel_t_ajustement1 = beta_actuel_t_ajustement1;
     else
        beta_actuel_t_ajustement1 = K_hauteur / 1000;
     end

     if beta_actuel_t_ajustement1 > K_hauteur
        p = 2.5;
        p_v1 = p;
        q=(p*(((t1 - tc)/(tc - t0))));
        q_v1 = q;
     else

%>>>>>>>>> estimation de p et q par ajustement de la vitesse reconstruite <<<<<<<<<<<%
%>>>>>>>>> sur la vitesse curviligne en instant intermidaire t_ajustement <<<<<<<<<<<%

        p = abs((log(beta_actuel_t_ajustement1/K_hauteur))/((log((t_ajustement - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement)/(t1 - tc)))));
        p_v1 = p;
        q=(p*(((t1 - tc)/(tc - t0))));
        q_v1 = q;
     end
     
%%%%%%%%%%   mémorisation de l'instant t_ajustement actuel
     t_ajustement_actuel = t_ajustement;
     Vcurv_t_ajustement1_actuel = Vcurv_t_ajustement1;
%%%%%%%%%%   definition de l'instant t_ajustement suivant
     t_ajustement = ( tc + (1*( t1 - tc )/ 2) ) ;%2;
     indice_t_ajustement_sup = find(XX >= t_ajustement);    %H%
     indice_t_ajustement_inf = find(XX <= t_ajustement);    %H%
     ind_t = 1;
     t_ajustement_s = XX(indice_t_ajustement_sup(ind_t));
     t_ajustement_i = XX(indice_t_ajustement_inf(ind_t));
     Vcurv_t_ajustement1_s = YY(indice_t_ajustement_sup(ind_t));
     Vcurv_t_ajustement1_i = YY(indice_t_ajustement_inf(ind_t));
     if (t_ajustement_s - t_ajustement_i) ~= 0
        Vcurv_t_ajustement1 = Vcurv_t_ajustement1_i + ( (Vcurv_t_ajustement1_s  -  Vcurv_t_ajustement1_i) * ((t_ajustement - t_ajustement_i) / (t_ajustement_s - t_ajustement_i)) );
     else
        Vcurv_t_ajustement1 = Vcurv_t_ajustement1_i;
     end
     %%% delta_temps = 0.01;%pas;
     delta_temps = pas;
     Derivee_Vsig_t_ajustement = ( (YY(indice_t_ajustement_sup(1))) - (YY((indice_t_ajustement_sup(1))-1)) ) / ( delta_temps );
%,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,%

%.............................................%
%%%%%%%%%% Calcul d'une autre estimation de p et q de l'impulsion beta_actuelle en tenant compte
%%%%%%%%%% d'une estimée de l'impulsion Beta suivante

     if j < nombre_total_impulsions_beta
        t0_suiv = Suite_Intervalles_Beta_corrige(j+1,1,2);
        tc_suiv = Suite_Intervalles_Beta_corrige(j+1,2,2);
        t1_suiv = Suite_Intervalles_Beta_corrige(j+1,3,2);
        K_hauteur_suiv = Suite_Intervalles_Beta_corrige(j+1,2,3);

        %%%% estimation des paramètres p et q de l'impulsion suivante
        if  (tc_suiv - t0_suiv) > (t1_suiv - tc_suiv)
            q_suiv = 2;
            p_suiv = (q_suiv*(((tc_suiv - t0_suiv)/(t1_suiv - tc_suiv))));
            if p_suiv > 4
               p_suiv = 4;
               q_suiv =( p_suiv * (((t1_suiv - tc_suiv)/(tc_suiv - t0_suiv))) );
            end
        else
            p_suiv = 2;
            q_suiv =( p_suiv * (((t1_suiv - tc_suiv)/(tc_suiv - t0_suiv))) );
            if q_suiv > 4
               q_suiv = 4;
               p_suiv = (q_suiv*(((tc_suiv - t0_suiv)/(t1_suiv - tc_suiv))));
            end
        end
        %%%% Calcul de la valeur estimee beta_suivante_t_ajustement1
        beta_suivante_estimee_t_ajustement1 = K_hauteur_suiv * ((((t_ajustement - t0_suiv)/(tc_suiv - t0_suiv))^p_suiv)*(((t1_suiv - t_ajustement)/(t1_suiv - tc_suiv))^q_suiv));
        
        %%%% Calcul d'une autre estimation de beta_act_ajustement1
        beta_actuelle_estimee_pour_nouveau_t_ajustement1 = Vcurv_t_ajustement1 - beta_suivante_estimee_t_ajustement1;
        p_v2 = abs((log(beta_actuelle_estimee_pour_nouveau_t_ajustement1/K_hauteur))/((log((t_ajustement - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement)/(t1 - tc)))));
        q_v2 = (p_v2*(((t1 - tc)/(tc - t0))));
        
        %%%% retien de la moyenne des deux estimations faites en considérant
        %%%% l'impulsion précédente et une estimée de l'impulsion suivante
        p = ( (1 * p_v1) + (1 * p_v2) ) / 2;
        q = ( (1 * q_v1) + (1 * q_v2) ) / 2;

     else
        p = p;
        q=(p*(((t1 - tc)/(tc - t0))));
     end
%.............................................%

%>>>>>>>>>>>>>> contrôle des resultats obtenus pour p et q <<<<<<<<<<<<<<<%
%+++++++++++++++++++++++++++++++++++++++++++++%
%%%%%%%%%% eviter l'apparition de minimums locaux intermidiaires

     if (j > 1)%&(j < nombre_total_impulsions_beta)
    
        if ( K_hauteur < K_hauteur_preced )
           t_test = tc - ((tc - t0)/8);
        else
           t_test = t0 + ((tc - t0)/8);
        end
     
        rapport_t_test1 = (K_hauteur_preced)*(q_preced/(t1_preced - tc_preced))*(((t1_preced-t0_preced)/(tc_preced-t0_preced))^p_preced)*(((t1_preced-t_test)/(t1_preced-tc_preced))^(q_preced-1))*( 1 / ((K_hauteur)*(((t1-t_test)-((t_test-t0)*(t1-tc)/(tc-t0)))/((tc-t0)*(t1-tc)))) );
        rapport_t_test2 = (K_hauteur)*(p/(tc - t0))*(((t1-t0)/(t1-tc))^q)*(((t_test-t0)/(tc-t0))^(p-1)) * ( 1 / ((K_hauteur_preced)*(((t_test-t0_preced)-((t1_preced-t_test)*(tc_preced-t0_preced)/(t1_preced-tc_preced)))/((t1_preced-tc_preced)*(tc-t0_preced)))) );

        if (p > rapport_t_test1) & ( K_hauteur < K_hauteur_preced )
           p_v2 = rapport_t_test1 * ( K_hauteur_preced / K_hauteur );
           p = ((p_v2 * 3) + (p * 2)) / 5;
           q=(p*(((t1 - tc)/(tc - t0))));
              

        elseif ( q_preced > rapport_t_test2 )& ( K_hauteur > K_hauteur_preced )
           q_preced_v2 = rapport_t_test2 * ( K_hauteur_preced / K_hauteur );
           q_preced = ((q_preced_v2 * 3) + (q_preced * 2)) / 5;
%           q_preced = ((q_preced_v2 * 5) + (q_preced * 0)) / 5;
%           q_preced = ((q_preced_v2 * 1) + (q_preced * 1)) / 2;
        
           p_preced = (q_preced*(((tc_preced - t0_preced)/(t1_preced - tc_preced))));
           fonctions_beta = [fonctions_beta(1:j-2,:); t0_preced, tc_preced, t1_preced, p_preced, q_preced, K_hauteur_preced, rong_extremum_preced];
           beta_preced_t_ajustement1 = K_hauteur_preced*(((((t_ajustement_actuel - t0_preced)/(tc_preced - t0_preced))^p_preced)*(((t1_preced - t_ajustement_actuel)/(t1_preced - tc_preced))^q_preced)));
           beta_actuel_t_ajustement1 = (Vcurv_t_ajustement1_actuel) - (beta_preced_t_ajustement1);
           if beta_actuel_t_ajustement1 > 0
              beta_actuel_t_ajustement1 = beta_actuel_t_ajustement1;
           else
              beta_actuel_t_ajustement1 = K_hauteur / 1000;
           end
           if beta_actuel_t_ajustement1 > K_hauteur
              p = p;
              q = q;
           else
              p = abs((log(beta_actuel_t_ajustement1/K_hauteur))/((log((t_ajustement_actuel - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement_actuel)/(t1 - tc)))));
              q=(p*(((t1 - tc)/(tc - t0))));
           end
        
        end
   
     end

%+++++++++++++++++++++++++++++++++++++++++++++%

%---------------------------------------------%
%%%%%%%%%% éviter l'apparition de maximums locaux intermidiaires
     if (j > 1)&(j < nombre_total_impulsions_beta)
        K_hauteur_max = max( K_hauteur , K_hauteur_preced );
        K_hauteur_min = min( K_hauteur , K_hauteur_preced );
        if ( p < 2 + (K_hauteur_min / K_hauteur_max) ) %& ( K_hauteur < K_hauteur_preced )
           p = p + (( ( 2 + (K_hauteur_min / K_hauteur_max) ) - p ) /2 ) + max(0 , (( ( 2 + (K_hauteur_min / K_hauteur_max) ) - q_preced ) /2.5 ) );%/3 ) );%/4 ) );
           q=(p*(((t1 - tc)/(tc - t0))));
        end
        if ( q_preced < 2 + (K_hauteur_min / K_hauteur_max) ) %& ( K_hauteur > K_hauteur_preced )
           q_preced = q_preced + (( ( 2 + (K_hauteur_min / K_hauteur_max) ) - q_preced ) /2 ) + max(0, (( ( 2 + (K_hauteur_min / K_hauteur_max) ) - p ) /3 ) );%/3 ) );%/4 ) );
           p_preced =(q_preced*(((tc_preced - t0_preced)/(t1_preced - tc_preced))));
           fonctions_beta = [fonctions_beta(1:j-2,:); t0_preced, tc_preced, t1_preced, p_preced, q_preced, K_hauteur_preced, rong_extremum_preced];           
           beta_preced_t_ajustement1 = K_hauteur_preced*(((((t_ajustement_actuel - t0_preced)/(tc_preced - t0_preced))^p_preced)*(((t1_preced - t_ajustement_actuel)/(t1_preced - tc_preced))^q_preced)));
           beta_actuel_t_ajustement1 = (Vcurv_t_ajustement1_actuel) - (beta_preced_t_ajustement1);
           if beta_actuel_t_ajustement1 > 0
              beta_actuel_t_ajustement1 = beta_actuel_t_ajustement1;
           else
              beta_actuel_t_ajustement1 = K_hauteur / 1000;
           end
           p = abs((log(beta_actuel_t_ajustement1/K_hauteur))/((log((t_ajustement_actuel - t0)/(tc - t0)))+(((t1 - tc)/(tc - t0))*log((t1 - t_ajustement_actuel)/(t1 - tc)))));
           q=(p*(((t1 - tc)/(tc - t0))));
        end

     end
%---------------------------------------------%

%+++++++++++++++++++++++++++++++++++++++++++++%
%%%%%%%%%% éviter l'apparition de minimums locaux intermidiaires

     if (j > 1)%&(j < nombre_total_impulsions_beta)
    
        if ( K_hauteur < K_hauteur_preced )
           t_test = tc - ((tc - t0)/8);%8);
        else
           t_test = t0 + ((tc - t0)/8);%8);
        end
        dimH1 = (((t_test-t0)/(tc-t0))^(p-1))*(((t1-t_test)/(t1-tc))^(q-1));
        dimH2 = (((t1_preced-t_test)/(t1_preced-tc_preced))^(q_preced-1))*(((t_test-t0_preced)/(tc_preced-t0_preced))^(p_preced-1));
        rapport_t_test1 = (K_hauteur_preced)*(((q_preced*(t_test-t0_preced))-(p_preced*(t1_preced-t_test)))/((tc_preced-t0_preced)*(t1_preced-tc_preced)))*(((t_test-t0_preced)/(tc_preced-t0_preced))^(p_preced-1))*(((t1_preced-t_test)/(t1_preced-tc_preced))^(q_preced-1)) * ( 1 / ((K_hauteur)*( (((t1-t_test)-((t_test-t0)*(t1-tc)/(tc-t0)))/((tc-t0)*(t1-tc))) * (dimH1) ) ) );
        rapport_t_test2 = (K_hauteur)*(((p*(t1-t_test))-(q*(t_test-t0)))/((t1-tc)*(tc-t0)))*(((t1-t_test)/(t1-tc))^(q-1))*(((t_test-t0)/(tc-t0))^(p-1)) * ( 1 / ((K_hauteur_preced)*( (((t_test-t0_preced)-((t1_preced-t_test)*(tc_preced-t0_preced)/(t1_preced-tc_preced)))/((t1_preced-tc_preced)*(tc_preced-t0_preced))) * (dimH2) ) ) );

        if (p > rapport_t_test1) & ( K_hauteur < K_hauteur_preced )
           p_v2 = rapport_t_test1;
           p = ((p_v2 * 3) + (p * 2)) / 5;
%           p = ((p_v2 *5) + (p * 0)) / 5;
%           p = ((p_v2 *2.5) + (p * 2.5)) / 5;
       
           q=(p*(((t1 - tc)/(tc - t0))));
        elseif ( q_preced > rapport_t_test2 )& ( K_hauteur > K_hauteur_preced )
           q_preced_v2 = rapport_t_test2 ;
%           q_preced = rapport_t_test2 ;
           q_preced = ((q_preced_v2 * 3) + (q_preced * 2)) / 5;
%           q_preced = ((q_preced_v2 * 2.5) + (q_preced * 2.5)) / 5;
%           q_preced = ((q_preced_v2 * 5) + (q_preced * 0)) / 5;
        
           p_preced = (q_preced*(((tc_preced - t0_preced)/(t1_preced - tc_preced))));
           fonctions_beta = [fonctions_beta(1:j-2,:); t0_preced, tc_preced, t1_preced, p_preced, q_preced, K_hauteur_preced, rong_extremum_preced];
        end
   
     end
%+++++++++++++++++++++++++++++++++++++++++++++%

%.............................................%
%%%%%%%%%% filtrage linéaire des paramètres calculés de l'impulsion bêta
     if p < 2
        p_v2 = 2;%1.5;
        p = ((1.2*p_v2) + (3.8*p))/5;
        q=(p*(((t1 - tc)/(tc - t0))));
     elseif q < 2
        q_v2 = 2;%1.5;
        q = ((1.2*q_v2) + (3.8*q))/5;
        p=(q*(((tc - t0)/(t1 - tc))));
     end

     beta_preced_t_ajustement1 = K_hauteur_preced*(((((t_ajustement_actuel - t0_preced)/(tc_preced - t0_preced))^p_preced)*(((t1_preced - t_ajustement_actuel)/(t1_preced - tc_preced))^q_preced)));
     beta_actuel_t_ajustement1 = (Vcurv_t_ajustement1_actuel) - (beta_preced_t_ajustement1);

     memoire_points_passage_beta_methode2 = [memoire_points_passage_beta_methode2; t_ajustement_actuel, Vcurv_t_ajustement1_actuel, beta_preced_t_ajustement1, beta_actuel_t_ajustement1];
%..............................................%
%%%%%%%%%% mémorisation de la valeur de l'impulsion bêta traitée à l'instant t_ajustement1
     beta_preced_t_ajustement1 = K_hauteur*((((t_ajustement - t0)/(tc - t0))^p)*(((t1 - t_ajustement)/(t1 - tc))^q));
%%%%%%%%%% mémorisation des paramètres calculés de l'impulsion bêta traité (-> précédente)
     K_hauteur_preced = K_hauteur;
     q_preced = q;
     p_preced = p;
     t0_preced = t0;
     tc_preced = tc;
     t1_preced = t1;
     
%%%%%%%%%%   enregistrement des fonctions beta et de leurs parametres
     rong_extremum = j + 1;
     rong_extremum_preced = rong_extremum;
     
     fonctions_beta = [fonctions_beta; t0, tc, t1, p, q, K_hauteur, rong_extremum];
     
%%%%%%%%%%  méthode p = constante = 2     
     p_methode_constante = 3.5;
     q_methode_constante=(p_methode_constante*(((t1 - tc)/(tc - t0))));
     fonctions_beta_p_constante = [fonctions_beta_p_constante; t0, tc, t1, p_methode_constante, q_methode_constante, K_hauteur, rong_extremum];

    
 end

 beta_actuel_t_ajustement1 = 0;
 memoire_points_passage_beta_methode2 = [memoire_points_passage_beta_methode2; t_ajustement, Vcurv_t_ajustement1, beta_preced_t_ajustement1, beta_actuel_t_ajustement1];


%============       repartition des fonctions beta chevauchées sur les intervalles de temps        ============%
%==============================================================================================================%

%%%%%%%%%%   pour la première méthode
nbr_fb = size(fonctions_beta_ajustee_moyenne,1);

temps1 = fonctions_beta_ajustee_moyenne(1,1);
temps2 = fonctions_beta_ajustee_moyenne(1,2);
rong_fonction1 = 0;
rong_fonction2 = 1;

K_amplitude1 = 0;
K_amplitude2 = fonctions_beta_ajustee_moyenne(1,6);
indice_data_1 = Suite_Extremums_Ordonnee(1,1);
indice_data_2 = Suite_Extremums_Ordonnee(2,1);

tableau_affectation_methode_1 = [temps1, temps2, rong_fonction1, rong_fonction2, K_amplitude1, K_amplitude2, indice_data_1, indice_data_2];

for j = 2 : nbr_fb
    temps1 = fonctions_beta_ajustee_moyenne(j,1);
    temps2 = fonctions_beta_ajustee_moyenne(j,2);
    rong_fonction1 = j-1;
    rong_fonction2 = j;
    rong_extremum1_liste = fonctions_beta_ajustee_moyenne(rong_fonction1,7);
    rong_extremum2_liste = fonctions_beta_ajustee_moyenne(rong_fonction2,7);
    
    K_amplitude1 = K_amplitude2;
    K_amplitude2 = fonctions_beta_ajustee_moyenne(j,6);
    indice_data_1 = Suite_Extremums_Ordonnee(j,1);
    indice_data_2 = Suite_Extremums_Ordonnee(j+1,1);

    tableau_affectation_methode_1 = [tableau_affectation_methode_1; temps1, temps2, rong_fonction1, rong_fonction2, K_amplitude1, K_amplitude2, indice_data_1, indice_data_2];
end

temps1 = fonctions_beta_ajustee_moyenne(nbr_fb,2);
temps2 = fonctions_beta_ajustee_moyenne(nbr_fb,3);
rong_fonction1 = nbr_fb;
rong_fonction2 = 0;

K_amplitude1 = K_amplitude2;
K_amplitude2 = 0;
indice_data_1 = Suite_Extremums_Ordonnee(nbr_fb+1,1);
indice_data_2 = Suite_Extremums_Ordonnee(nbr_fb+2,1);

tableau_affectation_methode_1 = [tableau_affectation_methode_1; temps1, temps2, rong_fonction1, rong_fonction2, K_amplitude1, K_amplitude2, indice_data_1, indice_data_2];

%%%%%%%%%%   pour la deuxième méthode
nbr_fb = size(fonctions_beta,1);

temps1 = fonctions_beta(1,1);
temps2 = fonctions_beta(1,2);
rong_fonction1 = 0;
rong_fonction2 = 1;

K_amplitude1 = 0;
K_amplitude2 = fonctions_beta(1,6);
indice_data_1 = Suite_Extremums_Ordonnee(1,1);
indice_data_2 = Suite_Extremums_Ordonnee(2,1);

tableau_affectation_methode_2 = [temps1, temps2, rong_fonction1, rong_fonction2, K_amplitude1, K_amplitude2, indice_data_1, indice_data_2];

for j = 2 : nbr_fb
    temps1 = fonctions_beta(j,1);
    temps2 = fonctions_beta(j,2);
    rong_fonction1 = j-1;
    rong_fonction2 = j;
    rong_extremum1_liste = fonctions_beta(rong_fonction1,7);
    rong_extremum2_liste = fonctions_beta(rong_fonction2,7);
    indice_data_1 = Suite_Extremums_Ordonnee(j,1);
    indice_data_2 = Suite_Extremums_Ordonnee(j+1,1);

    K_amplitude1 = K_amplitude2;
    K_amplitude2 = fonctions_beta(j,6);

    tableau_affectation_methode_2 = [tableau_affectation_methode_2; temps1, temps2, rong_fonction1, rong_fonction2, K_amplitude1, K_amplitude2, indice_data_1, indice_data_2];
end

temps1 = fonctions_beta(nbr_fb,2);
temps2 = fonctions_beta(nbr_fb,3);
rong_fonction1 = nbr_fb;
rong_fonction2 = 0;

K_amplitude1 = K_amplitude2;
K_amplitude2 = 0;
indice_data_1 = Suite_Extremums_Ordonnee(nbr_fb+1,1);
indice_data_2 = Suite_Extremums_Ordonnee(nbr_fb+2,1);

tableau_affectation_methode_2 = [tableau_affectation_methode_2; temps1, temps2, rong_fonction1, rong_fonction2, K_amplitude1, K_amplitude2, indice_data_1, indice_data_2];






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fonctions_beta_sortie = fonctions_beta_ajustee_moyenne_m;
%fonctions_beta_sortie = fonctions_beta;

tableau_affectation_intervalles_sortie = tableau_affectation_methode_1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ajout_impulsion_finale_nulle == 'oui'
   
   nbr_fct_beta_exacte = size( fonctions_beta_sortie , 1 );
   
   t0_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,1);
   tc_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,2);
   t1_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,3);
   p_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,4);
   q_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,5);
   K_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,6);
   rong_extremum_derniere_fct_beta = fonctions_beta(nbr_fct_beta_exacte,7);
   
   t0_fct_beta_immobile = tc_derniere_fct_beta;
   tc_fct_beta_immobile = t1_derniere_fct_beta;
   t1_fct_beta_immobile = tc_fct_beta_immobile + ( (tc_fct_beta_immobile - t0_fct_beta_immobile) / 3 );
   p_fct_beta_immobile = 2;
   q_fct_beta_immobile = (p_fct_beta_immobile * (((t1_fct_beta_immobile - tc_fct_beta_immobile)/(tc_fct_beta_immobile - t0_fct_beta_immobile))) );
   K_fct_beta_immobile = 0.000001;
   rong_extremum_fct_beta_immobile = rong_extremum_derniere_fct_beta + 1;
   
   
   
   fct_beta_immobile = [t0_fct_beta_immobile, tc_fct_beta_immobile, t1_fct_beta_immobile, p_fct_beta_immobile, q_fct_beta_immobile, K_fct_beta_immobile, rong_extremum_fct_beta_immobile ];
   
   
   fonctions_beta_sortie = [fonctions_beta_sortie; fct_beta_immobile];
   %tableau_affectation_intervalles_sortie = [tableau_affectation_intervalles_sortie; ]
   
   
    
end





          %HB% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %HB%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          %HB% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %HB%


