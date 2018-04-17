function [IMPULSION_BETA, ENTRAINEMENT,vect_X0XCX1] = copie_extract_beta_mot_H6(V,DPV,points,T,t,RAP_lim,Sensibilite_plus, pas)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% la fonction extract_beta calcule les  5 parametres caracteristiques des differents beta             %%%%% 
%%%%% et les deux paramètres de la composante d'entrainement                                              %%%%%
%%%%% contenues dans un signal de vitesse défini                                                          %%%%%
%%%%% les parametres d'entree etant le signal de vitesse curviligne  filtré, lui meme                     %%%%%
%%%%% les paramètres de sortie sont les caractéristiques des Betas et de la composante d'entrainement     %%%%%
%%%%% constituant le signal de vitesse curviligne à savoir le t0, tc, t1, K_hauteur, p, q et V_ini, V_fin %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%HB%----------------------------------derivee seconde de la vitesse
delta_t = diff(t);
delta_DPV = diff(DPV);
DSV = (delta_DPV)./(delta_t);
DSV = [0; DSV];

%------------------------------------------------------------------

%------------------------------------longueur du vecteur vitesse
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
           v_max = [v_max; (round(XX(i-1)/pas) + 1) XX(i-1) YY(i-1)];
           extrema = [extrema; (round(XX(i-1)/pas) + 1) XX(i-1) YY(i-1) 2];
        end
     end
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
        end
     end
    
   i = i+1;   
end

%-------------------------------------- le dernier points
v_min = [v_min; (round(XX(i-1)/pas) + 1) XX(l) YY(l)];
extrema = [extrema; (round(XX(i-1)/pas) + 1) XX(l) YY(l) 0];                      
%--------------------------------------------------------------------------


p = 2; 
i = 1;
 
%--------------------- Elimination des points voisins ( bruit )
extremum = [extremum; extrema(1,:)];
for i = 1 : size(extrema,1)
    extremum = [extremum; extrema(i,:)];
    hh = extremum;
end   
   
j = 1;
while (j <= length(extremum)-2)

%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%
   if (extremum(j,4) == 0)
      uv = 1;
      while (extremum(j+uv,4) ~= 2)&(j+uv <= length(extremum)-1)
          uv = uv + 1;
      end    
      if (extremum(j+uv,4) == 2)&(j+uv <= length(extremum)-1)
         uvp = 1;
         while (extremum(j+uv+uvp,4) ~= 0)&(j+uv+uvp <= length(extremum))
            uvp = uvp + 1;
         end    
         if (extremum(j+uv+uvp,4) == 0)&(j+uv+uvp <= length(extremum))
            X0 = [X0; extremum(j,:)];
            XC = [XC; extremum(j+uv,:)];
            X1 = [X1; extremum(j+uv+uvp,:)];
         end
      end

      j = j+uv+uvp - 1;
   end    
%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%H%          
   j = j + 1;
end  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vect_X0XCX1 = [X0, XC, X1];

                %HB% %%%%%%%%%%%%%% Detection des points de distortion %%%%%%%%%%%%%%% %HB%
           %HB% %%%%%%%%%%%%%% Detection des points d'inflexion particuliers %%%%%%%%%%%%%%% %HB%

nbr_intervalles_Sensibilite_plus = size( Sensibilite_plus, 1);

Rap_lim1 = RAP_lim(1);      % seuil de detection grossiere des points de double inflexion en montée 
Rap_lim2 = RAP_lim(2);      % seuil de detection grossiere des points de double inflexion en descente 
           
 j = 1;
 nbr_int_init = length(XC(:,1));
 nbr_int_fin = nbr_int_init;
 l = length(YY);
 l_DPV = length(DPV);
 delta_i = 1;%2;
 delta_i2 = 2; %4; %1;

 while j <= nbr_int_fin,
   panier_X5 = [];
   panier_X6 = [];
 
     DPV_MMAX_LOCAL = 0;
     DPV_MMIN_LOCAL = 0;
     marque_fin_DPV = 0;
          
     for ih = X0(j,1)-X0(1,1)+1 : X1(j,1)-X0(1,1)

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
               end
            end
         end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %HB% Detection des points de dissymétrie  ,   codés [6]%
%          if [(DPV(ih) >= 0)&(DPV(ih)> DPV(ih-delta_i))&(DPV(ih) > DPV(ih+delta_i))]
%             DPV_max_local = DPV(ih);
%             if DPV_max_local > DPV_MMAX_LOCAL
%                DPV_MMAX_LOCAL = DPV_max_local;
%                i_DPV_MMAX_LOCAL = ih;
%             end
%          elseif [(DPV(ih) <= 0)&(DPV(ih)<= DPV(ih-delta_i))&(DPV(ih) <= DPV(i+delta_i))]|[(ih ==  l_DPV - 1 )&(DPV(ih) <= DPV(ih-delta_i))]
%             DPV_min_local = DPV(ih);
%             if DPV_min_local < DPV_MMIN_LOCAL
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
% %         if (RAP_symetrie < 0.4)&(RAP_DPV1 > 3)
%         if ( (RAP_symetrie < 0.4)&(RAP_DPV1 > 3) ) & ( (ih+X0(1,1)-1 > X0(j,1))&(ih+X0(1,1)-1 < X1(j,1))&(ih+X0(1,1)-1 ~= XC(j,1)) )
% 
%            long_demi_intervalle_impulsion = i_DPV_MMAX_LOCAL - (X0(j,1) - X0(1,1)+1);
%            i_fin_impulsion = i_DPV_MMAX_LOCAL + round((long_demi_intervalle_impulsion));
%            panier_X6 = [i_fin_impulsion+X0(1,1)  XX(i_fin_impulsion) YY(i_fin_impulsion) 6];
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
%             RAP_DPV2 = DPV(i_DPV_MMIN_LOCAL) / DPV(i_DPV_MMIN_LOCAL - delta_i2-1)
%          else
%             RAP_DPV2 = 1;
%          end
% %          if (RAP_symetrie > 10)&(RAP_DPV2 > 6)
%          if ( (RAP_symetrie > 10)&(RAP_DPV2 > 6) ) & ( (ih+X0(1,1)-1 > X0(j,1))&(ih+X0(1,1)-1 < X1(j,1))&(ih+X0(1,1)-1 ~= XC(j,1)) )
%             if (((X1(j,1)-X0(1,1)+1) - i_DPV_MMIN_LOCAL)*pas) == 0
%                long_demi_intervalle_impulsion = 3;
%             else
%                long_demi_intervalle_impulsion = (X1(j,1)-X0(1,1)+1) - i_DPV_MMIN_LOCAL;
%             end
%             i_debut_impulsion = i_DPV_MMIN_LOCAL - round((long_demi_intervalle_impulsion));
%             panier_X6 = [i_debut_impulsion+X0(1,1)  XX(i_debut_impulsion) YY(i_debut_impulsion) 6];
%          end
%      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mettre dans l'ordre les indices des nouveaux points détectés %
     panier_X5X6 = [];
     panier_X0XCX1_P = [];
     if (size(panier_X5,1) > 0)&(size(panier_X6,1) > 0)
        marque = 0;
        for jalon = 1 : size(panier_X5,1),
            if (panier_X5(jalon,1) > panier_X6(1))&(marque == 0)
               marque = 1;
               panier_X5X6 = [panier_X5X6; panier_X6];
            end
               panier_X5X6 = [panier_X5X6;panier_X5(jalon,:)];
        end
     elseif (size(panier_X5,1) == 0)&(size(panier_X6,1) > 0)
            panier_X5X6 = panier_X6;
     elseif (size(panier_X5,1) > 0)&(size(panier_X6,1) == 0)
            panier_X5X6 = panier_X5;
     end

% insérer les nouveaux points détectés dans les vecteurs des extremums

        reduct = 2.5;%4; %8;
        ind_X1_initial_pere = X1(j,1);

        longeur_TX0 = size(X0,1);
        nouveau_X0 = X0(j,:);

        Tab_X0 = [X0(1:j-1,:); nouveau_X0];
        Tab_XC = [XC(1:j-1,:)];
        Tab_X1 = [X1(1:j-1,:)];
        
        le_reel_size_de_panier_X5X6 = 0;                     % initialisation
        for jalon = 1 : size(panier_X5X6,1),
          ind_nouveau_X1 = panier_X5X6(jalon,1);
          ind_nouveau_X0 = nouveau_X0(1);
          if (abs(ind_nouveau_X1 - ind_nouveau_X0) >= 3)&(abs(ind_nouveau_X1 - ind_X1_initial_pere) >= 3)
            nouveau_X1 = panier_X5X6(jalon,:);
            Tab_X1 = [Tab_X1; nouveau_X1];
            if ((nouveau_X1(1) >= XC(j,1))&(nouveau_X0(1) <= XC(j,1)))
               panier_X0XCX1_P = [panier_X0XCX1_P; nouveau_X0 XC(j,:) nouveau_X1];
               Tab_XC = [Tab_XC; XC(j,:)];
            else

               if nouveau_X1(4) == 5
                  indice_nouveau_XC = round( nouveau_X1(1) - ((nouveau_X1(1) - nouveau_X0(1))/reduct) );
               else
                  indice_nouveau_XC = round((nouveau_X0(1) + nouveau_X1(1))/2);
               end

               
               %%% code des centres des nouveaux intervalles   [7] %%%
               
               nouveau_XC = [indice_nouveau_XC  XX(indice_nouveau_XC-X0(1,1)) YY(indice_nouveau_XC-X0(1,1)) 7];
               panier_X0XCX1_P = [panier_X0XCX1_P; nouveau_X0 nouveau_XC nouveau_X1];
               Tab_XC = [Tab_XC; nouveau_XC];
            end
            nouveau_X0 = panier_X5X6(jalon,:);
            Tab_X0 = [Tab_X0; nouveau_X0];

          
          le_reel_size_de_panier_X5X6 = le_reel_size_de_panier_X5X6 + 1;    % incrementation
          end

        end
        nouveau_X1 = X1(j,:);
        Tab_X1 = [Tab_X1; nouveau_X1];
        if ((nouveau_X1(1) >= XC(j,1))&(nouveau_X0(1) <= XC(j,1)))
           panier_X0XCX1_P = [panier_X0XCX1_P; nouveau_X0 XC(j,:) nouveau_X1];
           Tab_XC = [Tab_XC; XC(j,:)];
        else
           if nouveau_X1(4) == 5 
              indice_nouveau_XC = round( nouveau_X0(1) + ((nouveau_X1(1) - nouveau_X0(1))/reduct) );
           else
              indice_nouveau_XC = round((nouveau_X0(1) + nouveau_X1(1))/2);
           end

           nouveau_XC = [indice_nouveau_XC  XX(indice_nouveau_XC-X0(1,1)) YY(indice_nouveau_XC-X0(1,1)) 7];
           panier_X0XCX1_P = [panier_X0XCX1_P; nouveau_X0 nouveau_XC nouveau_X1];
           Tab_XC = [Tab_XC; nouveau_XC];
        end

        Tab_X0 = [Tab_X0; X0(j+1:longeur_TX0,:)];
        Tab_XC = [Tab_XC; XC(j+1:longeur_TX0,:)];
        Tab_X1 = [Tab_X1; X1(j+1:longeur_TX0,:)];

        X0 = Tab_X0;
        XC = Tab_XC;
        X1 = Tab_X1;


%%%%%%%%%%%%%%%%%%
        j = (j - 1) + 1 + le_reel_size_de_panier_X5X6 ;
        nbr_int_fin = nbr_int_fin + le_reel_size_de_panier_X5X6 ;
        
        j = j + 1;
 end
        %HB% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %HB%
  
vect_X0XCX1 = [X0, XC, X1];
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Définir la vitesse d'entrainement et l'impulsion beta  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%                 dans chaque intervalle                 %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 nombre_total_intervalle = length(XC(:,1));
 j=1 ;
 while (j <= nombre_total_intervalle),
       T0 = X0(j,2);
       TC = XC(j,2);
       T1 = X1(j,2);
       tct0 = TC - T0;
	   t1tc = T1 - TC;
           
    %----- Définir les paramètres de la vitesse d'entrainement ------%

       Vini = X0(j,3);
       Vfin = X1(j,3);
       T1_T0 = (X1(j,2) - X0(j,2));
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%      
       if T1_T0 == 0
          T1_T0 = 1;
       end
       %%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       aVL = (X1(j,3) - X0(j,3))/T1_T0;
       bVL = X0(j,3);
       tu = -6 * (X1(j,3) - Vini)/((T1_T0)^3);

    %---------- chercher les valeurs des paramétres p et q ----------%
    %----- à partir de données relevées sur le profile vitesse ------%
    %------- V(tchoix) = ychoix à un instant donné t = tchoix -------%
       if X1(j,1)-X0(j,1)+1 >= 3

          indice_tchoix1 = find(round(XX/pas) >= X0(j,1) + ((XC(j,1)- X0(j,1))/3));%1.3));
          tchoix1 = XX(indice_tchoix1(1));
          ychoix1 = YY(indice_tchoix1(1));

          tchoix1t0 = tchoix1 - T0;
          t1tchoix1 = T1 - tchoix1;

    %--------- Définir les paramètres de l'impulsion beta -----------%
 
          VCMtc = (tu*((((tct0)^3)/3) - ((T1-T0)*((tct0)^2)/2))) + Vini;             
          VCMtchoix1 = (tu*((((tchoix1t0)^3)/3) - ((T1-T0)*((tchoix1t0)^2)/2))) + Vini;             
          
          hauteur1 = abs(XC(j,3) - VCMtc);

          beta_choix1 = abs((YY(indice_tchoix1(1)) - VCMtchoix1)/hauteur1);

          if (beta_choix1 == 0)|(beta_choix1 == 1)
             beta_choix1 = YY(indice_tchoix1(1));
          end
          %% if t1tchoix1 == 0
          %%    t1tchoix1 = 1;
          %% end
          %% if tchoix1t0 == 0
          %%    tchoix1t0 = 1;
          %% end
          %% if t1tc == 0
          %%    t1tc = 1;
          %% end
          %% if tct0 == 0
          %%    tct0 = 1;
          %% end

          if t1tchoix1 == 0
             t1tchoix1 = pas/4;
          end
          if tchoix1t0 == 0
             tchoix1t0 = 3*pas/4;
          end
          if t1tc == 0
             t1tc = pas/2;
          end
          if tct0 == 0
             tct0 = pas/2;
          end

    %----------------------------------------------------------------%      
          %% if (tct0 == 0)
          %%    valeur_tct0 = 1.0000e-6;
          %% else
          %%    valeur_tct0 = tct0;
          %% end
          %% if (t1tc == 0)
          %%    valeur_t1tc = 1.0000e-6;
          %% else
          %%    valeur_t1tc = t1tc;
          %% end

          labez = ( tu*((t1tc)^2)*tct0/hauteur1 ) * ( log(t1tchoix1/t1tc) );
          p = abs( (log(beta_choix1) + labez)/((log(tchoix1t0/tct0))+((t1tc/tct0)*log(t1tchoix1/t1tc))) );
          %% labez = ( tu*((valeur_t1tc)^2)*valeur_tct0/hauteur1 ) * ( log(t1tchoix1/valeur_t1tc) );
          %% p = abs( (log(beta_choix1) + labez)/((log(tchoix1t0/valeur_tct0))+((valeur_t1tc/valeur_tct0)*log(t1tchoix1/valeur_t1tc))) );
                    
          p = real(p);
 
          %% if (p==Inf)|(p < 0)|(p==NaN)|(p==0)
          %%    p=2;
          %% end
          %% if (p < 0.7)
          %%    p = 0.7
          %% end
          %% if (p > 8)
          %%    p = 8
          %% end
          
    %------------------ cas d'un nouveau centre tc ------------------%
          if (XC(j,4) == 7 )
             if tct0 > t1tc
                p = 3;  q = (p*((t1tc/tct0)));
             else
                q = 3;  p = (q*((tct0/t1tc)));
             end
          else   
             q=(p*((t1tc/tct0))) - (tu*((t1tc)^2)*tct0/hauteur1); %HB% + ((Vini*x1xc)/hauteur1);

               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TesT
             %%%%%%%%%%% verification des valeurs des paramètres p et q %%%%%%%%%%

             Tt_betaMax = ((q * T0) + (p * T1))/(p + q);
             if (Tt_betaMax < (T1 - 1))&(Tt_betaMax > (T0 + 1))
                Tt_betaMax = Tt_betaMax;
             else
                Tt_betaMax = TC;
             end

             ttMt0 = Tt_betaMax - T0;
             t1ttM = T1 - Tt_betaMax;
             VPorteusePx_a_ttM = (tu*((((ttMt0)^3)/3) - ((T1-T0)*(((ttMt0)^2)/2)))) + Vini;
             Vimpulsive_betaPx_a_ttM = hauteur1*exp(log((((ttMt0./tct0)^p)*((t1ttM./t1tc)^q))));
    	     Vreconstruite_a_ttM = Vimpulsive_betaPx_a_ttM + VPorteusePx_a_ttM;
             V_reelle_a_ttM = YY( (round(Tt_betaMax/pas)) - X0(1,1)+2);
             if Vreconstruite_a_ttM > ((1.5)* V_reelle_a_ttM)
                indice_tchoix = find(round(XX/pas)>=(XC(j,1) + Tt_betaMax)/2);
                tchoix = XX(indice_tchoix(1));
                ychoix = YY(indice_tchoix(1));
                tchoixt0 = tchoix - T0;
                t1tchoix = T1 - tchoix;
                VCMtchoix = (tu*((((tchoixt0)^3)/3) - ((T1-T0)*((tchoixt0)^2)/2))) + Vini;             
                beta_choix = abs((YY(indice_tchoix(1)) - VCMtchoix)/hauteur1);
      
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
                if (beta_choix == 0)|(beta_choix == 1)
                    beta_choix = YY(indice_tchoix(1));
                end
                if t1tchoix == 0
                   t1tchoix = 1;
                end
                if tchoixt0 == 0
                   tchoixt0 = 1;
                end
                if t1tc == 0
                   t1tc = 1;
                end
                if tct0 == 0
                   tct0 = 1;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                p = abs(log(beta_choix)/(log((tchoixt0/tct0))+((t1tc/tct0)*log(t1tchoix/t1tc))));
                if (p == Inf)|(p < 0)
                   p = 2;
                end
                q = (p*((t1tc/tct0))) - (tu*((t1tc)^2)*tct0/hauteur1);

             end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TesT
               
          end          
          if (p==Inf)|(p==NaN)|(q==Inf)|(q==NaN)
             p = 2;
             q = 2;
          end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       else
          hauteur1 = 0;
          hauteur = 0;
          p = 2; q = 2;
          t1tc = 1; tct0 = 1;
       end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
       P=[P; p];
       Q=[Q; q];
       delta(j)=X1(j,2)-X0(j,2); 
       DELTA=[DELTA; delta(j)];
       HAUTEUR = [HAUTEUR; hauteur1];
       ENTRAINEMENT = [ENTRAINEMENT; tu, Vini, Vini, Vfin];

       
       j=j+1;
 end

 IMPULSION_BETA = [XC(:,2) DELTA P Q HAUTEUR X0(:,2) X1(:,2) XC(:,3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%