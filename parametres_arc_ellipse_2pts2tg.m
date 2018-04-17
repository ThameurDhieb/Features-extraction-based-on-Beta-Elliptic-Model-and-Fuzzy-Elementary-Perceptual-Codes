function [aa,bb,xx0,yy0,tteta,drap,erreur_MAX,ordre_ok,hA] = parametres_arc_ellipse_2pts2tg(A,B,X0,Y0,Teta_G,Teta_P,M_2points,M_2KMH,KMH_deb, DAT,num_fig,RRR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cette fonction permet de calculer les paramètres de l'équation de l'ellipse modélisant un segment %%%
%%%%%%%%% de la trajectoire dont le grand axe ait la direction de la tangente au point limite  %%%%%%%%%%
%%%%%% correspondant au maximum de vitesse et soit tangente aussi à la trajectoire au dexième %%%%%%%%%%%
%%%%%%%%%%%%%%         point limite du segment correspondant au minimum de vitesse         %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


erreur_MAX = 0;
ordre_ok = 0;
ang_dev_lim_sup = (pi/2) + (pi/11);%(pi/10);%(pi/8);
ang_dev_lim_sup_Max = (pi/2) + (pi/4);

%ang_dev_lim_inf = ( atan(B/A) );
bb_max = ( max(A,B) )  * ( 100 );%50 );%100 );%200 );

if A == 0
   A = 1.0000e-10;
end
if bb_max == B
   bb_max = bb_max + 1.0000e-10;
end
ang_dev_lim_inf = ( atan( ( (bb_max^2 / (bb_max - B)) - (bb_max - B) ) / A) ); 

ang_dev_inf = pi/24;%pi/12;

 xMG = M_2points(1);
 yMG = M_2points(2);
 xMP = M_2points(3);
 yMP = M_2points(4);
 
 KMH_G = M_2KMH(1);
 KMH_P = M_2KMH(2);
 indDAT_G = KMH_G - KMH_deb + 1;
 indDAT_P = KMH_P - KMH_deb + 1;
 indDAT_1 = min(indDAT_G, indDAT_P);
 indDAT_2 = max(indDAT_G, indDAT_P);
 u = 1;
 L_DAT = length(DAT);
 
ang_dev = Teta_G - Teta_P;
UU = 2;
if abs(ang_dev) > ang_dev_lim_sup_Max;
   if (indDAT_1 == indDAT_G)&((indDAT_P - UU) > 0)
      Teta_P = DAT(indDAT_P - UU);
   elseif (indDAT_1 == indDAT_P)&((indDAT_P + UU) <= L_DAT)
      Teta_P = DAT(indDAT_P + UU);
   end
   ang_dev = Teta_G - Teta_P;
   hhh = 0;
   while ( abs(ang_dev) < ang_dev_lim_inf ) & ( abs(ang_dev) < pi/2 )&(hhh < 4)
      Teta_P = ( Teta_P + DAT(indDAT_P) )/2;
      ang_dev = Teta_G - Teta_P;
      hhh = hhh + 1;
   end
end

ang_dev = Teta_G - Teta_P;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                   % Traitement des cas particuliers d'inflection rapide (non détectée) de la %
                              % trajectoire avant ou aprés un extremum de vitesse %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UU = 2;
if ( abs(indDAT_G - indDAT_P) > 2 * UU ) 
   if indDAT_1 == indDAT_P
      if (indDAT_G - UU > 0)&(indDAT_P + UU <= L_DAT)
         Teta_G_intermediaire = DAT(indDAT_G - UU);
         Teta_P_intermediaire = DAT(indDAT_P + UU);
      else
         Teta_G_intermediaire = Teta_G;
         Teta_P_intermediaire = Teta_P;
      end
   else
      if (indDAT_P - UU > 0)&(indDAT_G + UU <= L_DAT)
         Teta_G_intermediaire = DAT(indDAT_G + UU);
         Teta_P_intermediaire = DAT(indDAT_P - UU);
      else
         Teta_G_intermediaire = Teta_G;
         Teta_P_intermediaire = Teta_P;
      end

   end
   ang_dev_intermediaire1 = (Teta_G - Teta_P_intermediaire);
   ang_dev_intermediaire2 = (Teta_G_intermediaire - Teta_P);
   if ( abs(ang_dev_intermediaire1) > abs(ang_dev) )
                                                              % cas particulier ' inflexion rapide de la trajectoire 
      Teta_P = Teta_P_intermediaire;                          % apres un minimum de vitesse

   elseif ( abs(ang_dev_intermediaire2) > abs(ang_dev) )
                                                              % cas particulier ' inflexion rapide de la trajectoire 
      Teta_G = Teta_G_intermediaire;                          % apres un maximum de vitesse
      [A,B,X0,Y0,Teta_G] = parametres_quart_droit(xMG,yMG,xMP,yMP,Teta_G);
   end
end

UU1 = 4;
if abs(indDAT_G - indDAT_P) > UU1 + UU
 UU1 = 4;
 poursuite = 1;
 while ( poursuite == 1 )&( UU1 < abs(indDAT_G - indDAT_P) )
   ang_dev = Teta_G - Teta_P;
   if indDAT_1 == indDAT_P
      if (indDAT_P + UU1 <= L_DAT)
         Teta_P_intermediaire = DAT(indDAT_P + UU1);
      else
         Teta_P_intermediaire = Teta_P;
      end
   else
      if (indDAT_P - UU1 > 0)
         Teta_P_intermediaire = DAT(indDAT_P - UU1);
      else
         Teta_P_intermediaire = Teta_P;
      end

   end
   ang_dev_intermediaire1 = (Teta_G - Teta_P_intermediaire);
   if ( abs(ang_dev_intermediaire1) > abs(ang_dev) )
                                                              % cas particulier ' inflexion rapide de la trajectoire 
      Teta_P = Teta_P_intermediaire;                          % apres un minimum de vitesse
      poursuite = 1;
   else
      Teta_P = Teta_P;
      poursuite = 0;
   end
   UU1 = UU1 + UU;
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[haP]= angle_reyon_ellipse(xMP,yMP,xMG,yMG,Teta_G);

if haP >= (3*pi/2)
   if ( ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev <= (ang_dev_lim_inf) ) )
      Teta_P = Teta_G - ( 1.3 * ang_dev_lim_inf );
   elseif ( abs( ang_dev ) >= ang_dev_lim_sup_Max )
      Teta_P = Teta_G - ( 1.3 * ( atan(B/A) ) );
   end
elseif (haP <= (3*pi/2))&(haP >= (pi))
   if ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev >= - (ang_dev_lim_inf) )
      Teta_P = Teta_G + ( 1.3 * ang_dev_lim_inf );
   elseif ( abs( ang_dev ) >= ang_dev_lim_sup_Max )
      Teta_P = Teta_G + ( 1.3 * ( atan(B/A) ) );
   end
elseif (haP <= (pi))&(haP >= (pi/2))
   if ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev <= (ang_dev_lim_inf) )
      Teta_P = Teta_G - ( 1.3 * ang_dev_lim_inf );
   elseif ( abs( ang_dev ) >= ang_dev_lim_sup_Max )
      Teta_P = Teta_G - ( 1.3 * ( atan(B/A) ) ); 
   end
elseif (haP <= (pi/2))
   if ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev >= - (ang_dev_lim_inf) )
      Teta_P = Teta_G + ( 1.3 * ang_dev_lim_inf );
   elseif ( abs( ang_dev ) >= ang_dev_lim_sup_Max )
      Teta_P = Teta_G + ( 1.3 * ( atan(B/A) ) );
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ang_dev = Teta_G - Teta_P;
alfa = (Teta_G + (pi/2)) - Teta_P;


if (abs(ang_dev) > ang_dev_lim_sup)
   recip = (sign(ang_dev)) * ang_dev_lim_sup;
   alfa = recip + (pi/2);
elseif ( abs(ang_dev) < ang_dev_lim_inf ) & ( abs(ang_dev) < pi/2 )
   alfa = (pi/2) + ((sign(ang_dev)) * ang_dev_inf);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

teta_M2_Ref1 = (alfa - (pi/2));

if haP >= (3*pi/2)
   if abs( ang_dev ) <= (pi/2)
      tgM2 = - abs( tan(teta_M2_Ref1) );
   else
      tgM2 =   abs( tan(teta_M2_Ref1) );
   end
elseif (haP <= (3*pi/2))&(haP >= (pi))
   if abs( ang_dev ) <= (pi/2)
      tgM2 =   abs( tan(teta_M2_Ref1) );
   else
      tgM2 = - abs( tan(teta_M2_Ref1) );
   end
elseif (haP <= (pi))&(haP >= (pi/2))
   if abs( ang_dev ) <= (pi/2)
      tgM2 = - abs( tan(teta_M2_Ref1) );
   else
      tgM2 =   abs( tan(teta_M2_Ref1) );
   end
elseif (haP <= (pi/2))
   if abs( ang_dev ) <= (pi/2)
      tgM2 =   abs( tan(teta_M2_Ref1) );
   else
      tgM2 = - abs( tan(teta_M2_Ref1) );
   end
end

if tgM2 == 0
   tgM2 = 1.0000e-10;
end

ABsurtgM2 = (A * B) / tgM2;
Expression1 = (A^2) - (2 * ABsurtgM2);
Expression2 = (A^2) - (ABsurtgM2);

if Expression2 == 0
   Expression2 = 1.0000e-10;
end
if B == 0
   B = 1.0000e-10;
end

ALFA = Expression1 / (Expression2^2);
BETA = (Expression1 / ( Expression2 * B))^2;
if ALFA == 0
   ALFA = 1.0000e-10;
end
if BETA == 0
   BETA = 1.0000e-10;
end


racine_SIGMA = ABsurtgM2 / ( Expression2 );


aa = sqrt( 1 / ALFA );
bb = sqrt( 1 / BETA );

if ( (bb < B)&( abs( ang_dev ) < (pi/2) ) ) | ( (bb > B)&( abs( ang_dev ) > (pi/2) ) )
   tgM2 = - tgM2;
   
   ABsurtgM2 = (A * B) / tgM2;
   Expression1 = (A^2) - (2 * ABsurtgM2);
   Expression2 = (A^2) - (ABsurtgM2);
   
   if Expression2 == 0
      Expression2 = 1.0000e-10;
   end

   ALFA = Expression1 / (Expression2^2);
   BETA = (Expression1 / ( Expression2 * B))^2;

   if ALFA == 0
      ALFA = 1.0000e-10;
   end
   if BETA == 0
      BETA = 1.0000e-10;
   end

   racine_SIGMA = ABsurtgM2 / ( Expression2 );

   aa = sqrt( 1 / ALFA );
   bb = sqrt( 1 / BETA );
end   
   
% détermination du signe de l'ordonné y0_Ref1 du centre de l'ellipse
% dans le refereciel Ref1 : (centre quart, teta_G) sachant que son absice
% x0_Ref1 = 0

if haP >= (3*pi/2)
   if abs( ang_dev ) <= (pi/2)
      y0_Ref1 = - bb * racine_SIGMA;
   else
      y0_Ref1 = bb * racine_SIGMA;
   end
elseif (haP <= (3*pi/2))&(haP >= (pi))
   if abs( ang_dev ) <= (pi/2)
      y0_Ref1 = - bb * racine_SIGMA;
   else
      y0_Ref1 = bb * racine_SIGMA;
   end
elseif (haP <= (pi))&(haP >= (pi/2))
   if abs( ang_dev ) <= (pi/2)
      y0_Ref1 = bb * racine_SIGMA;
   else
      y0_Ref1 = - bb * racine_SIGMA;
   end
elseif (haP <= (pi/2))
   if abs( ang_dev ) <= (pi/2)
      y0_Ref1 = bb * racine_SIGMA;
   else
      y0_Ref1 = - bb * racine_SIGMA;
   end
end

% détermination des coordonnées (xc_Ref0,yc_Ref0) du centre de l'ellipse 
% dans le referenciel Ref0 : (O(0,0),x,y) 

teta_M1_Ref0 = Teta_G;
tg_teta_M1_Ref0 = tan( teta_M1_Ref0 );
teta_M1_Ref0 = atan(tg_teta_M1_Ref0);

% Deux possibilités (solutions):
xc_Ref0_1 = X0 + ( y0_Ref1 * sin( teta_M1_Ref0 ) );
yc_Ref0_1 = Y0 - ( y0_Ref1 * cos( teta_M1_Ref0 ) );
xc_Ref0_2 = X0 - ( y0_Ref1 * sin( teta_M1_Ref0 ) );
yc_Ref0_2 = Y0 + ( y0_Ref1 * cos( teta_M1_Ref0 ) );


% Détermination de la bonne solution (coicidence avec les points MP et MG): 
dist_MCMG_1 = sqrt( ((xMG - xc_Ref0_1)^2) + ((yMG - yc_Ref0_1)^2) );
dist_MCMG_2 = sqrt( ((xMG - xc_Ref0_2)^2) + ((yMG - yc_Ref0_2)^2) );
diff_1 = abs(bb - dist_MCMG_1);
diff_2 = abs(bb - dist_MCMG_2);

if diff_1 == min(diff_1, diff_2)
   xx0 = xc_Ref0_1;
   yy0 = yc_Ref0_1;
else
   xx0 = xc_Ref0_2;
   yy0 = yc_Ref0_2;
end

tteta = teta_M1_Ref0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%HBH% [ok,erreur_MAX,ordre_ok] = verification_des_fuites7(a,b,x0,y0,teta,M_5points,M_5KMH); 
%HBH% if (ok == 1)%|(erreur_MAX < 0.05)
            drap = 1;            
            
            [haG]= angle_reyon_ellipse(xMG,yMG,xx0,yy0,tteta);
            [haP]= angle_reyon_ellipse(xMP,yMP,xx0,yy0,tteta);
            hA = [haG, haP];
            
%HBH% else
%HBH%    drap = 0;
%HBH%    a = 0; b = 0; x0 = 0; y0 = 0; teta = 0;
%HBH%    HA = [0, 0, 0];
%HBH% end
        
ok = 1;
erreur_MAX = 0;
ordre_ok = 1;



        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
