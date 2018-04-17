function [a,b,x0,y0,teta,drap,erreur_MAX,ordre_ok,hA] = parametres_quart_ellipse_oblique(A,B,X0,Y0,Teta_G,Teta_P,M_2points,M_2KMH,KMH_deb, DAT,num_fig,RRR)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cette fonction permet de calculer les paramètres de l'équation de l'ellipse modélisant un segment %%%
%%%%%%%%%%%%%%%%    de la trajectoire en appliquant une projection oblique à un cercle    %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%    de rayon R que l'on calculera telque l'ellipse obtenue soit tangent   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%          au segment de la trajectoire en ses deux point limites          %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


erreur_MAX = 0;
ordre_ok = 0;
ang_dev_lim_sup = (pi/2) + (pi/10);%(pi/11);%(pi/30);%(pi/8);
ang_dev_lim_sup = (pi/2) + (pi/11);%(pi/30);%(pi/8);
ang_dev_lim_sup_Max = (pi/2) + (pi/4);

if A == 0
   A = 1.0000e-10;
end
ang_dev_lim_inf = atan(B/A);
ang_dev_inf = pi/24;% pi/36; % pi/12;


 xMG = M_2points(1);
 yMG = M_2points(2);
 xMP = M_2points(3);
 yMP = M_2points(4);
 
 Teta_G_init = Teta_G;
 
 KMH_G = M_2KMH(1);
 KMH_P = M_2KMH(2);
 indDAT_P = KMH_P - KMH_deb + 1;
 indDAT_G = KMH_G - KMH_deb + 1;
 indDAT_1 = min(indDAT_G, indDAT_P);
 indDAT_2 = max(indDAT_G, indDAT_P);
 u = 1;
 L_DAT = length( DAT );

ang_dev = Teta_G - Teta_P;

UU = 2;
if abs(ang_dev) > ang_dev_lim_sup_Max
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
UU = 2;

if ( abs(indDAT_G - indDAT_P ) > 2 * UU )

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
%      Teta_P = DAT(indDAT_P + 3);%Teta_G_intermediaire;
      [A,B,X0,Y0,Teta_G] = parametres_quart_droit(xMG,yMG,xMP,yMP,Teta_G);
   end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[haP]= angle_reyon_ellipse(xMP,yMP,xMG,yMG,Teta_G);

if haP >= (3*pi/2)
   if ( ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev <= (ang_dev_lim_inf) ) )
      Teta_P = Teta_G - ( 1.3 * ang_dev_lim_inf );
   end
elseif (haP <= (3*pi/2))&(haP >= (pi))
   if ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev >= - (ang_dev_lim_inf) )
      Teta_P = Teta_G + ( 1.3 * ang_dev_lim_inf );
   end
elseif (haP <= (pi))&(haP >= (pi/2))
   if ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev <= (ang_dev_lim_inf) )
      Teta_P = Teta_G - ( 1.3 * ang_dev_lim_inf );
   end
elseif (haP <= (pi/2))
   if ( abs( ang_dev ) <= (pi/2) ) & ( ang_dev >= - (ang_dev_lim_inf) )
      Teta_P = Teta_G + ( 1.3 * ang_dev_lim_inf );
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ang_dev = Teta_G - Teta_P;

if Teta_G < Teta_P
   alfa = (Teta_G + (pi/2)) - Teta_P;
else
   alfa = (Teta_G - (pi/2)) - Teta_P;
end


if abs(ang_dev) > ang_dev_lim_sup;
   recip = (sign(ang_dev)) * ang_dev_lim_sup;
   alfa = recip + (pi/2);
elseif ( abs(ang_dev) < ang_dev_lim_inf ) & ( abs(ang_dev) < pi/2 )
   alfa = (pi/2) + ((sign(ang_dev)) * ang_dev_inf);
end

tgMG = tan( Teta_G );

if tgMG >= 0
   if (abs(ang_dev) <= (pi/2))&(X0 <= xMP)
      x0 = X0 + (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 + (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   elseif (abs(ang_dev) <= (pi/2))&(X0 >= xMP)
      x0 = X0 - (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 - (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   elseif (abs(ang_dev) >= (pi/2))&(X0 <= xMP)
      x0 = X0 - (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 - (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   elseif (abs(ang_dev) >= (pi/2))&(X0 >= xMP)
      x0 = X0 + (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 + (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   else
      x0 = X0;
      y0 = Y0;
   end
else
   if (abs(ang_dev) <= (pi/2))&(X0 <= xMP)
      x0 = X0 + (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 - (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   elseif (abs(ang_dev) <= (pi/2))&(X0 >= xMP)
      x0 = X0 - (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 + (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   elseif (abs(ang_dev) >= (pi/2))&(X0 <= xMP)
      x0 = X0 - (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 + (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   elseif (abs(ang_dev) >= (pi/2))&(X0 >= xMP)
      x0 = X0 + (B * abs( (tan(alfa)) * (cos(Teta_G)) ) );
      y0 = Y0 - (B * abs( (tan(alfa)) * (sin(Teta_G)) ) );
   else
      x0 = X0;
      y0 = Y0;
   end
end

R = sqrt( ((x0 - xMP)^2) + ((y0 - yMP)^2) );

K = B / ( abs(R * cos(alfa)) );

gama = (1/2) * atan( ((2 * K)/(1 - (K^2))) * sin(alfa) );

ab1 = (1 + (K^2)) * (R^2) / 2;
ab2 = (1 - (K^2)) * (R^2) / 2;
ab3 = abs( ab2 * cos(2 * gama) );
ab4 = abs( (K * (R^2) * sin(alfa)) * sin(2 * gama) );

hb1 = ab1 + ab3 + ab4;
hb2 = ab1 - ab3 - ab4;
a = sqrt( hb1 );
b = sqrt( hb2 );

ab5 = cos(alfa) * sin(gama);
if K == 0
   K = 1.0000e-10;
end
ab6 = (cos(gama) / K) + (sin(alfa) * sin(gama));
if ab6 == 0
   ab6 = 1.0000e-10;
end

teta = ( atan(ab5 / ab6) ) + Teta_G;


if K > 1
   teta = teta + (pi/2);
end

 
%HBH% [ok,erreur_MAX,ordre_ok] = verification_des_fuites7(a,b,x0,y0,teta,M_5points,M_5KMH); 
%HBH% if (ok == 1)%|(erreur_MAX < 0.05)
            drap = 1;            

            [haG]= angle_reyon_ellipse(xMG,yMG,x0,y0,teta);
            [haP]= angle_reyon_ellipse(xMP,yMP,x0,y0,teta);
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
