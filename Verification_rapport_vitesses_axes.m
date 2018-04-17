function Verification_rapport_vitesses_axes(IMPULSION_BETA, ENTRAINEMENT, param_trajectoire, V, DAT, KMH_deb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Vérification du rapport ( amplitudes de vitesse / axes elliptiques ) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
for  i = 1 : nbr_seg
     delta_T = T1(i) - T0(i);
     DELTA_T = [DELTA_T; delta_T];
end


A1 = param_trajectoire(:,1);
B1 = param_trajectoire(:,2);
B2 = param_trajectoire(:,3);
TETA = param_trajectoire(:,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vmax_HB = [];
Kc = IMPULSION_BETA(:,8);
for  i = 1 : nbr_seg
     Vmax_HB = [Vmax_HB; Kc(i); Kc(i)];
end

Vmin_HB = [];
for  i = 1 : nbr_seg
     Vmin_HB = [Vmin_HB; Vini(i); Vfin(i)];
end

A2 = param_trajectoire(:,5);
a_HB = [];
for  i = 1 : nbr_seg
     a_HB = [a_HB; A1(i); A2(i)];
end

b_HB = [];
for  i = 1 : nbr_seg
     b_HB = [b_HB; B1(i); B2(i)];
end

Alpha_HB = [pi/2];
if nbr_seg >= 2
for  i = 2 : nbr_seg
     alpha_HB = abs(TETA(i) - TETA(i-1)) / 2;
     Alpha_HB = [Alpha_HB; alpha_HB; alpha_HB];
end
end
Alpha_HB = [Alpha_HB; pi/2];


teta_p1 = param_trajectoire(:,6);
teta_p2 = param_trajectoire(:,7);
teta_p1_HB = [];
for  i = 1 : nbr_seg
     teta_p1_HB = [teta_p1_HB; teta_p1(i)];
end
teta_p2_HB = [];
for  i = 1 : nbr_seg
     teta_p2_HB = [teta_p2_HB; teta_p2(i)];
end


Alpha_HB2 = [];
for  i = 1 : (nbr_seg)
     alpha_HB2_1 = abs(TETA(i) - teta_p1_HB(i));
     alpha_HB2_2 = abs(TETA(i) - teta_p2_HB(i));
     Alpha_HB2 = [Alpha_HB2; alpha_HB2_1; alpha_HB2_2];
end


DeltaT_HB = [];
for  i = 1 : nbr_seg
     DeltaT_HB = [DeltaT_HB; Tc(i) - T0(i); T1(i) - Tc(i)];
end



a_quart_HB_1 = param_trajectoire(:,8);
b_quart_HB_1 = param_trajectoire(:,11);
a_quart_HB_2 = param_trajectoire(:,14);
b_quart_HB_2 = param_trajectoire(:,17);
a_HB = []; b_HB = [];
for  i = 1 : (nbr_seg)
     a_HB = [a_HB; a_quart_HB_1(i); a_quart_HB_2(i)];
     b_HB = [b_HB; b_quart_HB_1(i); b_quart_HB_2(i)];
end

Test_HB_1 = [];
for i = 1 : (2 * nbr_seg)
    produit_1 = Vmax_HB(i) * b_HB(i);
    produit_2 = Vmin_HB(i) * a_HB(i);

    facteur_HB = sin( Alpha_HB(i) ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha_HB(i) ) ) ) );
    facteur_HB2 = sin( Alpha_HB2(i) ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha_HB2(i) ) ) ) );
    produit_3_1 = produit_2 * abs( facteur_HB );
    produit_3_2 = produit_2 * abs( facteur_HB2 );
    diff_produit_1_produit_3_1 = abs(produit_1 - produit_3_1);
    diff_produit_1_produit_3_2 = abs(produit_1 - produit_3_2);
    if min(diff_produit_1_produit_3_1, diff_produit_1_produit_3_2) == diff_produit_1_produit_3_1
       produit_3 = produit_3_1;
    else
       produit_3 = produit_3_2;
    end
    
    Test_HB_1 = [ Test_HB_1; produit_1 produit_2 produit_3 ];
end


a_oblique_HB_1 = param_trajectoire(:,9);
b_oblique_HB_1 = param_trajectoire(:,12);
a_oblique_HB_2 = param_trajectoire(:,15);
b_oblique_HB_2 = param_trajectoire(:,18);
a_HB = []; b_HB = [];
for  i = 1 : (nbr_seg)
     a_HB = [a_HB; a_oblique_HB_1(i); a_oblique_HB_2(i)];
     b_HB = [b_HB; b_oblique_HB_1(i); b_oblique_HB_2(i)];
end

Test_HB_2 = [];
for i = 1 : (2 * nbr_seg)
    produit_1 = Vmax_HB(i) * b_HB(i);
    produit_2 = Vmin_HB(i) * a_HB(i);
    facteur_HB = sin( Alpha_HB(i) ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha_HB(i) ) ) ) );
    facteur_HB2 = sin( Alpha_HB2(i) ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha_HB2(i) ) ) ) );
    produit_3_1 = produit_2 * abs( facteur_HB );
    produit_3_2 = produit_2 * abs( facteur_HB2 );
    diff_produit_1_produit_3_1 = abs(produit_1 - produit_3_1);
    diff_produit_1_produit_3_2 = abs(produit_1 - produit_3_2);
    if min(diff_produit_1_produit_3_1, diff_produit_1_produit_3_2) == diff_produit_1_produit_3_1
       produit_3 = produit_3_1;
    else
       produit_3 = produit_3_2;
    end
    
    Test_HB_2 = [ Test_HB_2; produit_1 produit_2 produit_3 ];
end

Mp_KMH_HB = [];
Mg_KMH_HB = [];
for  i = 1 : (nbr_seg)
     Mp_KMH_HB = [Mp_KMH_HB; param_trajectoire(i,20); param_trajectoire(i,22)];
     Mg_KMH_HB = [Mg_KMH_HB; param_trajectoire(i,21); param_trajectoire(i,21)];
end


a_2pts2tg_HB_1 = param_trajectoire(:,10);
b_2pts2tg_HB_1 = param_trajectoire(:,13);
a_2pts2tg_HB_2 = param_trajectoire(:,16);
b_2pts2tg_HB_2 = param_trajectoire(:,19);
a_HB = []; b_HB = [];
for  i = 1 : (nbr_seg)
     a_HB = [a_HB; a_2pts2tg_HB_1(i); a_2pts2tg_HB_2(i)];
     b_HB = [b_HB; b_2pts2tg_HB_1(i); b_2pts2tg_HB_2(i)];
end

Test_HB_3_1 = [];
Test_HB_3 = [];

for i = 1 : (2 * nbr_seg)
    produit_1 = Vmax_HB(i) * b_HB(i);
    produit_2 = Vmin_HB(i) * a_HB(i);
    facteur_HB = sin( Alpha_HB(i) ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha_HB(i) ) ) ) );
    facteur_HB2 = sin( Alpha_HB2(i) ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha_HB2(i) ) ) ) );
    produit_3_1 = produit_2 * abs( facteur_HB );
    produit_3_2 = produit_2 * abs( facteur_HB2 );
    diff_produit_1_produit_3_1 = abs(produit_1 - produit_3_1);
    diff_produit_1_produit_3_2 = abs(produit_1 - produit_3_2);
    if min(diff_produit_1_produit_3_1, diff_produit_1_produit_3_2) == diff_produit_1_produit_3_1
       produit_3 = produit_3_1;
    else
       produit_3 = produit_3_2;
    end
    
    Test_HB_3_1 = [ Test_HB_3_1; produit_1 produit_2 produit_3 ];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    UU = 4;
    if ( abs(produit_3 - produit_1) / max(produit_3 , produit_1) )* 100 > 30
       
       indDAT_G = Mg_KMH_HB(i) - KMH_deb + 1;
       indDAT_P = Mp_KMH_HB(i) - KMH_deb + 1;
       indDAT_p_UU = indDAT_P + (UU*(sign(indDAT_G - indDAT_P))) ;
       Vmin = V( indDAT_p_UU );
       Alpha = abs( TETA( round((i+0.1)/2) ) - DAT(indDAT_p_UU) );
       produit_2_UU = Vmin * a_HB(i);
       facteur_HB_UU = sin( Alpha ) /( sin( atan( (a_HB(i)/b_HB(i)) * tan( Alpha ) ) ) );
       produit_3_UU = produit_2_UU * abs( facteur_HB_UU );
       diff_produit_1_produit_3 = abs(produit_1 - produit_3);
       diff_produit_1_produit_3_UU = abs(produit_1 - produit_3_UU);
       if min(diff_produit_1_produit_3, diff_produit_1_produit_3_UU) == diff_produit_1_produit_3
          produit_3 = produit_3;
       else
          produit_3 = produit_3_UU;
       end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Test_HB_3 = [ Test_HB_3; produit_1 produit_2 produit_3 ];

end

Test_HB_1
pause;
Test_HB_2
pause;
Test_HB_3
pause;
pause;

