function Verification_rapport_vitesses_axes_Beta_chevauchees(fonctions_beta, param_trajectoire, tableau_affectation_intervalles, V, DAT, KMH_deb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Vérification du rapport ( amplitudes de vitesse / axes elliptiques ) %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbr_intervalles = size(tableau_affectation_intervalles, 1);
Vmax_HB = [];
Vmin_HB = [];
Mp_KMH_HB = [];
Mg_KMH_HB = [];
for i = 1 : nbr_intervalles

    K1 = tableau_affectation_intervalles(i, 5);
    K2 = tableau_affectation_intervalles(i, 6);
    K_max = max( K1 , K2 );
    K_min = min( K1 , K2 );
    Vmax_HB = [Vmax_HB; K_max];
    Vmin_HB = [Vmin_HB; K_min];
    if K_max == K1
       Mg_KMH_HB = [Mg_KMH_HB; tableau_affectation_intervalles(i,7)];
       Mp_KMH_HB = [Mp_KMH_HB; tableau_affectation_intervalles(i,8)];
    else
       Mg_KMH_HB = [Mg_KMH_HB; tableau_affectation_intervalles(i,8)];
       Mp_KMH_HB = [Mp_KMH_HB; tableau_affectation_intervalles(i,7)];
    end

end

A = param_trajectoire(:,1);
B = param_trajectoire(:,2);
TETA = param_trajectoire(:,3);
TETA_p = param_trajectoire(:,4);

a_quart_e_HB = param_trajectoire(:,5);
b_quart_e_HB = param_trajectoire(:,8);
a_oblique_HB = param_trajectoire(:,6);
b_oblique_HB = param_trajectoire(:,9);
a_2pts2tg_HB = param_trajectoire(:,7);
b_2pts2tg_HB = param_trajectoire(:,10);


Alpha_HB2 = [];
for i = 1 : nbr_intervalles
    alpha_hb = abs(TETA(i) - TETA_p(i));
    Alpha_HB2 = [Alpha_HB2; alpha_hb];
end

K1_HB = tableau_affectation_intervalles(:, 5);
K2_HB = tableau_affectation_intervalles(:, 6);
Alpha_HB = Alpha_HB2;
if nbr_intervalles > 4
   i = 2;
   while  i <= nbr_intervalles - 1
       if (K1_HB(i) > K2_HB(i)) & (K1_HB(i+1) < K2_HB(i+1))
          alpha_HB = abs(TETA(i + 1) - TETA(i)) / 2;
          Alpha_HB(i) = alpha_HB;
          Alpha_HB(i+1) = alpha_HB;
          i = i + 1;
        end
        i = i + 1;
   end
end

a_HB = a_quart_e_HB;
b_HB = b_quart_e_HB;
Test_HB_1 = [];
for i = 1 : nbr_intervalles
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

    
a_HB = a_oblique_HB;
b_HB = b_oblique_HB;
Test_HB_2 = [];
for i = 1 : nbr_intervalles
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

a_HB = a_2pts2tg_HB;
b_HB = a_2pts2tg_HB;
Test_HB_3_1 = [];
Test_HB_3 = [];

for i = 1 : nbr_intervalles
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
    
    Test_HB_3_1 = [ Test_HB_3; produit_1 produit_2 produit_3 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    UU = 4;
    if ( abs(produit_3 - produit_1) / max(produit_3 , produit_1) )* 100 > 30
       indDAT_G = Mg_KMH_HB(i) - KMH_deb + 1;
       indDAT_P = Mp_KMH_HB(i) - KMH_deb + 1;
       indDAT_p_UU = indDAT_P + (UU*(sign(indDAT_G - indDAT_P))) ;
       if (indDAT_p_UU == 0 )
          indDAT_p_UU = 1;
       elseif indDAT_p_UU > length(DAT)
          indDAT_p_UU = length(DAT);
       end
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


Test_HB_Final = [];
for i = 1 : nbr_intervalles
    produit_1_Test2 = Test_HB_2(i,1);
    produit_3_Test2 = Test_HB_2(i,3);
    produit_1_Test3 = Test_HB_3(i,1);
    produit_3_Test3 = Test_HB_3(i,3);
    rap_diff_produit_1_produit_3_Test2 = ( abs(produit_1_Test2 - produit_3_Test2) ) / ( max(produit_1_Test2 , produit_3_Test2) );
    rap_diff_produit_1_produit_3_Test3 = ( abs(produit_1_Test3 - produit_3_Test3) ) / ( max(produit_1_Test3 , produit_3_Test3) );
    if rap_diff_produit_1_produit_3_Test2 < rap_diff_produit_1_produit_3_Test3
       ligne_tab_Test = Test_HB_2(i,:);
    else
       ligne_tab_Test = Test_HB_3(i,:);
    end

    Test_HB_Final = [ Test_HB_Final; ligne_tab_Test ];
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Test_HB_1
pause;
Test_HB_2
pause;
Test_HB_3
pause
Test_HB_Final
pause;
pause;

