function [max_min_x_y, DAT, points_filtre] = Filtrage_conventionnel_plus_Dat_et_Detection_lim_max_min(points, rayon_filtre_traject, sigma_p_filtre_traject)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction Filtrage_conventionnel_plus_Detection_lim_max_min permet de ......  %%%%%
%%%%%                                                                                 %%%%%
%%%%% sur une ecriture manuscrite arbe en ligne ou squeletisée                        %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_teta_LB = 0.2;%0.1;
min_teta_LB = - max_teta_LB;

n_saut = 1; % 2; 3;

k1 = 1;
k2 = size(points,1);

for i = 1 : 7
%    rayon_filtre_traject = 3%16%8%16;
%    sigma_p_filtre_traject = 3.5%5%5.5%2%3.5%2;
    if ((rayon_filtre_traject * 2) + 1) > k2
%       rayon = max( round( k2 / 2 ) , 1);
       rayon_filtre_traject = k2;
    end
    [points,XYpoint] = filtre_lineaire_1(rayon_filtre_traject, sigma_p_filtre_traject, 1, size(points,1), points);
end

points_filtre = points;

max_x_filtre = max(points_filtre(:,1));
min_x_filtre = min(points_filtre(:,1));
max_y_filtre = max(points_filtre(:,2));
min_y_filtre = min(points_filtre(:,2));

max_min_x_y = [max_x_filtre , min_x_filtre , max_y_filtre , min_y_filtre];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Angle d'inclinaison de la tangente au tracé, DAT. %%%%%%%%%%%%

[DAT, memDAT, type_point, Ttemps] = Direction_Angulaire_continue_de_la_Tangente_DAT_Abd_AL_KARIM(k1, k2, points);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% ******************************************************* %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

