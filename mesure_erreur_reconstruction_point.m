function [dist_Erreur] = mesure_erreur_reconstruction_point(xM,yM,parametre_ellipse)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction mesure_erreur_reconstruction_point permet de calculer                              %%%%%
%%%%% l'erreur de reconstruction en distance d'un point donné d'un segment de la trajectoire         %%%%%
%%%%% par rapport à son image sur l'arc d'ellipse modélisant ce segment et                           %%%%%
%%%%% défini par les parmètres de son équation.                                                      %%%%%
%%%%% Erreur de reconstruction représente la différence entre les deux rayons qui divergent          %%%%%
%%%%% du centre de l'ellipse et atteignent respectivement l'arc et la trajectoire.                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = parametre_ellipse(1);
b = parametre_ellipse(2);
x0 = parametre_ellipse(3);
y0 = parametre_ellipse(4);
teta = parametre_ellipse(5);

[ha]= angle_reyon_ellipse(xM,yM,x0,y0,teta);

tg_ha = tan(ha);
               
XX = sqrt(((a^2)*(b^2))/((b^2)+((a*tg_ha)^2)));
YY = XX * tg_ha;
               
RO = sqrt((XX^2) + (YY^2));

dist_MO = sqrt( ((xM - x0)^2) + ((yM - y0)^2) );

dist_Erreur = abs( dist_MO - RO );

        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %HB%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%HB%

