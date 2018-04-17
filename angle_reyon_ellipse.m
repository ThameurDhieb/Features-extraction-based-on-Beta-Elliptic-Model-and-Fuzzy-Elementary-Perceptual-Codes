function [ha]= angle_reyon_ellipse(xM,yM,x0,y0,teta)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% La fonction angle_reyon_ellipse permet de calculer entre [zéro et 2*pi] l'angle d'inclinaison    %%%%%
%%%%% 'ha' d'un segment (rayon d'ellipse) [O M] par rapport à la droite inclinée d'un angle teta       %%%%%
%%%%% et passant par le centre O                                                                       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 XXM = ((xM-x0)*cos(teta)) + ((yM-y0)*sin(teta));
 YYM = ((yM-y0)*cos(teta)) + ((x0-xM)*sin(teta));

 if XXM == 0
     hai = pi/2;
 else
     hai = atan(abs(YYM)/abs(XXM));
 end
 
 if ((sign(XXM) == 1)|(sign(XXM) == 0))&((sign(YYM) == 1)|(sign(YYM) == 0))
     ha = hai;
 elseif (sign(XXM) == -1)&(sign(YYM) == 1)
     ha = pi - hai;
 elseif ((sign(XXM) == -1)|(sign(XXM) == 0))&((sign(YYM) == -1)|(sign(YYM) == 0))
     ha = hai + pi;
 elseif (sign(XXM) == 1)&(sign(YYM) == -1)
     ha = (2*pi) - hai;
 else
     ha = 0;
 end
 
