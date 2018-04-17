function [a,b,x0,y0,teta,drap,erreur_MAX,ordre_ok,HA] = parametres_quart_ellipse(M_2points,M_2KMH,teta,num_fig,RRR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cette fonction permet de calculer les paramètres de l'équation du quart d'ellipse modélisant  %%%%%
%%%%%%%%%  un segment de la trajectoire dont le grand axe ait la direction de la tangente au   %%%%%%%%%%
%%%%%%%%%                 point limite correspondant au maximum de vitesse                      %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


erreur_MAX = 0;
ordre_ok = 0;

xMG = M_2points(1);
yMG = M_2points(2);
xMP = M_2points(3);
yMP = M_2points(4);


if teta == (pi / 2)
   teta = (pi / 2) - (1.0000e-10);
elseif teta == - (pi / 2)
   teta = -(pi / 2) + (1.0000e-10);
elseif teta == 0
   teta = 0 - (1.0000e-10);
elseif teta == - pi
   teta = (-pi) + (1.0000e-10);
end


alfa1 = teta + (pi/2);
alfa2 = teta;


AA1 = tan( alfa1 );
CC1 = yMG - (AA1*xMG);

AA2 = tan( alfa2 );
CC2 = yMP - (AA2*xMP);

AA = [AA1 -1;
     AA2 -1];

% K = [x0; y0];

BB = [-CC1; -CC2];

K = inv(AA) * BB;

x0 = K(1);
y0 = K(2);

a = sqrt( ((xMP - x0)^2) + ((yMP - y0)^2) );
b = sqrt( ((xMG - x0)^2) + ((yMG - y0)^2) );

% if a == max(a, b)
%    teta = teta;
%    a = a; b = b;
% end
% if b == max(a, b)
%    teta = teta + (pi/2);
%    c = a; a = b; b = c;
% end


%HBH% [ok,erreur_MAX,ordre_ok] = verification_des_fuites7(a,b,x0,y0,teta,M_5points,M_5KMH); 
%HBH% if (ok == 1)
            drap = 1;            
            
            [haG]= angle_reyon_ellipse(xMG,yMG,x0,y0,teta);
            [haP]= angle_reyon_ellipse(xMP,yMP,x0,y0,teta);
            HA = [haG, haP];

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
