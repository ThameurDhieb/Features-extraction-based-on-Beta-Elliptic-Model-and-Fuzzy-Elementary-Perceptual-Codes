% la fonction accelerationx calcule la composante en acceleration suivant l'axe des X 
% les paramètres d'entrée sont le signal de vitesse suivant l'axe des X et la composante temps 


function AX=accelerationx(VX,t)

%-----------differentiel du temps

dt=diff(t);
%-----------differentiel du signal vitesse suivant X
dDVX=diff(VX);


AX=(dDVX)./(dt);
AX=[0;0; AX(2:length(AX))];

return;
