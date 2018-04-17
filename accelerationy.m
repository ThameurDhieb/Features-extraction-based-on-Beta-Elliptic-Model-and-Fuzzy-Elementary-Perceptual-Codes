% la fonction accelerationx calcule la composante en acceleration suivant l'axe des Y 
% les paramètres d'entrée sont le signal de vitesse suivant l'axe des Y et la composante temps 



function AY=accelerationy(VY,t)


%-----------differentiel du temps

dt=diff(t);
%-----------differentiel du signal vitesse suivant Y

dDVY=diff(VY);


AY=(dDVY)./(dt);
AY=[0;0; AY(2:length(AY))];

return;
