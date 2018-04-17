function VY=vitessey(Y,t)



%calcul de la vitesse suivant Y

% les entrées sont la composante suivant Y et le temps t
% la sortie est la composante de vitesse suivant Y, VY


dt=diff(t);
dDy=diff(Y);

VY=(dDy)./(dt);
VY=[0; VY];
VY=[VY(1:length(VY)-1);0];

return;
