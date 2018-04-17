function VX =vitessex(X,t)
%calcul de la vitesse suivant X

% les entrées sont la composante suivant X et le temps t
% la sortie est la composante de vitesse suivant X, VX



%X=points(:,1);
dt=diff(t);
dDx=diff(X);

VX=(dDx)./(dt);
VX=[0; VX];
VX=[VX(1:length(VX)-1);0];

return;

