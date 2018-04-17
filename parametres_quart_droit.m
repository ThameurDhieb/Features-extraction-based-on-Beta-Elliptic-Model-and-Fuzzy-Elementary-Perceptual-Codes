function [a,b,x0,y0,teta] = parametres_quart_droit(xMG,yMG,xMP,yMP,teta)


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
