function [DAT, memDAT, type_point, Ttemps] = Direction_Angulaire_continue_de_la_Tangente_DAT_Abd_AL_KARIM(indice_deb, indice_fin, points)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Angle d'inclinaison de la tangente au tracé, DAT. %%%%%%%%%%%%

k1 = indice_deb;
k2 = indice_fin;

DAT = [];
type_point = [];
Ttemps = [];
u = 1;%4%1%4%2; %2;
dk = (2*u) + 1;
Lp = length(points);
pas = 0.01;

cotion_preced = 1;

for k = k1 : k2,
    if (k - u < 1)&(k + u  < Lp)
        if (points(k + u,1) - points(1,1)) == 0
           denum= 1.0e-7;
        else
           denum = (points(k + u,1) - points(1,1));
        end
%        cotion = (points(k + u,2) - points(1,2)) / (points(k + u,1) - points(1,1));
        cotion = (points(k + u,2) - points(1,2)) / denum;
        dat_k = atan( cotion );
    elseif (k + u > Lp)
%        cotion = (points(Lp,2) - points(k - u,2)) / (points(Lp,1) - points(k - u,1));
        if (points(Lp,1) - points(k-1 ,1)) == 0
           denum = 1.0e-7;
        else
           denum = (points(Lp,1) - points(k-1 ,1));
        end
%        cotion = (points(Lp,2) - points(k-1 ,2)) / (points(Lp,1) - points(k-1 ,1));
        cotion = (points(Lp,2) - points(k-1 ,2)) / denum;
        
        dat_k = atan( cotion );        
    elseif (k + u < Lp)&(k - u > 1)
        num = (points(k + u,2) - points(k - u,2));
        denum = (points(k + u,1) - points(k - u,1));
%         if denum == 0
%            denum = 1/100;%10000;
%         end
        if denum == 0
           denum = 1.0e-7;
           cotion = abs( num / denum ) * sign(cotion_preced);
        else
           cotion = num / denum;
        end
        cotion_preced = cotion;

        dat_k = atan( cotion );
        
    end

    DAT = [DAT; dat_k];
%    type_point = [type_point, 0];
    Ttemps = [Ttemps,k*pas];
    
    
end
memDAT = DAT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eviter les discontinuités dans la variation de l'angle d'inclinaison de la tangente DAT

tour_ajout = 0;
L_DAT = length(DAT);

for i = 2 : L_DAT,
    tour_ajout = fix((DAT(i) - DAT(i-1))/(pi))*pi;
    DAT(i) = DAT(i) - tour_ajout;

    diff_angle_1 = abs( (DAT(i)) - DAT(i-1) );
    diff_angle_2 = abs( (DAT(i) + pi) - DAT(i-1) );
    diff_angle_3 = abs( (DAT(i) - pi) - DAT(i-1) );
    
    x = [diff_angle_1, diff_angle_2, diff_angle_3];
    diff_angle_min = min(x);
    if diff_angle_min == diff_angle_1
        DAT(i) = DAT(i);
        hu = 0;
    elseif diff_angle_min == diff_angle_2
        DAT(i) = DAT(i) + pi;
        hu = pi;
    elseif diff_angle_min == diff_angle_3
        DAT(i) = DAT(i) - pi;
        hu = -pi;
    end
end
memDAT2 = DAT;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% recherche des minimums locaux de rayon de courbure %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Diff_DAT = [];
for i = 2 : L_DAT,
    Diff_DAT = [Diff_DAT; DAT(i) - DAT(i-1)];
end
Diff_DAT = [0 ;Diff_DAT];
for i = 2 : L_DAT-1,
    if ( Diff_DAT(i) < Diff_DAT(i-1) )&( Diff_DAT(i) < Diff_DAT(i+1) )
       type_pnt = 1;                                       % minimum de vitesse ou de rayon de courbure
    else
       type_pnt = 0;
    end
    type_point = [type_point; type_pnt];
end
type_point = [0;type_point;0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sens_preced = 0;

for i = 2 : L_DAT,
    
    sens_actuel = sign( DAT(i) - DAT(i-1) );
    if (sens_actuel ~= sens_preced)&(sens_preced ~= 0)           %sign( DAT(i-1) - DAT(i-2) )
       if (type_point(i) == 1)|(type_point(i-1) == 1)            % minimum de vitesse ou de rayon de courbure
          indice_k = i + k1 - 1;
          [continuite,inversion] = verification_changement_sens_parcour_Abd_AL_KARIM(indice_k,points);
          if (continuite == 0)&(inversion == 1)
             DAT(i) = DAT(i) + (sens_preced*(pi));
             sens_actuel = sign( DAT(i) - DAT(i-1) );
          else
             DAT(i) = DAT(i);
          end          
       end
    end
    sens_preced = sens_actuel;

end





