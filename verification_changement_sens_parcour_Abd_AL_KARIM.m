function [continuite,inversion] = verification_changement_sens_parcour_Abd_AL_KARIM(k_Mextremum,points)

continuite = 1;
inversion = 0;

%%%%%%%%%  verification du changement de sens de rotation  %%%%%%%%%
w = 3;
u = 2; %2;
L = length(points);
%k_Mextremum = liste_extremum(rong_liste_extremum1,2); 

k_M1 = k_Mextremum - w;
if (k_M1 >0)&(k_Mextremum+w - 1 < L)

   if ((k_Mextremum-w) - u < 1)
      denum = (points(2,1) - points(1,1));
      if  denum == 0
         denum = 1.0e-8;
      end
       
%      dat_k = atan( (points(2,2) - points(1,2)) / (points(2,1) - points(1,1)) );
      dat_k = atan( (points(2,2) - points(1,2)) / denum );
   elseif ((k_Mextremum-w) + u > L)
      denum = (points(L,1) - points((k_Mextremum-w) - u,2));
      if denum == 0
         denum = 1.e-8;
      end
%      dat_k = atan( (points(L,2) - points((k_Mextremum-w) - u,2)) / (points(L,1) - points((k_Mextremum-w) - u,2)) );
      dat_k = atan( (points(L,2) - points((k_Mextremum-w) - u,2)) / denum );
   else
      denum = (points((k_Mextremum-w) + u,1) - points((k_Mextremum-w) - u,1));
      if denum == 0
         denum = 1.e-8;
      end

%      dat_k = atan( (points((k_Mextremum-w) + u,2) - points((k_Mextremum-w) - u,2)) / (points((k_Mextremum-w) + u,1) - points((k_Mextremum-w) - u,1)) );
      dat_k = atan( (points((k_Mextremum-w) + u,2) - points((k_Mextremum-w) - u,2)) / denum );

      if (points((k_Mextremum-w) + u,1) - points((k_Mextremum-w) - u,1)) < 0
         dat_k = dat_k + pi;
      end
   end

   datM1 = dat_k;

   X_M1 = points(k_Mextremum-w,1);                  % M1 un point en amont du point extremum Mk
   Y_M1 = points(k_Mextremum-w,2);
   X_M2 = points(k_Mextremum-w + 1,1);              % M2 un point en amont du point extremum Mk 
   Y_M2 = points(k_Mextremum-w + 1,2);              % qui vient apres M1

   X_M3 = points(k_Mextremum,1);                    % M3 le point extremum Mk
   Y_M3 = points(k_Mextremum,2);


   X_M4 = points(k_Mextremum+w - 1,1);               % M4 un point en aval du point extremum Mk 
   Y_M4 = points(k_Mextremum+w - 1,2);

   long_segM1M2 = sqrt( ((Y_M2 - Y_M1)^2) + ((X_M2 - X_M1)^2) );
   long_segM1M3 = sqrt( ((Y_M3 - Y_M1)^2) + ((X_M3 - X_M1)^2) );
   long_segM1M4 = sqrt( ((Y_M4 - Y_M1)^2) + ((X_M4 - X_M1)^2) );

   alfat_M1M2 = atan( (Y_M2 - Y_M1)/(X_M2 - X_M1) );
   if (X_M2 - X_M1) < 0
      alfat_M1M2 = alfat_M1M2 + pi;
   end
   alfat_M1M3 = atan( (Y_M3 - Y_M1)/(X_M3 - X_M1) );
   if (X_M3 - X_M1) < 0
      alfat_M1M3 = alfat_M1M3 + pi;
   end
   alfat_M1M4 = atan( (Y_M4 - Y_M1)/(X_M4 - X_M1) );
   if (X_M4 - X_M1) < 0
      alfat_M1M4 = alfat_M1M4 + pi;
   end

   Ang_diff_M1M2 = alfat_M1M2 - datM1;
   Ang_diff_M1M3 = alfat_M1M3 - datM1;
   Ang_diff_M1M4 = alfat_M1M4 - datM1;

   X_projete_M2 = long_segM1M2 * cos(Ang_diff_M1M2);
   X_projete_M3 = long_segM1M3 * cos(Ang_diff_M1M3);
   X_projete_M4 = long_segM1M4 * cos(Ang_diff_M1M4);

   if ((X_projete_M4 > X_projete_M3)&(X_projete_M3 > X_projete_M2)&(X_projete_M2 > 0))|((X_projete_M4 < X_projete_M3)&(X_projete_M3 < X_projete_M2)&(X_projete_M2 < 0))
      continuite = 1;
      inversion = 0;
   else
      continuite = 0;
      inversion = 1;
   end

end



