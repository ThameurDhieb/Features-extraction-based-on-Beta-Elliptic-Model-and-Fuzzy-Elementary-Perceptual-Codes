function [points,XYpoint] = filtre_lineaire_1(rayon, sigma_p, indice_deb, indice_fin, points)

%rayon = 4;
%rayon = 5;
%sigma_p = 1.2; %2; % 1.5;%1;%2;

XYpoint = [];

poids_courant = 1/(sigma_p*2.5066);
poids = [poids_courant];
somme_poids = poids_courant;

for i = 2 : rayon
    poids_p = poids_courant * exp(- ((1 - i)/sigma_p)^2 / 2);
    poids = [ poids; poids_p];%exp(- i^2 / 2);
    somme_poids = somme_poids + (poids_p*2);
end

for jk = indice_deb : indice_fin
    Xpointt = points(jk,1) * poids(1);
    Ypointt = points(jk,2) * poids(1);
    
    for i = 2 : rayon
        indice_voisin_gch = jk - i + 1;
        if indice_voisin_gch < 1
           indice_voisin_gch = 1;
        end
        indice_voisin_drt = jk + i - 1;
        if indice_voisin_drt > length(points)
           indice_voisin_drt = length(points);
        end
        
        Xpointt = Xpointt + ( points(indice_voisin_gch,1)* poids(i) ) + ( points(indice_voisin_drt,1)* poids(i) );
        Ypointt = Ypointt + ( points(indice_voisin_gch,2)* poids(i) ) + ( points(indice_voisin_drt,2)* poids(i) );
    end
    Xpointt = Xpointt / somme_poids;
    Ypointt = Ypointt / somme_poids;
    
    XYpoint = [XYpoint; Xpointt Ypointt];
end

for jk = indice_deb : indice_fin

    points(jk,1) = XYpoint(jk - indice_deb + 1, 1);
    points(jk,2) = XYpoint(jk - indice_deb + 1, 2);

end


% for jk = indice_deb : indice_fin
%     Xpointt = points(jk,1) * poids(1);
%     Ypointt = points(jk,2) * poids(1);
%     
%     for i = 2 : rayon
%         indice_voisin_gch = jk - i + 1;
%         if indice_voisin_gch < 1
%            indice_voisin_gch = 1;
%         end
%         indice_voisin_drt = jk + i - 1;
%         if indice_voisin_drt > length(points)
%            indice_voisin_drt = length(points);
%         end
%         
%         Xpointt = Xpointt + ( points(indice_voisin_gch,1)* poids(i) ) + ( points(indice_voisin_drt,1)* poids(i) );
%         Ypointt = Ypointt + ( points(indice_voisin_gch,2)* poids(i) ) + ( points(indice_voisin_drt,2)* poids(i) );
%     end
%     Xpointt = Xpointt / somme_poids;
%     Ypointt = Ypointt / somme_poids;
%     points(jk,1) = Xpointt;
%     points(jk,2) = Ypointt;
%     
%     XYpoint = [XYpoint; Xpointt Ypointt];
% end
