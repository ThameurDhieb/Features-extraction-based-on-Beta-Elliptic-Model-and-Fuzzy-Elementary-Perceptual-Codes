function [Xpoint] = filtre_lineaire(rayon, sigma_p, points)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cette fonction permet d'appliquer un filtrage lineaire sur le vecteur       %%%%%
%%%%% de données 'points' en coulissant une fenêtre gaussiènne de rayon 'rayon'   %%%%%
%%%%% d'écart type 'sigma_p'. les données filtées sont retourner dans 'Xpoint'    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


poids_courant = 1/(sigma_p*2.5066);
poids = [poids_courant];
somme_poids = poids_courant;

for i = 2 : rayon
    poids_p = poids_courant * exp(- ((1 - i)/sigma_p)^2 / 2);
    poids = [ poids; poids_p];%exp(- i^2 / 2);
    somme_poids = somme_poids + (poids_p*2);
end

Xpoint = [points(1)];

for jk = 2 : length(points)-1
    Xpointt = points(jk) * poids(1);
    for i = 2 : rayon
        indice_voisin_gch = jk - i + 1;
        if indice_voisin_gch < 1
           indice_voisin_gch = 1;
        end
        indice_voisin_drt = jk + i - 1;
        if indice_voisin_drt > length(points)
           indice_voisin_drt = length(points);
        end
        
        Xpointt = Xpointt + ( points(indice_voisin_gch)* poids(i) ) + ( points(indice_voisin_drt)* poids(i) );
    end
    Xpointt = Xpointt / somme_poids;
    Xpoint = [Xpoint; Xpointt];
end

Xpoint = [Xpoint; points(length(points),1)];

