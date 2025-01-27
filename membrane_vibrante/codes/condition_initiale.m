function [U0, U1] = condition_initiale(x,y,A,r0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluation des conditions initiales.
%
%   [U0, U1] = condition_initiale(x,y,A,r0)
%       
% entrées : x,y : vecteur contenant les coordonnées (x,y) de chacun des
%                 points du maillage.
%           A  :  amplitude de la fonction cosinus en entrée
%           r0 :  rayon sur lequel la fonction est non nulle.
% sorties : U0 : valeurs du déplacement prescrit en chacun des points du
%                maillage
%           U1 : condition initiale en vitesse (que l'on prendra nulle).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%condition initiale en déplacement
U0=(0.5*A*(1+cos(pi*sqrt(x.^2+y.^2)/r0))).*(sqrt(x.^2+y.^2) <= r0);

%condition initiale en vitesse
U1 = zeros(size(x));

