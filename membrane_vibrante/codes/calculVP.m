function [OM,V]=calculVP(nom_maillage,Nvp,T,sigma)

%Calcul des matrices de masse et de raideur
[M,K]=calculKM(nom_maillage,T,sigma);

%Calcul des valeurs propres
%Attention à l'option 'smallestabs' de eigs pour commencer par la plus petite vp
[Vi,D] = eigs(K,M,Nvp,'smallestabs');

OM=sqrt(diag(D));

% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% %Rajout des zéros au bord pour représentation
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh(nom_maillage);
REFINT = find(Refneu==0);%Noeuds interieurs
V=zeros(Nbpt,Nvp);
V(REFINT,:)=Vi;


