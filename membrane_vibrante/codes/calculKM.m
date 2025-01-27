function [M,K]=calculKM(nom_maillage,T,sigma)

% ----------------------
% calcul des matrices EF
% ----------------------

[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh(nom_maillage);


% declarations
% ------------
K = sparse(Nbpt,Nbpt); % matrice de rigidite
M = sparse(Nbpt,Nbpt); % matrice de masse

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  
  % calcul des matrices elementaires du triangle l 
  
   [Kel]=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:));
           
   [Mel]=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),Coorneu(Numtri(l,3),:));
    
    % On fait l'assemblage de la matrice globale
  for i=1:3
       for j=1:3
           K(Numtri(l,i),Numtri(l,j)) = K(Numtri(l,i),Numtri(l,j)) + Kel(i,j);
           M(Numtri(l,i),Numtri(l,j)) = M(Numtri(l,i),Numtri(l,j)) + Mel(i,j);
       end
   end

end % for l


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%ELIMINATION

% On suppose que les noeuds Dirichlet ont une reference 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REFINT = find(Refneu==0);

% On elimine les lignes et les colonnes associees a des noeuds Dirichlet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = K(REFINT,REFINT);
M = M(REFINT,REFINT);

%Dimensionnement
K=T*K;
M=sigma*M;

