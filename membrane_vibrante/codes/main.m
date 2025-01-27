%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%                   MS 204                    %%%%%%%%%%%%%%%                   
%%%%%%%%%%%%%%%                     PC6                     %%%%%%%%%%%%%%%
%%%%%%%%    Elements finis et schema de Newmark sous Matlab        %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% constantes physiques
sigma = 1;
T = 1;

%Prise en main des maillages
[Nbpt05,Nbtri05,Coorneu05,Refneu05,Numtri05,Reftri05]=lecture_msh('membraneh05.msh');
[Nbpt02,Nbtri02,Coorneu02,Refneu02,Numtri02,Reftri02]=lecture_msh('membraneh02.msh');
[Nbpt01,Nbtri01,Coorneu01,Refneu01,Numtri01,Reftri01]=lecture_msh('membraneh01.msh');

%Représentation sous matlab à l aide de la fonction trimesh
figure(1)
trimesh(Numtri05,Coorneu05(:,1),Coorneu05(:,2));
figure(2)
trimesh(Numtri02,Coorneu02(:,1),Coorneu02(:,2));
figure(3)
trimesh(Numtri01,Coorneu01(:,1),Coorneu01(:,2));

% Nbpt correspond au nombre de noeuds total N du maillage
% Pour connaitre le nb de points intérieurs il faut utiliser les références
% liées à chaque noeud. Le maillage a été défini de telle sorte que
% les points intérieurs aient la référence 0 et les points de la bordure
% ont pour référence 1.
REFINT05 = find(Refneu05==0);%permet de trouver les indices des noeuds intérieurs
REFINT02 = find(Refneu02==0);%permet de trouver les indices des noeuds intérieurs
REFINT01 = find(Refneu01==0);%permet de trouver les indices des noeuds intérieurs
%la taille de chacun de ces vecteurs (par ex length(REFINT05))
%permet donc de trouver rapidement le nb de noeuds intérieurs N_0 demandés.

% --------------------------------------------------------------
% calcul des matrices de masse et de raideur pour chaque maillage
% --------------------------------------------------------------

[M05,K05]=calculKM('membraneh05.msh',T,sigma);
[M02,K02]=calculKM('membraneh02.msh',T,sigma);
[M01,K01]=calculKM('membraneh01.msh',T,sigma);

% --------------------------------------------------------------
% Fréquences propres et déformées modales
% --------------------------------------------------------------

Nvp=45;%nombre de valeurs propres à calculer
[OM05,V05]=calculVP('membraneh05.msh',Nvp,T,sigma);
[OM02,V02]=calculVP('membraneh02.msh',Nvp,T,sigma);
[OM01,V01]=calculVP('membraneh01.msh',Nvp,T,sigma);

%Représentation de qq deformees modales
%exemple 1 : premier mode; maillage h=05
figure(4)
trimesh(Numtri05,Coorneu05(:,1),Coorneu05(:,2),V05(:,1));
%exemple 2 : 14e mode; maillage h=02
figure(5)
trimesh(Numtri02,Coorneu02(:,1),Coorneu02(:,2),V02(:,14));

%Convergence des valeurs propres selon le raffinement du maillage
figure(6)
plot(OM05,'ks')
hold on
plot(OM02,'ro')
plot(OM01,'bd')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Maillage utilise pour le calcul a partir d ici : h=0.2
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri]=lecture_msh('membraneh02.msh');
[M,K]=calculKM('membraneh02.msh',T,sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --------------------------------------------------------------
% Condition initiale
% --------------------------------------------------------------
A=2;
r0=0.5;
[U0,V0] = condition_initiale(Coorneu(:,1),Coorneu(:,2),A,r0);

%Représenter la condition initiale pour vérifier la fonction
trimesh(Numtri,Coorneu(:,1),Coorneu(:,2),U0)

% --------------------------------------------------------------
% Intégration en temps - schéma de Newmark
% --------------------------------------------------------------

%Constantes pour l intégration temporelle
delta_t=0.01;%pas de temps
Tps_final=10;%temps total de la simulation
Tps_initial = 0;
alpha = 1/delta_t;%fréquence d échantillonage
N_t = (Tps_final-Tps_initial)/delta_t; % le nombre d iterations necessaires


% Matrice Amortissement EF
% -------------------------
a_1 =0;
a_2 =0;
C = a_1*M+a_2*K;

%Allocation memoire pour les sorties
U=zeros(Nbpt,N_t);%deplacement
V=zeros(Nbpt,N_t);%vitesse
A=zeros(Nbpt,N_t);%acceleration

%Constantes du schéma de Newmark
gamma = 0.5;
beta =1/4;

%Initialisation
% on initialise U et V a partir des conditions initiales
% calculées précédemment.
U(:,1)=U0;
V(:,1)=V0;

%Tous les calculs suivants n ont lieu que sur les points interieurs
REFINT = find(Refneu==0);

%calcul de l acceleration initiale
A(REFINT,1)=M\(-C*V0(REFINT) - K*U0(REFINT));
%calcul de la matrice S (intermediaire de calcul)
S = M + gamma*delta_t*C + beta*delta_t^2*K;

%Boucle en temps
for k=2:N_t
    %Prediction
    Vstar = V(REFINT,k-1) + (1-gamma)*delta_t*A(REFINT,k-1);
    Ustar = U(REFINT,k-1) + delta_t*V(REFINT,k-1) + (.5-beta)*delta_t^2*A(REFINT,k-1);
    %Calcul de l acceleration   
    A(REFINT,k) = S\(-C*Vstar - K*Ustar);
    %Correction
    V(REFINT,k) = Vstar + delta_t*gamma*A(REFINT,k);
    U(REFINT,k) = Ustar + delta_t^2*beta*A(REFINT,k);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Animation de la solution
figure(7)
for uu=1:499
trisurf(Numtri,Coorneu(:,1),Coorneu(:,2),U(:,uu))
axis([-2 2 -2 2 -2 2])
view(-20,30)
pause(0.05)
end