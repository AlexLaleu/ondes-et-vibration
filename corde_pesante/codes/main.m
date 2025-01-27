% Cours MEC_4MS04_TA
%
% PC 5
%
% Corde tendue soumise à la gravité : effet d une tension variable
%
% Projection des équations du mouvement sur les modes propres d une corde
% sans gravité
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Nombre de modes
nmod = 100;

%Nombre de modes a tracer
Nmax = 15 ;

% Points d espace pour lesquels seront calculés les modes
x=0:0.002:1;
Nx = length(x) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construction des matrices de masse et de raideur

% matrice de masse
M = eye(nmod);

% Matrice de raideur
K = zeros (nmod,nmod);
for n1 = 1:nmod ;
    for m1 = 1:nmod ;
        if m1 == n1
            K(m1,n1) = n1^2*pi^2 + amn(m1,n1) + bmn(m1,n1);
        else
            K(m1,n1) = amn(m1,n1) + bmn(m1,n1);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des fréquences et modes propres et tri croissant des frequences
% VC matrice contenant les vecteurs propres
% DC matrice contenant les valeurs propres dans la diagonale

[VC,DC] = eig(K,M);  % [A]*V = a*[B]*V ssi [V,a] = eig([A],[B]) ici on a w^2[M]*V = K[V]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tri des fréquences propres par ordre croissant
% et rangement des vecteurs propres correspondants
%
% FR contient à présent les fréquences (et non pas
% le carré des pulsations propres comme c est le
% cas pour VC)

[FR, ind] = sort (sqrt(diag(DC))/(2*pi)) ;
VCR = VC (:,ind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vérification de l orthogonalité des modes propres
% Calcul des raideurs et masses modales
% CAD : coeffs des matrices (diagonales) de raideur et de masse
% dans la base des modes propres.

KG = VCR'*K*VCR;
MG = VCR'*M*VCR;

eps = 1e-10;
for n1 = 1:nmod 
    for m1 = 1:nmod 
        if n1 ~= m1
            if abs(KG(n1,m1)) > eps
                %disp('Propblème sur l orthogonalité provenant de K');
            end
            if abs(MG(n1,m1)) > eps
                %disp('Problème sur l orthogonalité provenant de M');
            end
        end
    end
end

figure (1)
imagesc(KG)

figure(2)
imagesc(MG)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calcul des modes propres de la corde avec gravité -> Psi
%
% On calcule ici une matrice Psi de taille (Nx,Nmod)
% Chaque colonne représente un mode, exprimé en chaque point x
% Ce calcul se fait à partir des modes Phi de la corde sans gravité
% qu il convient aussi de calculer

Phi = zeros(Nx,nmod);
for n = 1:Nx
    for m = 1:nmod
        Phi(n,m) = sqrt(2)*sin(m*pi*x(n));
    end
end

PSI = Phi*VCR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tracé des Nmax premiers modes propres et comparaison avec les modes de corde

figure (3);
clf;
for n1 = 1: Nmax
    subplot(5,ceil(Nmax/5),n1); 
    plot (x,PSI(:,n1),'linewidth',2)
    hold on;
    plot(x,sqrt(2)* sin (n1*pi*x),'r');
    title(['mode =', num2str(n1), ', FR =', num2str(FR(n1)), 'Forig = ', num2str(n1/2)] );
    title(['mode =', num2str(n1), ', FR =', num2str(FR(n1))]);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation temporelle
% Condition initiale : \sqrt(2)*sin(pi x), c est à dire une déformation
% selon le premier mode dans la base phi
% Le vecteur q0 a donc 1 dans son premier élément et des zéros
% sur tout le reste de la colonne

q0 = zeros(nmod,1);
q0(1)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r0 est l expression de cette condition initiale dans la base
% modale de la corde avec tension variable


r0 = VCR\q0;  % car q = [VCR]*r où r est la contribution de q décomposé sur la base modale de la corde pesante.
%  on aurait pu mettre r0 = VCR'*M*q0;     % car VCR'*M*VCR = Id et q = [VCR]*r
% ou encore r0 = inv(VCR)*q0;

%%%%%%%%%%%%%%%
% Vecteur temps
t=0:0.05:20;
Nt=length(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Déplacement dans la base modale en fonction du temps
% On utilise ici un cosinus car la condition initiale est donnée
% sur le déplacement.
% r est ici une matrice de taille [nmod, Nt]
%
% La forme générale de la solution est r_n(t) = r0_n cos(2*pi*F_n*t)
% On pourra utiliser la fonction repmat en étape intermédiaire:
% r0m = repmat(r0, 1, Nt);

r0m = repmat(r0, 1, Nt);
COS = zeros(nmod,Nt);
for i=1:nmod
    for j=1:Nt
        COS(i,j) = cos(FR(i)*2*pi*t(j));
    end
end
r = COS .* r0m ; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Déplacement "physique"
% y est ici une matrice de taille [Nx, Nt]
% 

y = zeros (Nx,Nt) ;

y= PSI*r;

%%%%%%%%%%%%%
% Animation !
%%%%%%%%%%%%%

h = figure(4);
set(gcf,'DoubleBuffer','On','BackingStore','Off');
plot([0 0], [0 1]);
set(gca,'Xlim', [-10 10]);
hlin = line(x,zeros(size(x)));
for i1=1:length(t);
    set(hlin,'XData',y(:,i1),'YData',1-x,'LineStyle','-','Marker','none','Color','r', 'LineWidth',3);
    drawnow;
    pause(0.1);
end