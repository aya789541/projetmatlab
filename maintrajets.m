%aya ben aicha maintrajets
clear;
close all;
clc;
%% Initialisation des paramètres
Fe = 4e6 ; %  la fréquence échantillonnage 
Modu = 4; %  les nombre des symboles pour la modulation
roll=0.35; %roll off coefficient
nb_symbole = log2 ( Modu ) ; % Nombre de bits par symboles
DS=1e6; % le  Debit des  symboles
Ts=1/DS; % détermine le temps symbole
nombre_symbole=5000; %Nombre de symboles a emettre par paquet
NFFT=512 ;%nombre de points pour calculer Le DFT
Span=8 ;%(Ts)
alpha=1; %le coefficient de multiplication au niveau de canal 
nb_symbole = log2(Modu); %le  nombre  de symbole  qui  doit  etre  génerer en fonction  de  modulation
Te = 1/Fe; % Temps  d'echantillonage
Fse = Ts/Te;% c'est le  facteur correspondant au  surechantillonnage
TP = 4*Ts ;
Nombre_bit = nombre_symbole*nb_symbole; % Nombre de bits
l = 21 ;%nombre d'échatillons
d = 1.3;
q=0:l-1;
%% le Filtre (RRC)
filt_g_m = rcosdesign(roll, 2*TP/Ts, Fse,'sqrt');
  %% Émetteur Génération des symboles i.i.d. BPSK
        symbole_generer = zeros(1,nombre_symbole);
        symbole = randi([0,1],1,Nombre_bit);
        symboles = exp(1i*(pi/Modu+(0:Modu-1)*pi/2));

for s=1:nb_symbole:Nombre_bit
        if symbole(s) == 0 && symbole(s+1) == 0
            symbole_generer(ceil(s/nb_symbole)) = symboles(1);
        elseif symbole(s) == 0 && symbole(s+1) == 1
            symbole_generer(ceil(s/nb_symbole)) = symboles(2);
        elseif symbole(s) == 1 && symbole(s+1) == 1
            symbole_generer(ceil(s/nb_symbole)) = symboles(3);
        else 
            symbole_generer(ceil(s/nb_symbole))= symboles(4);
        end
end
   
%%surechantilloner l'ensemble  des symboles avec FSE
Symbole_surechantilloner = upsample(symbole_generer,Fse);
S_obtenue = conv(Symbole_surechantilloner,filt_g_m,"same");%signal  obtenu aprés passage  par  le  filtre  g de mse
   

 
 
        %% Canal à 1 trajet 
       h=canal_c(d,l) ;%filtre de canal 
       Symbole_canal = conv(S_obtenue, h, "same");%passage des symboles par le canal h
    
% Vérification du décalage temporel pour d = -2.5 et d = 1.3
delay = round((10 + d) / Te); % Calcul du décalage temporel arrondi à l'entier le plus proche
Symbole_canal_shift = circshift(Symbole_canal, delay); % Décalage temporel du signal
% Comparaison des signaux décalé et original
if isequal(Symbole_canal_shift, Symbole_canal)
    disp("Le décalage temporel est correctement vérifié");
else
    disp("Le décalage temporel n'est pas correctement vérifié");
end
  fvtool(Symbole_canal);
% Vérification de la réponse en fréquence du filtre

  %% Canal à 2 trajet 

%% Canal à 2 trajets
a0 = 1;
a1 = 1;
d0 = 0;
d1 = 1;
h=canal_c_2trajets(a0,a1,l,d1,d0);
Symbole_canal_1 = conv(S_obtenue, h, "same");%passage des symboles par le canal h
figure;
subplot(3,1,1);
plot(q, real(h));
title('h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 1');
xlabel('n');
ylabel('Re(h[n])');
subplot(3,1,2);
plot(q, imag(h));
title('h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 1');
xlabel('n');
ylabel('Im(h[n])');

subplot(3,1,3);
plot(q, abs(h));
title('Module de h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 1');
xlabel('n');
ylabel('|h[n]|');

a0 = 1;
a1 = 1;
d0 = 0;
d1 = 2;
h=canal_c_2trajets(a0,a1,l,d1,d0);
figure;
subplot(3,1,1);
plot(q, real(h));
title('h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 2');
xlabel('n');
ylabel('Re(h[n])');

subplot(3,1,2);
plot(q, imag(h));
title('h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 2');
xlabel('n');
ylabel('Im(h[n])');

subplot(3,1,3);
plot(q, abs(h));
title('Module de h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 2');
xlabel('n');
ylabel('|h[n]|');
a0 = 1;
a1 = 1;
d0 = 0;
d1 = 3;
h=canal_c_2trajets(a0,a1,l,d1,d0);
figure;
subplot(3,1,1);
plot(q, real(h));
title('h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 3');
xlabel('n');
ylabel('Re(h[n])');

subplot(3,1,2);
plot(q, imag(h));
title('h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 3');
xlabel('n');
ylabel('Im(h[n])');

subplot(3,1,3);
plot(q, abs(h));
title('Module de h[n] avec a0 = 1, a1 = 1, d0 = 0, d1 = 3');
xlabel('n');
ylabel('|h[n]|');
