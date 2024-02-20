%% aya ben aicha
clear ;
close all ;
clc ;

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
h=[1 zeros(1,3) alpha];% filtre de canal de  legalisation
P =100;% Ordre des filtres ZF/MMSE
lambda = 1 ;%le facteur  qui sera  multiplier  avec  le filtre adapté
%% le Filtre (RRC)
filt_g_m = rcosdesign(roll, 2*TP/Ts, Fse,'sqrt');

%% Paramétres  qui permet  de  determiner la probabilité dlerreur binaire en fonction
%du rapport signal à bruit 
sig_variance=1 ;%  determine la variance  des symboles
energie_mse=sum(filt_g_m.*filt_g_m); % determine Energie  de MSE
%on determine la liste des energieB/N0 
liste_eb_nb = 0:0.5:20; 
EBno= 10.^( liste_eb_nb /20) ;
Sigma= sig_variance * energie_mse ./ ( nb_symbole * EBno ) ; % Variance du bruitcomplexe en bande de base
Taux_EB = zeros(size(EBno)); %  % résultats de taux deuureur binare
Probeurreur = qfunc ( sqrt ( 2*EBno ) ) ; % Probabilités d eurreurs





%%
for N=1:length(EBno)
     compteur_b = 0;%permet de compter les bits deurreur 
    compteur_e = 0;%permet de compter le nobre deuureur entre les symbole emis/recu
 
    while compteur_e < 100
 
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
      
 
        %% Canal  sur  lequel on va  effectuer  l'églisation
       Symbole_canal = conv(S_obtenue, h, "same");%passage des symboles par le canal h

        %% Récepteur   
gadapte = lambda*(conj(flip(filt_g_m))) ;%filtre adapté 
v_r = conv(Symbole_canal, gadapte,'same');
v_t = conv(filt_g_m,conv(gadapte,h,'same'),'same') ; %filtre equivalent à laquelle on va effectuer lechantillonage au rythme  NTS
Vn = v_t(1:Fse:length(v_t)) ;%echantillonage de vt au rythme nts
rn = v_r(1+2*TP/Te:Fse:length(v_r));% Echantillonage au rythme Ts de signal reçu
 
%% Matrice de convolution T du canal de taille P*(P+K-1)
k = length(Vn);
T = zeros(P, P+k-1);

for N = 1:P
    for j = 1:k
        T(N,j+N-1) = Vn(j);
    end
end

%%Construction des égaliseurs zf/mmse
symbole_apres_zf = [] ;%contient le résultat des symbole aprés avoir appliquer le filtrage zéroforcing
symbole_apres_MMSE = [] ;%contient le résultat des symbole aprés avoir appliquer le filtrage MMSE
% Egaliseur ZF--------------------------------------
W = pinv(T') ;
ed_egalisation = Filter_ZF_ed(P,k,T,W);
filtre_zero_forcing = W*ed_egalisation;
% Egalisation  final des symboles ZF
longeur=P+k-1;
for comp_zf=longeur:length(rn)
    symbole_apres_zf = [symbole_apres_zf filtre_zero_forcing'*rn(comp_zf:-1:comp_zf-P+1)'];  
end
% Egaliseur  MMSE---------------------------------
ed_egalisation=Filter_mmse_ed(P,k,T);
filter_MMSE = ((inv(T*T')+eye(P))*T)*ed_egalisation;

% Egalisation  final des symboles MMSE
for comp_mmse=longeur:length(rn)
    symbole_apres_MMSE = [symbole_apres_MMSE filter_MMSE'*rn(comp_mmse:-1:comp_mmse-P+1)'];
 
end

err=0;
for k=1:length(symbole_apres_zf)
   if(symbole_generer(k)~=symbole_apres_zf(k)) 
    err=err+1;%l'erreur binaire 
   end     
end
        
        compteur_e=compteur_e+err; %incrementer le compteur d'erreur
        compteur_b=compteur_b+2*nombre_symbole; %incrementer le nombre de bit envoyee

    end
    Taux_EB(N) = compteur_e/compteur_b ;
end


%% Graphes

 %dsp                                              
figure,
[dsp,f]=pwelch(S_obtenue,NFFT,0,NFFT,Fe,'centered');  
semilogy(f,dsp,'r')
xlabel("Frequency (HZ) ")
ylabel("Power/Frequency (db/Hz) ")
title("DSP  pratique ")
grid on

                                            
figure();
semilogy(liste_eb_nb,Probeurreur,'b');
xlabel(' SNR ');
ylabel('pb eurreur');
title("la probabilité d’erreur binaire en fonctiondu rapport signal à bruit")
grid on;


% constellation avant/aprés l'egalisation  en ZF
figure,
subplot(1,2,1)
plot(v_r, '*') ;
title("constellation avant Égalistion ZF");

subplot(1,2,2)
plot(symbole_apres_zf, '*') ;
title("constellation aprés Egalisation ZF")

% constellation avant/aprés l'egalisation  en MMSE
figure,
subplot(1,2,1)
plot(v_r, '*') ;
title("constellation avant Égalistion  MMSE");

subplot(1,2,2)
plot(symbole_apres_MMSE, '*') ;
title("constellation aprés Egalisation MMSE")

