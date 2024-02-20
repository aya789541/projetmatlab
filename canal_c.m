function h = canal_c(d,L)
fenetre_hann = hann(L)' ;
j = 0:L-1 ;
h = (sin(pi*(j - 10 - d))./(pi*(j - 10 - d))).*fenetre_hann ;%filtre

end
