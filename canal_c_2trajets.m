function h = canal_c_2trajets(a0,a1,L,d1,d0)
fenetre_hann = hann(L)';
s = 0:L-1;
h = a0*(sin(pi*(s-10-d0))./(pi*(s-10-d0))).*fenetre_hann + a1*(sin(pi*(s-10-d1))./(pi*(s-10-d1))).*fenetre_hann;


end
