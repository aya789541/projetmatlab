function result_ed_zf = Filter_ZF_ed(S,J,H,W)
   V_normaliser = zeros(S+J-1,1);
for m=1:S+J-1
    result_ed_zf = zeros(S+J-1,1);
    result_ed_zf(m,1) = 1;
    V_normaliser(m,1) = norm((H'*W)*result_ed_zf, 2).^2;
end
retard = find(V_normaliser == max(V_normaliser));
result_ed_zf = zeros(S+J-1,1);
if (retard == S+J-1)
    result_ed_zf(retard, 1) = 1;
else
    result_ed_zf(retard+1, 1) = 1;
end
end
