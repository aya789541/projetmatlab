function result_ed_mmse = Filter_mmse_ed(N,K,T)
V_normaliser = zeros(N+K,1);
for S=1:N+K-1
    result_ed_mmse = zeros(N+K-1,1);
    result_ed_mmse(S,1) = 1;
    V_normaliser(S,1) = (result_ed_mmse'*(T'*((inv(T*T')+eye(N))*T)))*result_ed_mmse;
end
RETARD_mm = find(V_normaliser == max(V_normaliser));
result_ed_mmse = zeros(N+K-1,1);
result_ed_mmse(RETARD_mm+1,1) = 1;
end
