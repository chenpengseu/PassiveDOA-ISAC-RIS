function x = MUSIConesnapshot(r, param)
target_num = param.K;
r = r(:);
N = length(r);

L = 8;
hankel_mat = zeros(L, N-L+1);
for idx = 1:L
    hankel_mat(idx,:) = param.vecH(r(idx:idx+N-L));
end

[U,D,V] = svd(hankel_mat);

D = abs(diag(D));
[~, sort_idx] = sort(D,'descend');

U1 = U(:, sort_idx(1:target_num));
U2 = U(:, sort_idx(1+target_num:end));

dic_mat = param.get_steer(param.cont_ang, L);

sp = param.vec(sum(abs(dic_mat).^2,1))./param.vec(sum(abs(dic_mat'*U2).^2, 2));
sp = sp/max(sp);
x = param.vecH(sp);
% figure; plot(ang_range, 10*log10(x)); hold on; stem(theta, zeros(3,1), 'BaseValue', -100);
