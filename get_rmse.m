function rmse = get_rmse(est_spectrum, theta_TR)
[pks, pks_idx] = findpeaks(est_spectrum(:, 2));
est_ang = zeros(length(theta_TR),1);
if ~isempty(pks_idx) 
    for idx = 1:length(theta_TR)
        [~,min_idx] = min(abs(theta_TR(idx)-est_spectrum(pks_idx,1)));
        est_ang(idx) = est_spectrum(pks_idx(min_idx),1);
    end
end
rmse = sqrt(norm(est_ang-theta_TR)^2/length(theta_TR));