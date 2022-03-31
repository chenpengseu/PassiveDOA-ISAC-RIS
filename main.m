clear all; close all;
param.vec = @(MAT) MAT(:);
param.vecH = @(MAT) MAT(:).';
rng(111);
param.theta_RS = 0;
param.d_E = 0.5;
param.M = 64; % the number of RIS elements
param.get_steer = @(theta, L) exp(1j*2*pi*[0:1:L-1].'*param.d_E*param.vecH(sind(theta+param.theta_RS)));
param.theta_AR = -10+rand(1);
param.K = 3;
SNR_dB = 20;
param.d_AT = 20;
param.d_TR = 30;
param.d_RS = 3;
param.d_AR = 5;
z = 1/(param.d_AT*param.d_TR*param.d_RS)*exp(1j*rand(param.K,1));
q = 1/(param.d_AR*param.d_RS)*exp(1j*rand(1));
param.theta_TR = [-25, 15, 30].'+rand(param.K,1);
param.N = 16; % the number of measurements

% generate the measurement matrix
a_tmp = param.get_steer(param.theta_AR, param.M);
A_tmp = a_tmp*a_tmp';
cvx_begin sdp quiet
	variable G_tilde(param.M,param.M) hermitian 
	minimize (trace(A_tmp*G_tilde))
	subject to
		G_tilde>=0
		for idx_m = 1:param.M
			G_tilde(idx_m, idx_m) == 1;
		end
cvx_end
[U_tilde, D_tilde] = eig(G_tilde);
D_tilde = diag(sqrt(real(diag(D_tilde))));
N_test = 1e4;
g_tilde = sqrt(1/2)*(randn(param.M,N_test)+1j*randn(param.M,N_test));
G_tmp = exp(1j*angle(U_tilde*D_tilde*g_tilde));
G = G_tmp';
param.G = G(1:param.N, :);
% test the measurement matrix
ang_mean = [-45:0.01:45].';
pow_g = abs(G*param.get_steer(ang_mean, param.M)).^2;
pow_g = mean(pow_g);
pow_g = pow_g/max(pow_g);
figure; 
plot(ang_mean, 10*log10(pow_g), 'LineWidth', 3);
hold on; stem(param.theta_AR, 0,'o', 'BaseValue', -8, 'LineWidth', 3, 'MarkerSize', 12);
xlabel('Spatial angle(deg)');
ylabel('Beamforming power (dB)');
legend('Optimized measurement matrix', 'AP direction');
grid on;

% add noise
r_q = param.G*param.get_steer(param.theta_AR, param.M)*q;
r_z = param.G*param.get_steer(param.theta_TR, param.M)*z;
noise_pow = norm(r_z)^2/10^(SNR_dB/10);
w = sqrt(1/2)*(randn(size(r_z))+1j*randn(size(r_z)));
noise_var = noise_pow/norm(w)^2;
w = sqrt(noise_var)*w;
r = r_z+r_q+w;

param.dic_ang = [-45:1:45].'; % for grid estimation
param.dic_mat = param.get_steer(param.dic_ang, param.M);

param.cont_ang = [-45:0.001:45].'; % for off-grid estimation
param.cont_dic  = param.get_steer(param.cont_ang, param.M);

% proposed method
b = param.G*param.get_steer(param.theta_AR, param.M);
rho =  sqrt(log(param.M)*param.M)*sqrt(noise_var);
cvx_begin sdp quiet
	variable est_xi(param.M, 1) complex
	variable u(param.M, 1) complex
	variable est_eta complex
	variable Z(param.M,param.M) hermitian toeplitz
	variable nu(1,1)
	minimize(quad_form(r-param.G*est_xi-est_eta*b,eye(length(r)))+rho/2*(nu+1/param.M*trace(Z)))
	subject to
	[Z, est_xi; est_xi', nu]>=0
	Z(:, 1) == u
cvx_end
est_x = MUSIConesnapshot(est_xi, param);
est_spectrum = [param.cont_ang, zeros(length(param.cont_ang), 1)];
est_spectrum(:,2) = abs(est_x);
rmse_propose = get_rmse(est_spectrum, param.theta_TR)

sp_tmp = 20*log10(abs(est_x));
sp_tmp = sp_tmp - max(sp_tmp);
figure;
plot(param.cont_ang, sp_tmp, 'LineWidth', 3);
legend_str{1} = 'Proposed method';
sp_tmp2 = 20*log10(abs(z));
sp_tmp2 = sp_tmp2 - max(sp_tmp2);
hold on; stem(param.theta_TR, sp_tmp2, 'BaseValue', min(sp_tmp)-10, 'LineWidth', 3, 'MarkerSize',12,'Marker', 'o');
legend_str{2} = 'Ground-truth DOA';
legend(legend_str);
axis([-45,45,min(sp_tmp)-10,10]);
grid on;
    
% CRLB
B=zeros(param.M,param.K);
for idx = 1:param.K
	B(:, idx) = 1j*2*pi*param.d_E*z(idx)*cosd(param.theta_TR(idx)+param.theta_RS)*(param.get_steer(param.theta_TR(idx), param.M).*[0:param.M-1].');
end
F = length(r_z)/noise_pow*[2*real(B'*param.G'*param.G*B), B'*param.G'*param.G*param.get_steer(param.theta_TR, param.M), B'*param.G'*param.G*param.get_steer(param.theta_AR, param.M);
	param.get_steer(param.theta_TR, param.M)'*param.G'*param.G*B, param.get_steer(param.theta_TR, param.M)'*param.G'*param.G*param.get_steer(param.theta_TR, param.M), param.get_steer(param.theta_TR, param.M)'*param.G'*param.G*param.get_steer(param.theta_AR, param.M);
	param.get_steer(param.theta_AR, param.M)'*param.G'*param.G*B, param.get_steer(param.theta_AR, param.M)'*param.G'*param.G*param.get_steer(param.theta_TR, param.M), param.get_steer(param.theta_AR, param.M)'*param.G'*param.G*param.get_steer(param.theta_AR, param.M)];
crlb_all = abs(diag(inv(F)));
crlb = rad2deg(sqrt(sum(crlb_all(1:param.K))/param.K)) 
