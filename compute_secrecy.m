function R = compute_secrecy(rho_in, B_in, As_in, N, M, ...
    s_u, s_e, sigma2, beta0, d_ant, lambda_c, D)
    R = 0;
       for n = 1:N
        du2 = norm(rho_in(:,n)-s_u)^2 + D^2;
        de2 = norm(rho_in(:,n)-s_e)^2 + D^2;

        a_u = exp(1j*2*pi*d_ant/lambda_c*(0:M-1).'* ...
              (norm(rho_in(:,n)-s_u)/sqrt(du2)));
        h_u = sqrt(beta0/du2) * a_u;

        a_e = exp(1j*2*pi*d_ant/lambda_c*(0:M-1).'* ...
              (norm(rho_in(:,n)-s_e)/sqrt(de2)));
        h_e = sqrt(beta0/de2) * a_e;

        sinr_u = real(h_u'*B_in{n}*h_u) / (real(h_u'*As_in{n}*h_u) + sigma2);
        sinr_e = real(h_e'*B_in{n}*h_e) / (real(h_e'*As_in{n}*h_e) + sigma2);
        R = R + max(0, log2(1+sinr_u) - log2(1+sinr_e));
    end
    R = R / N;
end

% function R = compute_secrecy(rho, B, As, N, M, ...
%     s_u, s_e, sigma2, beta0, d_ant, lambda_c, D)
% 
% R = 0;
% 
% for n = 1:N
% 
%     dx_u = rho(1,n) - s_u(1);
%     dy_u = rho(2,n) - s_u(2);
%     du2 = dx_u^2 + dy_u^2 + D^2;
% 
%     dx_e = rho(1,n) - s_e(1);
%     dy_e = rho(2,n) - s_e(2);
%     de2 = dx_e^2 + dy_e^2 + D^2;
% 
%     % Correct steering
%     sin_u = dx_u / sqrt(du2);
%     sin_e = dx_e / sqrt(de2);
% 
%     a_u = exp(1j*2*pi*d_ant/lambda_c*(0:M-1).' * sin_u);
%     a_e = exp(1j*2*pi*d_ant/lambda_c*(0:M-1).' * sin_e);
% 
%     h_u = sqrt(beta0/du2) * a_u;
%     h_e = sqrt(beta0/de2) * a_e;
% 
%     sinr_u = real(h_u'*B{n}*h_u) / (real(h_u'*As{n}*h_u) + sigma2);
%     sinr_e = real(h_e'*B{n}*h_e) / (real(h_e'*As{n}*h_e) + sigma2);
% 
%     R = R + max(0, log2(1+sinr_u) - log2(1+sinr_e));
% end
% 
% R = R / N;
% end