function avg_sec_rate = run_AO_scheme(M, D, N, Pmax, sigma2, beta0, ...
    d_ant, lambda_c, Vmax, Gamma_t, Gamma_e, ...
    s_u, s_t, s_e, rho_I, rho_F, scheme_type)

    steer = @(rho_uav, s_gnd) ...
        exp(1j*2*pi*d_ant/lambda_c*(0:M-1).'* ...
        (norm(rho_uav-s_gnd)/sqrt(norm(rho_uav-s_gnd)^2+D^2)));
    chan = @(rho_uav, s_gnd) ...
        sqrt(beta0/(norm(rho_uav-s_gnd)^2+D^2))*steer(rho_uav, s_gnd);

    % Init trajectory
    rho = zeros(2,N);
    for n = 1:N
        rho(:,n) = rho_I + (rho_F-rho_I)*n/(N+1);
    end

    B_s = cell(1,N); As_s = cell(1,N);
    for n = 1:N
        B_s{n} = (Pmax/(2*M))*eye(M);
        As_s{n} = (Pmax/(2*M))*eye(M);
    end

    for iter = 1:3  % Few AO iters for speed
        %% Beamforming
        Bl = B_s; Al = As_s;
        for s1 = 1:3
            for n = 1:N
                hu = chan(rho(:,n),s_u); he = chan(rho(:,n),s_e);
                at = steer(rho(:,n),s_t); ae = steer(rho(:,n),s_e);
                dt2 = norm(rho(:,n)-s_t)^2+D^2;
                de2 = norm(rho(:,n)-s_e)^2+D^2;
                Hu = hu*hu'; He = he*he';
                
                nul = real(trace(Hu*Bl{n})); dul = real(trace(Hu*Al{n}))+sigma2;
                nel = real(trace(He*Bl{n})); del = real(trace(He*Al{n}))+sigma2;
                
                c1u=1/(log(2)*(nul+dul)); c2u=nul/(log(2)*dul*(nul+dul));
                c1e=1/(log(2)*(nel+del)); c2e=nel/(log(2)*del*(nel+del));
                Rcu=log2(1+nul/dul)-c1u*nul+c2u*real(trace(Hu*Al{n}));
                Rce=log2(1+nel/del)-c1e*nel+c2e*real(trace(He*Al{n}));
                
                cvx_begin sdp quiet
                    cvx_solver mosek
                    variable Bv(M,M) hermitian semidefinite
                    variable Av(M,M) hermitian semidefinite
                    Rus=Rcu+c1u*real(trace(Hu*Bv))-c2u*real(trace(Hu*Av));
                    Res=Rce+c1e*real(trace(He*Bv))-c2e*real(trace(He*Av));
                    maximize(Rus-Res)
                    subject to
                        real(trace(Bv))+real(trace(Av))<=Pmax;
                        real(at'*(Bv+Av)*at)>=Gamma_t*dt2;
                        if Gamma_e < Inf
                            real(ae'*(Bv+Av)*ae)<=Gamma_e*de2;
                        end
                cvx_end
                if strcmp(cvx_status,'Solved')||strcmp(cvx_status,'Inaccurate/Solved')
                    B_s{n}=Bv; As_s{n}=Av;
                end
            end
            Bl=B_s; Al=As_s;
        end
        % Rank-1
        for n=1:N
            [V,Dg]=eig(B_s{n}); dv=real(diag(Dg));
            [mv,idx]=max(dv);
            if mv>0; bn=sqrt(mv)*V(:,idx); else; bn=zeros(M,1); end
            B_s{n}=bn*bn';
        end
        
        %% Trajectory (gradient projection)
        if ~strcmp(scheme_type,'straight')
            ss = 5;
            for tt = 1:10
                gr = zeros(2,N);
                for n = 2:N-1
                    for dim = 1:2
                        rp=rho; rm=rho; rp(dim,n)=rho(dim,n)+0.5; rm(dim,n)=rho(dim,n)-0.5;
                        Rp=0; Rm=0;
                        for nn=1:N
                            hup=chan(rp(:,nn),s_u); hep=chan(rp(:,nn),s_e);
                            Rp=Rp+max(0,log2(1+real(hup'*B_s{nn}*hup)/(real(hup'*As_s{nn}*hup)+sigma2))-log2(1+real(hep'*B_s{nn}*hep)/(real(hep'*As_s{nn}*hep)+sigma2)));
                            hum=chan(rm(:,nn),s_u); hem=chan(rm(:,nn),s_e);
                            Rm=Rm+max(0,log2(1+real(hum'*B_s{nn}*hum)/(real(hum'*As_s{nn}*hum)+sigma2))-log2(1+real(hem'*B_s{nn}*hem)/(real(hem'*As_s{nn}*hem)+sigma2)));
                        end
                        gr(dim,n)=(Rp-Rm);
                    end
                end
                rho_new=rho+ss*gr; rho_new(:,1)=rho_I; rho_new(:,N)=rho_F;
                for n=2:N
                    d=rho_new(:,n)-rho_new(:,n-1);
                    if norm(d)>Vmax; rho_new(:,n)=rho_new(:,n-1)+Vmax*d/norm(d); end
                end
                rho=rho_new;
            end
        end
    end

    % Final rate
    avg_sec_rate = 0;
    for n = 1:N
        hu=chan(rho(:,n),s_u); he=chan(rho(:,n),s_e);
        su=real(hu'*B_s{n}*hu)/(real(hu'*As_s{n}*hu)+sigma2);
        se=real(he'*B_s{n}*he)/(real(he'*As_s{n}*he)+sigma2);
        avg_sec_rate=avg_sec_rate+max(0,log2(1+su)-log2(1+se));
    end
    avg_sec_rate = avg_sec_rate/N;
end