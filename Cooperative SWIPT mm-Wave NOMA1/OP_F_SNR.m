clear all
close all
%%
Sim_times   = 1e6;              % Monte-Carlo simulation iterations
%
G          = 3;                % # of BS antenna
rho         = .3;               % power splitting ratio
alpha       = .3;               % time fraction for EH
%
SNR_dB      = -10:1:20;         % in dB
SNR         = 10.^(SNR_dB./10);
%naN     = (10^(-7))*1e6;        % naN = -100 dBm, BW = 1 MHz
%ncN     = (10^(-6))*1e6;        % naN = -90 dBm,  BW = 1 MHz
%naF     = (10^(-7))*1e6;
%ncF     = (10^(-6))*1e6;
naN     = 1;
naF     = 1;
%
RtF         = 1;                % for User F
RtN         = 1;                % for User N
R_on_off    = 1;                % for ON/OFF
rho2        = 2^RtF - 1;        % full time slot user F
rho1        = 2^RtN - 1;        % full time slot user N
rho2h       = 2^(2*RtF) - 1;    % half time slot user F
rho1h       = 2^(2*RtN) - 1;    % half time slot user N
rho0        = 2.^R_on_off - 1;  % SNR threshold for ON/OFF
% NOMA parameters
eta         = 0.7;              % energy conversion coefficient
pF          = .8;
pN          = 1 - pF;
theta       = pF/pN;
% Path-loss model
L           = 1e3;              % 30 dB path-loss at reference distance
epsilon     = 2;                % path-loss exponent
d0          = 1;                % reference distance
dSF         = 10;
dSN         = 3;
dNF         = dSF - dSN;
%
lSN         = L*((dSN/d0)^(-epsilon)); % lambda_SN
lNF         = L*((dNF/d0)^(-epsilon)); % lambda_NF
lSF         = L*((dSF/d0)^(-epsilon)); % lambda_SF
lrsi        = .1; % lambda_rsi
% Rician parameters
K_dB        = 3;                  % Rician K-factor in dB
K           = 10.^(K_dB./10);
% Rician Distribution parameters
nu          = sqrt(K/(K+1)*lrsi); % Noncentrality parameter
sigma       = sqrt(lrsi/2/(K+1)); % Scale parameter
%% Simulation
%
for ss = 1:length(SNR_dB)
    fprintf('SNR = %d dB \n',SNR_dB(ss))
    for rr = 1:length(rho)
        % Channel modeling
    for ii = 1:G
         hiF(:,ii)  = random('Rayleigh',sqrt(lSF/2),[1,Sim_times]);
         hiN(:,ii)  = random('Rayleigh',sqrt(lSN/2),[1,Sim_times]);
     
    end
         hNF     = random('Rayleigh',sqrt(lNF/2),[1,Sim_times]);
         hrsi    = random('rician',nu,sigma,[1,Sim_times]); % Rician fading
            
            % Channel gains
            giF     = abs(hiF).^2;
            giN     = abs(hiN).^2;
            gNF     = abs(hNF).^2;
            grsi    = abs(hrsi).^2;
            % SNRs
        snr_SN_xF = (1-rho(rr)).*pF.*SNR(ss).*giN./...
            ((1-rho(rr)).*pN.*SNR(ss).*giN ...
            + naN);
        snr_SN_xN = (1-rho(rr)).*pN.*SNR(ss).*giN/...
            (naN);
        snr_SF   = pF.*SNR(ss).*giF./(pN.*SNR(ss).*giF + naF);
        snr_NF   = eta.*SNR(ss).*giN.*giF.*rho(rr)/(naF);
        %
        OP_F_sim_temp = 0;
        % Find the best antenna for User F based on end-to-end SNR
        snrFe2e(:,:) = min(snr_SN_xF(:,:),max(snr_NF(:,:),snr_SF(:,:)));
        [snrFe2e_b,I] = max(snrFe2e,[],3);
        % count outage events
         for zz = 1:Sim_times
        % for F SWIPT CNOMA
        if (snrFe2e_b(zz) < rho2h)
            OP_F_sim_temp = OP_F_sim_temp + 1;
        end
         end
         OP_F_sim(ss) = OP_F_sim_temp/Sim_times;
        %count = snrFe2e < rho2h;
       % OP_F_sim(ss) = sum(count)/Sim_times;
        %% Analysis
        a1 = (1-rho(rr))*pF*SNR(ss)/(naN);
        a2 = (1-rho(rr))*pN*SNR(ss)/(naN);
        b1 = pF * SNR(ss) / (naF);
        b2 = pN * SNR(ss) / (naF);
        c  = eta*rho(rr)*SNR(ss)/(naF);
        mu_a = rho2h/(a1-a2*rho2h);
        mu_b = rho2h/(b1-b2*rho2h);
        theta = pF/pN;
        %
        if rho2h < theta
            OP_F_ana(ss) = 1 - exp(-mu_a/lSN - mu_b/lSF)...
                - (1 - exp(-mu_b/lSF)) ...
                * (exp(-mu_a/lSN) - rho2h/lSN/lNF/c*igamma(0,mu_a/lSN));
        else
            OP_F_ana(ss) = 1;
        end
     end
end
%% plot
semilogy(SNR_dB,OP_F_sim,'o:',...
    SNR_dB,OP_F_ana,'*-')
legend('Location','southwest','User F SWIPT-CNOMA (sim.)', 'User F CNOMA (ana.)')
%
xlabel('SNR (dB)')
ylabel('OP')
%
save data_OP_F_sim.dat OP_F_sim -ascii
save data_OP_F_ana.dat OP_F_ana -ascii

