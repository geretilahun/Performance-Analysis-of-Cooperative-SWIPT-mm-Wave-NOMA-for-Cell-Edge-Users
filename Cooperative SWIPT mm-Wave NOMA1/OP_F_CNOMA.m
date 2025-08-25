close all
clear all
%
SNR_dB = -10:1:20; % in dB
%
load data_OP_F_co_sim.dat
load data_OP_F_co_ana.dat
load data_OP_F_sim.dat
load data_OP_F_ana.dat

%% plot
semilogy(SNR_dB,data_OP_F_co_sim,'o:',...
    SNR_dB,data_OP_F_co_ana,'*-')
hold on 
semilogy(SNR_dB,data_OP_F_sim,'o:',...
    SNR_dB,data_OP_F_ana,'*-')
legend('Location','southwest','User F NOMA (sim.)', 'User F NOMA (ana.)','User F SWIPT-CNOMA (sim.)','User F SWIPT-CNOMA (sim.)')
%
xlabel('SNR (dB)')
ylabel('OP')