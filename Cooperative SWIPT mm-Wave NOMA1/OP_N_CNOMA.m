close all
clear all
%
SNR_dB = -10:1:20; % in dB
%
load data_OP_N_co_sim.dat
load data_OP_N_co_ana.dat
load data_OP_N_sim.dat
load data_OP_N_ana.dat
%% plot
semilogy(SNR_dB,data_OP_N_co_sim,'r>:',...
    SNR_dB,data_OP_N_co_ana,'r-')
hold on 
semilogy(SNR_dB,data_OP_N_sim,'o:',...
    SNR_dB,data_OP_N_ana,'*-')
legend('Location','southwest','User N NOMA (sim.)', 'User N NOMA (ana.)','User N SWIPT-CNOMA (sim.)','User F SWIPT-CNOMA (sim.)')
%
xlabel('SNR (dB)')
ylabel('OP')