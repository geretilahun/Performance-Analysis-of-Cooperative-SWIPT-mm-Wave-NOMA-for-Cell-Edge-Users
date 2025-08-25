close all
clear all
%
SNR_dB = -10:1:20; % in dB
%
load data_OP_F_co_sim.dat
load data_OP_F_co_ana.dat
load data_OP_F_sim.dat
load data_OP_F_ana.dat
%
load data_OP_N_co_sim.dat
load data_OP_N_co_ana.dat
load data_OP_N_sim.dat
load data_OP_N_ana.dat
%% plot
semilogy(SNR_dB,data_OP_F_co_sim,'r>:',...
    SNR_dB,data_OP_F_co_ana,'r-')
hold on 
semilogy(SNR_dB,data_OP_F_sim,'o:',...
    SNR_dB,data_OP_F_ana,'*-')
hold on
semilogy(SNR_dB,data_OP_N_co_sim,'b<:',...
    SNR_dB,data_OP_N_co_ana,'b-')
hold on 
semilogy(SNR_dB,data_OP_N_sim,'gs:',...
    SNR_dB,data_OP_N_ana,'g-')
legend('Location','southwest','User F NOMA (sim.)','User F NOMA (ana.)','User F SWIPT-CNOMA (sim.)','User F SWIPT-CNOMA (ana.)','User N NOMA (sim.)', 'User N NOMA (ana.)','User N SWIPT-CNOMA (sim.)','User N SWIPT-CNOMA (ana.)')
%
xlabel('SNR (dB)')
ylabel('OP')