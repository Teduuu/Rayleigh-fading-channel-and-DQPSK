clc;
M = csvread('dbpsk_data.csv');
snr = M(:,1);
sim_ber = M(:,2);

N = csvread('dbpsk_data2.csv');
snr2 = N(:,1);
the_ber = N(:,2);

figure(1);
semilogy(snr,sim_ber,'-o');


hold on;

semilogy(snr2,the_ber);
grid on;
grid minor

title('DBPSK BER Wm = 2\pi*200');
xlabel('SNR_dB');
ylabel('BER');
legend('Simulated BER','Theoretical BER')
hold off;