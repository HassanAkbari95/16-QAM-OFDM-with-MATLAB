close all;clear;clc;
%##################################
data=randi([0 1],1,9600);
ii=1;
BER_row=1; %for BER rows
pilot_data1=zeros(1,52);
pilot_data=zeros(1,64);
rx_cycle_extention_signal=zeros(1,64);
synched_sig1=zeros(1,59);
synched_signal=zeros(1,48);
BER=zeros(100,25);
ber=zeros(1,25);
errors=0;
 
for frame=1:100
framed_data=data(ii:ii+95);
ii=ii+96;
 
% Convolutionally encoding data 
constlength=7;
polynomial = [171 133];    % Polynomial
trellis = poly2trellis(constlength, polynomial);
coded_data = convenc(framed_data, trellis);
 
%Interleaving coded data
coded_data_length=length(coded_data);
each_4bit_for_QAM=coded_data_length/4;
matrix=reshape(coded_data,each_4bit_for_QAM,4);
inter_leaved_ddata = matintrlv(matrix',2,2); % Interleave.
 
% Binary to decimal conversion
dec=bi2de(inter_leaved_ddata','left-msb');
 
%16-QAM Modulation
M=16;
y = qammod(dec,M);
 
% Pilot insertion
pilot=3+3j;
 
ll=1;
for jj=(1:13:52)
    
    pilot_data1(jj)=pilot;
 
    for kk=(jj+1:jj+12)
        pilot_data1(kk)=y(ll);
        ll=ll+1;
    end
end
 
pilot_data1=pilot_data1';   % size of pilt_data =52
pilot_data(1:52)=pilot_data1(1:52);    % upsizing to 64
pilot_data(13:64)=pilot_data1(1:52);   % upsizing to 64
 
for i=1:52
    
    pilot_data(i+6)=pilot_data1(i);
    
end
 
% IFFT
ifft_signal=ifft(pilot_data',64);
 
% Adding Cyclic Extension
cycle_extension_data=zeros(80,1);
cycle_extension_data(1:16)=ifft_signal(49:64);
for i=1:64
    
    cycle_extension_data(i+16)=ifft_signal(i);
    
end
 
% Channel
 
 BER_column=1;
for SNR=0:2:48
 
ofdm_sig=awgn(cycle_extension_data,SNR,'measured'); % Adding white Gaussian Noise
 
%RX
%Removing Cyclic Extension
for i=1:64
    rx_cycle_extention_signal(i)=ofdm_sig(i+16);
end
 
% FFT
fft_signal=fft(rx_cycle_extention_signal,64);
 
% Pilot Synch%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:52
    synched_sig1(i)=fft_signal(i+6);
end
 
mm=1;
 
for i=(1:13:52)
        
    for j=(i+1:i+12)
        synched_signal(mm)=synched_sig1(j);
        mm=mm+1;
    end
end
 
% Demodulation
demodulated_data= qamdemod(synched_signal,16);
 
% Decimal to binary conversion
binary=de2bi(demodulated_data','left-msb');
binary=binary';
 
% De-Interleaving
deinter_leaved_data = matdeintrlv(binary,2,2); % De-Interleave
deinter_leaved_data=deinter_leaved_data';
deinter_leaved_data=deinter_leaved_data(:)';
 
%Decoding data
decode_data =vitdec(deinter_leaved_data,trellis,5,'trunc','hard');  % decoding datausing veterbi decoder
rx_data=decode_data;
 
% Calculating BER
rx_data=rx_data(:)';
errors=nnz(xor(framed_data,rx_data));
 
BER(BER_row,BER_column)=errors/length(framed_data);
BER_column=BER_column+1;
 
 end 
 BER_row=BER_row+1;
end 
 
% Time averaging for optimum results
 
for col=1:25        %%%change if SNR loop Changed
    ber(1,col)=0;  
    for row=1:100
  
    
        ber(1,col)=ber(1,col)+BER(row,col);
    end
end
ber=ber./100; 
 
SNR=0:2:48;
semilogy(SNR,ber,'b',SNR,ber,'r*');
title('BER/SNR');ylabel('BER');xlabel('SNR (dB)');grid on
