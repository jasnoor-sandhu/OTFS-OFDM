clear
clc

%% Set Parameters

mod_method = 'QPSK';

% IFFT/FFT size
n_fft = 64;

% Size of cyclic prefix extension
n_cpe = 16;

% Target SNR (dB)
snr = 10;

% Number of channel taps (1 = no channel)
n_taps = 8;

% Channel estimation method: none,LS
ch_est_method = 'LS';

% Option to save plot to file
save_file = 0;

% Calculate modulation order from modulation method
mod_methods = {'BPSK', 'QPSK', '8PSK', '16QAM', '32QAM', '64QAM'} ;
mod_order = find(ismember(mod_methods, mod_method) );

%% Input data to binary stream
im = imread('iitbhu.bmp');
im_bin = dec2bin(im(:))';           %converts decimal numbers to binary representation
im_bin = im_bin(:);                 %concatenate all the binary representations into a single column vector

%% Binary Steam to Symbols
%Parse binary stream into mod_order bit symbols

sym_rem = mod(mod_order-mod(length(im_bin),mod_order),mod_order);          %remaining symbols
padding = repmat('0',sym_rem,1);                                           %padding is used to make the length of im_bin a multiple of mod_order
im_bin_padded = [im_bin;padding];                                          %concatinating the padding                                        
cons_data = reshape(im_bin_padded,mod_order,length(im_bin_padded)/mod_order)';      %reshaping s.t. each row have one symbol
cons_sym_id = bin2dec(cons_data);   %Binary to Decimal Conversion

%% symbol modulation
%bpsk
if mod_order == 1
    mod_ind = 2^(mod_order-1);                  %number of symbols, here 1
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;           %phase values
    in_phase = cos(n);
    quardrature = sin(n);
    symbol_book = (in_phase + quardrature*1i)';
end

%phase shift keying about unit circle
if mod_order==2 || mod_order == 3
    mod_ind = 2^(mod_order-1);
    n = 0:pi/mod_ind:2*pi-pi/mod_ind;
    in_phase = cos(n+pi/4);
    quardrature = sin(n+pi/4);
    symbol_book = (in_phase + quardrature*1i)';
end

%16QAM 64QAM modulation
if mod_order==4 || mod_order == 6
    mod_ind = sqrt(2^mod_order);
    in_phase = repmat(linspace(-1, 1, mod_ind) , mod_ind, 1) ;
    quadrature = repmat(linspace(-1,1, mod_ind)', 1, mod_ind) ;
    symbol_book = in_phase(:)+quadrature(:)*1i;
end

%32QAM modulation
%Generate 6*6 constellation and remove corners
if mod_order == 5
    mod_ind = 6;
    in_phase = repmat(linspace(-1, 1, mod_ind) , mod_ind, 1) ;
    quadrature = repmat(linspace(-1,1, mod_ind)', 1, mod_ind) ;
    symbol_book = in_phase(:) + quadrature(:)*1i;
    symbol_book = symbol_book([2:5 7:30 32:35]);
end

%modulate data according to symbol book
X = symbol_book(cons_sym_id+1);         %Mapping the data to constellation points

%% Use ifft to move to time domain

fft_rem = mod(n_fft - mod(length(X),n_fft),n_fft);              %number of remaining samples required to make the length of the signal X a multiple of n_fft.
X_padded = [X; zeros(fft_rem, 1)];                              %appends the padding
X_blocks = reshape(X_padded,n_fft,length(X_padded)/n_fft);      %reshape to perform ifft , make lenght a multiple of n_fft (serial to parallel shifter)
x = ifft(X_blocks);                                             %inverse fast fourirer transform (ifft)

%cyclic prefix insertion
x_cpe = [x(end-n_cpe+1:end,:);x];               %adds cyclic prefix of length n_cpe
x_s = x_cpe(:);                                 %transforms into column vector (Parallel to serial shifter)

%% Add AWGN

data_pwr = mean(abs(x_s.^2));                   %data power


noise_pwr = data_pwr/10^(snr/10);               %noise power
noise = normrnd(0,sqrt(noise_pwr/2),size(x_s))+ normrnd(0,sqrt(noise_pwr/2),size(x_s))*1i;      %generation of white gausssian noise
x_s_noise = x_s + noise;

%SNR of the noisy signal
snr_meas = 10* log10(mean(abs(x_s.^2))/mean(abs(noise.^2)));        %measured SNR

%% Apply fading channel
g = exp(-(0:n_taps-1));
g = g/norm(g);                                              %normalize fading coefficients
x_s_noise_fading = conv(x_s_noise,g,'same');                %convolution with fading channel

%% Use FFT to move to frequency domain

x_p = reshape(x_s_noise_fading,n_fft+n_cpe,length(x_s_noise_fading)/(n_fft+n_cpe));
x_p_cpr = x_p(n_cpe+1:end,:);                               %removing cyclic prefix n_cpe


X_hat_blocks = fft(x_p_cpr);                                %fast fourier transform

%% Estimate channel
if n_taps > 1
    switch(ch_est_method)
        case 'none'
        case 'LS'
            G = X_hat_blocks(:,1)./X_blocks(:,1);           %Least Squares channel estimation
            X_hat_blocks = X_hat_blocks./repmat(G,1,size(X_hat_blocks,2));
    end
end

%% Symbol Demodulation
%remove fft padding
X_hat = X_hat_blocks(:);
X_hat = X_hat(1:end-fft_rem);

%recover data from modulated symbols
rec_syms = knnsearch([real(symbol_book) imag(symbol_book)], [real(X_hat) imag(X_hat)])-1;

%parse to binary stream & remove symbol padding
rec_syms_cons = dec2bin(rec_syms);
rec_im_bin = reshape(rec_syms_cons', numel(rec_syms_cons),1);
rec_im_bin = rec_im_bin(1:end-sym_rem);
ber = sum(abs(rec_im_bin-im_bin)/length(im_bin));

%% Recover Image
rec_im = reshape(rec_im_bin, 8 , numel(rec_im_bin)/8);
rec_im = uint8(bin2dec(rec_im'));
rec_im = reshape(rec_im,size(im));

%% Gernerate Plots
%Transmit Constellation
subplot(2,2,1);
plot(X,'x','linewidth',2,'markersize',10);
xlim([-2 2]);
ylim([-2 2]);
xlabel('inPhase');
ylabel('Quardrature');
if n_taps>1
    title(sprintf('\\bfTransmit Constellation\n\\rm%s Modulation\n Multipath Channel Taps: %d',mod_method,n_taps));
else
    title(sprintf('\\bfTransmit Constellation\n\\rm%s Modulation',mod_method));
end
grid on;

%recovered constellation
subplot(2,2,2);
plot(X_hat(1:50:end),'x','markersize',3);
xlim([-2 2]);
ylim([-2 2]);
xlabel('inPhase');
ylabel('Quardrature');
if n_taps>1
    title(sprintf('\\bfRecieved Constellation\n\\rmmeasured SNR : %.2f dB\n Channel Estimation : %s',snr_meas,ch_est_method));
else
    title(sprintf('\\bfRecieved Constellation\n\\rmmeasured SNR : %.2f dB',snr_meas));
end
grid on;

%original Image
subplot(2,2,3);
imshow(im);
title('\bfTransmit Image');

%Recieved Image
subplot(2,2,4);
imshow(rec_im);
title(sprintf('\\bfRecovered Image\n\\rmBER = %.2g',ber));

%position Figure
set(gcf,'Position',[680 287 698 691]);

%save figure
%if save_file
%    print(sprintf('Plots/%s_%.0ffft_%.0fcpe_%.0fdB_%.0ftaps_%s',mod_method,n_fft,n_cpe,snr,n_taps,ch_est_method),'-dmeta');
%end
