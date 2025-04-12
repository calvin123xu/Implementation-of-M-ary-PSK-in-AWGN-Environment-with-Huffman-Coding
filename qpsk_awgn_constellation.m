clc; clear; close all;

%% create constellation diagram

%% 1. Generate Random Bitstream
num_bits = 1000;
bitstream = randi([0 1], 1, num_bits);

if mod(length(bitstream), 2) ~= 0
    bitstream = [bitstream, 0]; % Ensure even length
end

%% 2. QPSK Modulation (with Normalization)
bit_pairs = reshape(bitstream, 2, []).';
symbol_map = [-1+1j, -1-1j, 1+1j, 1-1j] / sqrt(2); % Normalize power
index = bi2de(bit_pairs, 'left-msb') + 1;
qpsk_symbols = symbol_map(index);

%% 3. Add AWGN (With Proper SNR)
EbN0_dB = 5; 
SNR_dB = EbN0_dB + 10*log10(2); % Convert Eb/N0 to SNR
received_symbols = awgn(qpsk_symbols, SNR_dB, 'measured');

%% 4. QPSK Demodulation (Corrected)
% Decision regions for each symbol
demodulated_bits = zeros(size(bit_pairs));
for i = 1:length(received_symbols)
    % Get real and imaginary parts
    I = real(received_symbols(i));
    Q = imag(received_symbols(i));
    
    % Decision making
    if I >= 0 && Q >= 0
        demodulated_bits(i,:) = [1 0];
    elseif I < 0 && Q >= 0
        demodulated_bits(i,:) = [0 0];
    elseif I < 0 && Q < 0
        demodulated_bits(i,:) = [0 1];
    else % I >= 0 && Q < 0
        demodulated_bits(i,:) = [1 1];
    end
end
demodulated_bits = reshape(demodulated_bits.', 1, []);

%% 5. Compute Bit Error Rate (BER)
num_errors = sum(bitstream ~= demodulated_bits);
BER = num_errors / length(bitstream);

%% 6. Display Results
disp(['Bit Error Rate (BER): ', num2str(BER)]);
figure;
scatter(real(received_symbols), imag(received_symbols), 'bo');
hold on;
scatter(real(symbol_map), imag(symbol_map), 'rx', 'LineWidth', 2);
legend('Received Symbols', 'Ideal Constellation');
xlabel('In-Phase');
ylabel('Quadrature');
title('QPSK Constellation with AWGN');
grid on;

disp('First 20 Transmitted Bits:');
disp(bitstream(1:20));

disp('First 20 Received Bits:');
disp(demodulated_bits(1:20));