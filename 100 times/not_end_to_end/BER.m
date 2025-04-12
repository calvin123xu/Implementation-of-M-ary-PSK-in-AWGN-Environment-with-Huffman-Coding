% Multi-PSK BER Comparison Plot
files = {
    'ber_huffman_bpsk.csv',   'BPSK';
    'ber_huffman_qpsk.csv',   'QPSK';
    'ber_huffman_m8psk.csv',  '8-PSK';
    'ber_huffman_m16psk.csv', '16-PSK';
    'ber_huffman_m32psk.csv', '32-PSK';
    'ber_huffman_m64psk.csv', '64-PSK'
};

colors = lines(size(files, 1));  % Generate distinguishable colors
figure;
hold on;
grid on;

for i = 1:size(files, 1)
    filename = files{i, 1};
    label = files{i, 2};
    
    T = readtable(filename);
    disp(['File: ', filename]);
    disp(T(1:5, :)); % Display first 5 rows
    semilogy(T.EbN0_dB, T.SimulatedBER, '-o', ...
        'DisplayName', label, ...
        'LineWidth', 1.8, ...
        'Color', colors(i,:), ...
        'MarkerSize', 6);
end

xlabel('Eb/N₀ (dB)');
ylabel('Simulated Bit Error Rate (BER)');
title('BER Comparison of Huffman-Coded M-PSK Schemes over AWGN');
legend('Location', 'southwest');
axis([0 20 1e-6 1]);
set(gca, 'YScale', 'log');  % Log scale for BER
grid on;
hold off;