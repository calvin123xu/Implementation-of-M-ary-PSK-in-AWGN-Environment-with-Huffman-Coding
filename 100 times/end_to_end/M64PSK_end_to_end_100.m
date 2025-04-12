function awgn_simulation_with_huffman_64psk()
    %% 1. Input Processing and Huffman Encoding
    filename = 'input.txt';
    text_data = fileread(filename);
    ascii_values = uint8(text_data);

    % Convert to binary stream (original bits)
    original_bitstream = de2bi(ascii_values, 8, 'left-msb')';
    original_bitstream = original_bitstream(:);

    % Huffman encoding
    symbols = bi2de(reshape(original_bitstream, 8, [])', 'left-msb');
    [unique_syms, ~, idx] = unique(symbols);
    counts = accumarray(idx, 1);
    prob = counts / sum(counts);
    dict = huffmandict(unique_syms, prob);
    encoded_data = huffmanenco(symbols, dict);

    %% 2. 64-PSK Simulation Parameters
    EbN0_dB = 0:20;
    ber = zeros(size(EbN0_dB));
    ser = zeros(size(EbN0_dB));

    %% 100 Times Averaging Loop
    for trial = 1:100
        tmp_ber = zeros(size(EbN0_dB));
        tmp_ser = zeros(size(EbN0_dB));

        % Ensure multiple of 6 bits
        if mod(length(encoded_data), 6) ~= 0
            encoded_data = [encoded_data; zeros(6 - mod(length(encoded_data), 6), 1)];
        end
        numBits = length(encoded_data);
        numSymbols = numBits / 6;

        angles = 2*pi*(0:63)/64 + pi/64;
        constellation = exp(1j*angles).';

        gray_mapping = [0  1  3  2  6  7  5  4 12 13 15 14 10 11  9  8 ...
                       24 25 27 26 30 31 29 28 20 21 23 22 18 19 17 16 ...
                       48 49 51 50 54 55 53 52 60 61 63 62 58 59 57 56 ...
                       40 41 43 42 46 47 45 44 36 37 39 38 34 35 33 32];

        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 6 * EbN0;
            noiseVar = 1/(2*EsN0);

            total_bit_errors = 0;
            total_symbol_errors = 0;
            rx_encoded_data = zeros(size(encoded_data));

            for sym_start = 1:6:numBits
                sym_end = sym_start+5;
                tx_bits = encoded_data(sym_start:sym_end);

                tx_sym = bi2de(tx_bits', 'left-msb');
                tx_sym_gray = gray_mapping(tx_sym + 1);
                txSignal = constellation(tx_sym_gray + 1);

                noise = sqrt(noiseVar)*(randn(1) + 1j*randn(1));
                rxSignal = txSignal + noise;

                [~, rx_sym_gray] = min(abs(rxSignal - constellation).^2);
                rx_sym_gray = rx_sym_gray - 1;
                rx_sym = find(gray_mapping == rx_sym_gray) - 1;
                rx_bits = de2bi(rx_sym, 6, 'left-msb')';
                rx_encoded_data(sym_start:sym_end) = rx_bits;

                total_symbol_errors = total_symbol_errors + any(tx_bits ~= rx_bits);
            end

            try
                rx_symbols = huffmandeco(rx_encoded_data, dict);
            catch
                rx_symbols = [];
            end

            if ~isempty(rx_symbols)
                rx_bitstream = reshape(de2bi(rx_symbols, 8, 'left-msb')', [], 1);
            else
                rx_bitstream = [];
            end

            min_len = min(length(original_bitstream), length(rx_bitstream));
            overlap_errors = sum(original_bitstream(1:min_len) ~= rx_bitstream(1:min_len));

            if length(original_bitstream) > length(rx_bitstream)
                excess_errors = length(original_bitstream) - length(rx_bitstream);
            else
                excess_errors = length(rx_bitstream) - length(original_bitstream);
            end

            total_bit_errors = overlap_errors + excess_errors;
            tmp_ber(i) = total_bit_errors / length(original_bitstream);
            tmp_ser(i) = total_symbol_errors / numSymbols;
        end

        ber = ber + tmp_ber;
        ser = ser + tmp_ser;
    end

    %% Averaging Results
    ber = ber / 100;
    ser = ser / 100;

    %% 4. Plot Results
    figure;
    semilogy(EbN0_dB, ber, 'bo-', 'LineWidth', 1.5, ...
             'DisplayName', 'End-to-End BER');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('End-to-End BER with Huffman + 64-PSK');
    legend;

    figure;
    semilogy(EbN0_dB, ser, 'rs--', 'LineWidth', 1.5, ...
             'DisplayName', 'Symbol Error Rate');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('Symbol Error Rate for 64-PSK');
    legend;

    %% 5. Save Results to CSV
    results_table = table(EbN0_dB(:), ber(:), ser(:), ...
        'VariableNames', {'EbN0_dB', 'EndToEndBER', 'HuffmanSER'});
    writetable(results_table, 'ber_huffman_64psk_end_to_end.csv');

    fprintf('\nResults saved to "ber_huffman_64psk_end_to_end.csv"\n');
end
