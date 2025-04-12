function awgn_simulation_with_huffman_8psk()
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

    %% 2. 8-PSK Simulation Parameters
    EbN0_dB = 0:20;
    total_ber = zeros(size(EbN0_dB));
    total_ser = zeros(size(EbN0_dB));
    num_trials = 100;

    for trial = 1:num_trials
        fprintf(trial+"\n");
        ber = zeros(size(EbN0_dB));
        ser = zeros(size(EbN0_dB));

        % Ensure multiple of 3 bits (8-PSK = 3 bits/symbol)
        if mod(length(encoded_data), 3) ~= 0
            encoded_data = [encoded_data; zeros(3 - mod(length(encoded_data), 3), 1)];
        end
        numBits = length(encoded_data);
        numSymbols = numBits / 3;

        % 8-PSK Constellation (Gray Coded)
        angles = 2*pi*(0:7)/8 + pi/8;
        constellation = exp(1j*angles).';

        % Gray code mapping for 3-bit
        gray_mapping = [0 1 3 2 6 7 5 4];

        %% 3. Main Simulation Loop
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 3 * EbN0;
            noiseVar = 1/(2*EsN0);

            total_symbol_errors = 0;
            rx_encoded_data = zeros(size(encoded_data));

            % Process 3-bit symbols
            for sym_start = 1:3:numBits
                sym_end = sym_start + 2;
                tx_bits = encoded_data(sym_start:sym_end);

                % Modulation
                tx_sym = bi2de(tx_bits', 'left-msb');
                tx_sym_gray = gray_mapping(tx_sym + 1);
                txSignal = constellation(tx_sym_gray + 1);

                % AWGN Channel
                noise = sqrt(noiseVar)*(randn(1) + 1j*randn(1));
                rxSignal = txSignal + noise;

                % Demodulation
                [~, rx_sym_gray] = min(abs(rxSignal - constellation).^2);
                rx_sym_gray = rx_sym_gray - 1;
                rx_sym = find(gray_mapping == rx_sym_gray) - 1;
                rx_bits = de2bi(rx_sym, 3, 'left-msb')';
                rx_encoded_data(sym_start:sym_end) = rx_bits;

                % Symbol error check
                total_symbol_errors = total_symbol_errors + any(tx_bits ~= rx_bits);
            end

            % Huffman Decoding
            try
                rx_symbols = huffmandeco(rx_encoded_data, dict);
            catch
                rx_symbols = [];
            end

            % Convert decoded symbols to bitstream
            if ~isempty(rx_symbols)
                rx_bitstream = reshape(de2bi(rx_symbols, 8, 'left-msb')', [], 1);
            else
                rx_bitstream = [];
            end

            %% BER Calculation (Handles Length Mismatches)
            min_len = min(length(original_bitstream), length(rx_bitstream));
            overlap_errors = sum(original_bitstream(1:min_len) ~= rx_bitstream(1:min_len));

            if length(original_bitstream) > length(rx_bitstream)
                excess_errors = length(original_bitstream) - length(rx_bitstream);
            else
                excess_errors = length(rx_bitstream) - length(original_bitstream);
            end

            total_bit_errors = overlap_errors + excess_errors;
            ber(i) = total_bit_errors / length(original_bitstream);
            ser(i) = total_symbol_errors / numSymbols;
        end

        total_ber = total_ber + ber;
        total_ser = total_ser + ser;
    end

    avg_ber = total_ber / num_trials;
    avg_ser = total_ser / num_trials;

    %% 4. Plot Results
    figure;
    semilogy(EbN0_dB, avg_ber, 'bo-', 'LineWidth', 1.5, ...
             'DisplayName', 'Average End-to-End BER');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('Average End-to-End BER with Huffman + 8-PSK (100 Trials)');
    legend;

    figure;
    semilogy(EbN0_dB, avg_ser, 'rs--', 'LineWidth', 1.5, ...
             'DisplayName', 'Average Huffman Symbol Error Rate');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('Average Huffman Symbol Error Rate for 8-PSK (100 Trials)');
    legend;

    %% 5. Save Results to CSV
    results_table = table(EbN0_dB(:), avg_ber(:), avg_ser(:), ...
        'VariableNames', {'EbN0_dB', 'EndToEndBER', 'HuffmanSER'});
    writetable(results_table, 'ber_huffman_8psk_end_to_end.csv');

    fprintf('\nAverage results saved to "ber_huffman_8psk_end_to_end.csv"\n');
end
