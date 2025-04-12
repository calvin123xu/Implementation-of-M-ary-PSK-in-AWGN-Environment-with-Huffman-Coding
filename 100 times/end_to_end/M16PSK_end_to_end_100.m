function awgn_simulation_with_huffman_16psk()
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

    %% 2. 16-PSK Simulation Parameters
    EbN0_dB = 0:20;
    ber = zeros(size(EbN0_dB));
    ser = zeros(size(EbN0_dB));

    % Repeat simulation 100 times and average
    numIterations = 100;
    for iter = 1:numIterations
        fprintf(iter+"\n");
        ber_iter = zeros(size(EbN0_dB));
        ser_iter = zeros(size(EbN0_dB));

        % Ensure multiple of 4 bits (16-PSK = 4 bits/symbol)
        if mod(length(encoded_data), 4) ~= 0
            encoded_data = [encoded_data; zeros(4 - mod(length(encoded_data), 4), 1)];
        end
        numBits = length(encoded_data);
        numSymbols = numBits / 4;

        % 16-PSK Constellation (Gray Coded)
        angles = 2*pi*(0:15)/16 + pi/16;
        constellation = exp(1j*angles).';

        % Gray code mapping
        gray_mapping = [0 1 3 2 6 7 5 4 12 13 15 14 10 11 9 8];

        %% 3. Main Simulation Loop
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 4 * EbN0;
            noiseVar = 1/(2*EsN0);

            total_symbol_errors = 0;
            rx_encoded_data = zeros(size(encoded_data));

            % Process 4-bit symbols
            for sym_start = 1:4:numBits
                sym_end = sym_start+3;
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
                rx_bits = de2bi(rx_sym, 4, 'left-msb')';
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
            ber_iter(i) = total_bit_errors / length(original_bitstream);
            ser_iter(i) = total_symbol_errors / numSymbols;
        end

        ber = ber + ber_iter;
        ser = ser + ser_iter;
    end

    % Average the results over all iterations
    ber = ber / numIterations;
    ser = ser / numIterations;

    %% 4. Plot Results
    figure;
    semilogy(EbN0_dB, ber, 'bo-', 'LineWidth', 1.5, ...
             'DisplayName', 'End-to-End BER');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('End-to-End BER with Huffman + 16-PSK');
    legend;

    figure;
    semilogy(EbN0_dB, ser, 'rs--', 'LineWidth', 1.5, ...
             'DisplayName', 'Huffman Symbol Error Rate');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('Huffman Symbol Error Rate for 16-PSK');
    legend;

    %% 5. Save Results to CSV
    results_table = table(EbN0_dB(:), ber(:), ser(:), ...
        'VariableNames', {'EbN0_dB', 'EndToEndBER', 'HuffmanSER'});
    writetable(results_table, 'ber_huffman_16psk_end_to_end.csv');

    fprintf('\nResults saved to "ber_huffman_16psk_end_to_end.csv"\n');
end
