function awgn_simulation_with_huffman_qpsk()
    %% 1. Input Processing and Huffman Encoding
    filename = 'input.txt';
    text_data = fileread(filename);
    ascii_values = uint8(text_data);
    
    original_bitstream = de2bi(ascii_values, 8, 'left-msb')';
    original_bitstream = original_bitstream(:);
    
    symbols = bi2de(reshape(original_bitstream, 8, [])', 'left-msb');
    [unique_syms, ~, idx] = unique(symbols);
    counts = accumarray(idx, 1);
    prob = counts / sum(counts);
    dict = huffmandict(unique_syms, prob);
    encoded_data = huffmanenco(symbols, dict);
    
    if mod(length(encoded_data), 2) ~= 0
        encoded_data = [encoded_data; zeros(2 - mod(length(encoded_data), 2), 1)];
    end
    numBits = length(encoded_data);
    numSymbols = numBits / 2;

    constellation = [1 + 1j, -1 + 1j, -1 - 1j, 1 - 1j];
    gray_mapping = [0 1 3 2];
    EbN0_dB = 0:20;
    
    total_ber = zeros(size(EbN0_dB));
    total_ser = zeros(size(EbN0_dB));
    num_repeats = 100;
    
    for repeat = 1:num_repeats
        fprintf(repeat+"\n");
        ber = zeros(size(EbN0_dB));
        ser = zeros(size(EbN0_dB));

        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 2 * EbN0;
            noiseVar = 1/(2*EsN0);

            total_symbol_errors = 0;
            rx_encoded_data = zeros(size(encoded_data));

            for sym_start = 1:2:numBits
                sym_end = sym_start + 1;
                tx_bits = encoded_data(sym_start:sym_end);

                tx_sym = bi2de(tx_bits', 'left-msb');
                tx_sym_gray = gray_mapping(tx_sym + 1);
                txSignal = constellation(tx_sym_gray + 1);

                noise = sqrt(noiseVar)*(randn(1) + 1j*randn(1));
                rxSignal = txSignal + noise;

                [~, rx_sym_gray] = min(abs(rxSignal - constellation).^2);
                rx_sym_gray = rx_sym_gray - 1;
                rx_sym = find(gray_mapping == rx_sym_gray) - 1;
                rx_bits = de2bi(rx_sym, 2, 'left-msb')';
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
            ber(i) = total_bit_errors / length(original_bitstream);
            ser(i) = total_symbol_errors / numSymbols;
        end

        total_ber = total_ber + ber;
        total_ser = total_ser + ser;
    end

    avg_ber = total_ber / num_repeats;
    avg_ser = total_ser / num_repeats;

    for i = 1:length(EbN0_dB)
        fprintf('Eb/N0 = %2d dB, Avg BER = %.2e, Avg Huffman SER = %.2e\n', EbN0_dB(i), avg_ber(i), avg_ser(i));
    end

    %% 4. Plot Results
    figure;
    semilogy(EbN0_dB, avg_ber, 'bo-', 'LineWidth', 1.5, ...
             'DisplayName', 'Average End-to-End BER');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('Average End-to-End BER with Huffman + QPSK');
    legend;

    figure;
    semilogy(EbN0_dB, avg_ser, 'rs--', 'LineWidth', 1.5, ...
             'DisplayName', 'Average Huffman Symbol Error Rate');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('Average Huffman Symbol Error Rate for QPSK');
    legend;

    %% 5. Save Results to CSV
    results_table = table(EbN0_dB(:), avg_ber(:), avg_ser(:), ...
        'VariableNames', {'EbN0_dB', 'EndToEndBER', 'HuffmanSER'});
    writetable(results_table, 'ber_huffman_qpsk_end_to_end.csv');
    
    fprintf('\nResults saved to "ber_huffman_qpsk_end_to_end.csv"\n');
end
