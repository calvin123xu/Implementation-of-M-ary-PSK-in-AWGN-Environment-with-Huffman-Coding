function awgn_simulation_with_huffman_bpsk_dual_plots()
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
    
    %% 2. BPSK Simulation Parameters
    EbN0_dB = 0:20;
    regular_ber_all = zeros(100, length(EbN0_dB));  % End-to-end bit error rate (100 trials)
    huffman_ser_all = zeros(100, length(EbN0_dB));  % Huffman symbol error rate (100 trials)
    
    % Store original Huffman codeword lengths
    codeword_lengths = cellfun(@length, dict(:,2));
    
    %% 3. Monte Carlo Simulation (100 times)
    for trial = 1:100
        fprintf('--- Trial %d ---\n', trial);
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            noiseVar = 1/(2*EbN0);
            
            % Initialize counters
            huffman_symbol_errors = 0;
            total_symbols_transmitted = 0;
            encoded_ptr = 1;
            symbol_idx = 1;
            all_rx_bits = [];

            while encoded_ptr <= length(encoded_data)
                if symbol_idx > length(symbols)
                    break;
                end

                current_symbol = symbols(symbol_idx);
                current_length = codeword_lengths(unique_syms == current_symbol);

                if encoded_ptr + current_length - 1 > length(encoded_data)
                    break; % Skip incomplete symbols
                end

                tx_bits = encoded_data(encoded_ptr : encoded_ptr + current_length - 1);

                % BPSK Modulation + AWGN
                txSignal = 2*tx_bits - 1;
                noise = sqrt(noiseVar)*randn(size(txSignal));
                rxSignal = txSignal + noise;

                % Demodulation
                rx_bits = rxSignal > 0;

                % Accumulate for end-to-end decoding
                all_rx_bits = [all_rx_bits; rx_bits];

                % Symbol error detection
                huffman_symbol_errors = huffman_symbol_errors + any(tx_bits ~= rx_bits);
                total_symbols_transmitted = total_symbols_transmitted + 1;

                % Advance pointers
                encoded_ptr = encoded_ptr + current_length;
                symbol_idx = symbol_idx + 1;
            end

            % Huffman decoding
            try
                decoded_symbols = huffmandeco(all_rx_bits, dict);
            catch
                decoded_symbols = [];
            end

            % Reconstruct bitstream
            if ~isempty(decoded_symbols)
                decoded_bits = de2bi(decoded_symbols, 8, 'left-msb')';
                decoded_bits = decoded_bits(:);
                min_len = min(length(decoded_bits), length(original_bitstream));
                bit_errors = sum(decoded_bits(1:min_len) ~= original_bitstream(1:min_len));
                regular_ber_all(trial, i) = bit_errors / min_len;
            else
                regular_ber_all(trial, i) = 1; % Total decoding failure
            end

            % Store symbol error rate
            huffman_ser_all(trial, i) = huffman_symbol_errors / total_symbols_transmitted;
        end
    end

    % Average results over 100 trials
    regular_ber = mean(regular_ber_all, 1);
    huffman_ser = mean(huffman_ser_all, 1);

    %% 4. Plotting
    figure;
    semilogy(EbN0_dB, regular_ber, 'bo-', 'LineWidth', 1.5, ...
             'MarkerFaceColor', 'b', 'DisplayName', 'End-to-End BER');
    hold on;
    semilogy(EbN0_dB, huffman_ser, 'rs--', 'LineWidth', 1.5, ...
             'MarkerFaceColor', 'r', 'DisplayName', 'Huffman Symbol Error Rate');

    % Theoretical BPSK BER
    theoretical_ber = 0.5 * erfc(sqrt(10.^(EbN0_dB/10)));
    semilogy(EbN0_dB, theoretical_ber, 'k-', 'LineWidth', 2, ...
             'DisplayName', 'Theoretical BPSK BER');

    grid on;
    xlabel('Eb/N0 (dB)');
    ylabel('Error Rate');
    title('End-to-End BER and Huffman SER with BPSK');
    legend('Location', 'southwest');
    axis([0 20 1e-6 1]);
    set(gca, 'FontSize', 12);

    %% 5. Save Results
    T = table(EbN0_dB(:), regular_ber(:), huffman_ser(:), theoretical_ber(:), ...
        'VariableNames', {'EbN0_dB', 'EndToEndBER', 'HuffmanSER', 'TheoreticalBER'});
    writetable(T, 'ber_huffman_bpsk_end_to_end.csv');
    fprintf('\nResults saved to "ber_huffman_bpsk_end_to_end.csv"\n');
end
