function awgn_simulation_with_huffman_32psk()
    %% 0. Initialize Averaging Variables
    num_runs = 100;
    EbN0_dB = 0:20;
    ber_total = zeros(size(EbN0_dB));
    ser_total = zeros(size(EbN0_dB));

    for run = 1:num_runs
        fprintf('\n--- Run %d/%d ---\n', run, num_runs);

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
        
        %% 2. 32-PSK Simulation Parameters
        ber = zeros(size(EbN0_dB));
        ser = zeros(size(EbN0_dB));
        
        % Ensure multiple of 5 bits (32-PSK = 5 bits/symbol)
        if mod(length(encoded_data), 5) ~= 0
            encoded_data = [encoded_data; zeros(5 - mod(length(encoded_data), 5), 1)];
        end
        numBits = length(encoded_data);
        numSymbols = numBits / 5;
        
        % 32-PSK Constellation (Gray Coded)
        angles = 2*pi*(0:31)/32 + pi/32;
        constellation = exp(1j*angles).';
        
        % Gray code mapping
        gray_mapping = [0  1  3  2  6  7  5  4 12 13 15 14 10 11  9  8 ...
                       24 25 27 26 30 31 29 28 20 21 23 22 18 19 17 16];
        
        %% 3. Main Simulation Loop
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 5 * EbN0;
            noiseVar = 1/(2*EsN0);
            
            total_symbol_errors = 0;
            rx_encoded_data = zeros(size(encoded_data));
            
            % Process 5-bit symbols
            for sym_start = 1:5:numBits
                sym_end = sym_start+4;
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
                rx_bits = de2bi(rx_sym, 5, 'left-msb')';
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
        
        % Accumulate results
        ber_total = ber_total + ber;
        ser_total = ser_total + ser;
    end

    %% Average over all runs
    ber_avg = ber_total / num_runs;
    ser_avg = ser_total / num_runs;

    %% 4. Plot Results
    figure;
    semilogy(EbN0_dB, ber_avg, 'bo-', 'LineWidth', 1.5, ...
             'DisplayName', 'Average End-to-End BER');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('Average End-to-End BER with Huffman + 32-PSK');
    legend;

    figure;
    semilogy(EbN0_dB, ser_avg, 'rs--', 'LineWidth', 1.5, ...
             'DisplayName', 'Average Huffman Symbol Error Rate');
    grid on;
    xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('Average Huffman Symbol Error Rate for 32-PSK');
    legend;

    %% 5. Save Results to CSV
    results_table = table(EbN0_dB(:), ber_avg(:), ser_avg(:), ...
        'VariableNames', {'EbN0_dB', 'EndToEndBER', 'HuffmanSER'});
    writetable(results_table, 'ber_huffman_32psk_end_to_end.csv');
    
    fprintf('\nAveraged results saved to "ber_huffman_32psk_end_to_end.csv"\n');
end
