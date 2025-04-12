function awgn_simulation_with_huffman_8psk_fixed()
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
    ber = zeros(100, length(EbN0_dB)); % Store BER for 100 iterations
    ser = zeros(100, length(EbN0_dB)); % Store SER for 100 iterations
    
    % Ensure multiple of 3 bits
    if mod(length(encoded_data), 3) ~= 0
        encoded_data = [encoded_data; zeros(3-mod(length(encoded_data),3),1)];
    end
    numBits = length(encoded_data);
    numSymbols = numBits / 3;
    
    % 8-PSK Constellation (Gray Coded)
    constellation = exp(1j*(2*pi*(0:7)/8 + pi/8)).';
    gray_mapping = [0 1 3 2 6 7 5 4]; % Proper Gray code ordering
    bit_mapping = de2bi(gray_mapping, 3, 'left-msb');
    
    %% 3. Main Simulation Loop (100 iterations)
    for j = 1:100
        fprintf(j+"\n");
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 3 * EbN0; % 3 bits/symbol
            noiseVar = 1/(2*EsN0);
            
            total_bit_errors = 0;
            total_symbol_errors = 0;
            
            % Process in symbol chunks
            for sym_start = 1:3:numBits
                sym_end = min(sym_start+2, numBits);
                tx_bits = encoded_data(sym_start:sym_end);
                
                % Symbol mapping (Gray coded)
                tx_sym = bi2de(tx_bits', 'left-msb');
                tx_sym_gray = gray_mapping(tx_sym + 1);
                txSignal = constellation(tx_sym_gray + 1);
                
                % AWGN Channel
                noise = sqrt(noiseVar)*(randn(1) + 1j*randn(1));
                rxSignal = txSignal + noise;
                
                % Demodulation (minimum distance)
                [~, rx_sym_gray] = min(abs(rxSignal - constellation).^2);
                rx_sym_gray = rx_sym_gray - 1;
                
                % Gray to binary
                rx_sym = find(gray_mapping == rx_sym_gray) - 1;
                rx_bits = de2bi(rx_sym, 3, 'left-msb')';
                
                % Error calculation
                symbol_errors = (rx_sym_gray ~= tx_sym_gray);
                bit_errors = sum(rx_bits ~= tx_bits);
                
                total_symbol_errors = total_symbol_errors + symbol_errors;
                total_bit_errors = total_bit_errors + bit_errors;
            end
            
            ser(j,i) = total_symbol_errors / numSymbols;
            ber(j,i) = total_bit_errors / numBits;
        end
    end
    
    %% 4. Averaging Results
    avg_ber = mean(ber);
    avg_ser = mean(ser);
    
    %% 5. Theoretical Calculations and Plotting
    theoretical_ser = 2*qfunc(sqrt(2*3*10.^(EbN0_dB/10))*sin(pi/8));
    theoretical_ber = theoretical_ser / 3; % Gray coding approximation
    
    figure;
    semilogy(EbN0_dB, avg_ber, 'bo-', 'DisplayName', 'Average Simulated BER');
    hold on;
    semilogy(EbN0_dB, theoretical_ber, 'r-', 'DisplayName', 'Theoretical BER');
    grid on; xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('8-PSK Performance with Huffman Coding');
    legend; axis([0 20 1e-6 1]);
    
    %% 6. Save BER Results to CSV
    T = table(EbN0_dB(:), avg_ber(:), theoretical_ber(:), ...
        'VariableNames', {'EbN0_dB', 'SimulatedBER', 'TheoreticalBER'});
    writetable(T, 'ber_huffman_m8psk_avg.csv');
    fprintf('\nAverage BER results saved to "ber_huffman_m8psk_avg.csv"\n');
end
