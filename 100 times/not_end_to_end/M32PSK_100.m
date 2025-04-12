function awgn_simulation_with_huffman_32psk()
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
    EbN0_dB = 0:20;
    ber = zeros(size(EbN0_dB));
    ser = zeros(size(EbN0_dB));
    
    % Ensure multiple of 5 bits (since 32-PSK carries 5 bits/symbol)
    if mod(length(encoded_data), 5) ~= 0
        encoded_data = [encoded_data; zeros(5-mod(length(encoded_data),5),1)];
    end
    numBits = length(encoded_data);
    numSymbols = numBits/5;
    
    % 32-PSK Constellation (Gray Coded)
    angles = 2*pi*(0:31)/32 + pi/32; % Add π/32 rotation for optimal detection
    constellation = exp(1j*angles).';
    
    % Gray code mapping for 32 symbols
    gray_mapping = [0  1  3  2  6  7  5  4 12 13 15 14 10 11  9  8 ...
                   24 25 27 26 30 31 29 28 20 21 23 22 18 19 17 16];
    bit_mapping = de2bi(gray_mapping, 5, 'left-msb');
    
    %% 3. Main Simulation Loop (100 Averaging Runs)
    runs = 100;
    for run = 1:runs
        fprintf(run+"\n");
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 5 * EbN0; % 5 bits/symbol for 32-PSK
            noiseVar = 1/(2*EsN0);
            
            total_bit_errors = 0;
            total_symbol_errors = 0;
            
            % Process in symbol chunks (5 bits at a time)
            for sym_start = 1:5:numBits
                sym_end = min(sym_start+4, numBits);
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
                rx_bits = de2bi(rx_sym, 5, 'left-msb')';
                
                % Error calculation
                symbol_errors = (rx_sym_gray ~= tx_sym_gray);
                bit_errors = sum(rx_bits ~= tx_bits);
                
                total_symbol_errors = total_symbol_errors + symbol_errors;
                total_bit_errors = total_bit_errors + bit_errors;
            end
            
            ser(i) = ser(i) + total_symbol_errors / numSymbols;
            ber(i) = ber(i) + total_bit_errors / numBits;
        end
    end
    
    % Average the results
    ber = ber / runs;
    ser = ser / runs;
    
    for i = 1:length(EbN0_dB)
        fprintf('Eb/N0 = %2d dB, BER = %.2e, SER = %.2e\n', EbN0_dB(i), ber(i), ser(i));
    end

    %% 4. Theoretical Calculations and Plotting
    % Approximate theoretical SER for 32-PSK
    theoretical_ser = 2*qfunc(sqrt(2*5*10.^(EbN0_dB/10))*sin(pi/32));
    
    % Approximate theoretical BER for Gray-coded 32-PSK
    theoretical_ber = theoretical_ser / 5; % Gray coding approximation
    
    figure;
    semilogy(EbN0_dB, ber, 'bo-', 'DisplayName', 'Simulated BER');
    hold on;
    semilogy(EbN0_dB, theoretical_ber, 'r-', 'DisplayName', 'Theoretical BER');
    grid on; xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('32-PSK Performance with Huffman Coding');
    legend('Location', 'southwest'); axis([0 20 1e-6 1]);
    
    figure;
    semilogy(EbN0_dB, ser, 'mo-', 'DisplayName', 'Simulated SER');
    hold on;
    semilogy(EbN0_dB, theoretical_ser, 'c-', 'DisplayName', 'Theoretical SER');
    grid on; xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('32-PSK Symbol Error Rate');
    legend('Location', 'southwest'); axis([0 20 1e-6 1]);
    
    %% 7. Save BER Results to CSV
    T = table(EbN0_dB(:), ber(:), theoretical_ber(:), ...
        'VariableNames', {'EbN0_dB', 'SimulatedBER', 'TheoreticalBER'});
    writetable(T, 'ber_huffman_m32psk.csv');
    fprintf('\nBER results saved to "ber_huffman_m32psk.csv"\n');
end
