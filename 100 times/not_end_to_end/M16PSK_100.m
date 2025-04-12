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
    ber_all = zeros(100, length(EbN0_dB)); % Store 100 runs
    ser_all = zeros(100, length(EbN0_dB));
    
    % Ensure multiple of 4 bits
    if mod(length(encoded_data), 4) ~= 0
        encoded_data = [encoded_data; zeros(4-mod(length(encoded_data),4),1)];
    end
    numBits = length(encoded_data);
    numSymbols = numBits / 4;
    
    % 16-PSK Constellation (Gray Coded)
    constellation = exp(1j*(2*pi*(0:15)/16 + pi/16)).';
    gray_mapping = [0 1 3 2 6 7 5 4 12 13 15 14 10 11 9 8]; 
    bit_mapping = de2bi(gray_mapping, 4, 'left-msb');
    
    %% 3. Main Simulation Loop over 100 runs
    for run = 1:100
        fprintf(run+"\n");
        ber = zeros(size(EbN0_dB));
        ser = zeros(size(EbN0_dB));
        
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            EsN0 = 4 * EbN0;
            noiseVar = 1/(2*EsN0);
            
            total_bit_errors = 0;
            total_symbol_errors = 0;
            
            for sym_start = 1:4:numBits
                sym_end = min(sym_start+3, numBits);
                tx_bits = encoded_data(sym_start:sym_end);
                
                tx_sym = bi2de(tx_bits', 'left-msb');
                tx_sym_gray = gray_mapping(tx_sym + 1);
                txSignal = constellation(tx_sym_gray + 1);
                
                noise = sqrt(noiseVar)*(randn(1) + 1j*randn(1));
                rxSignal = txSignal + noise;
                
                [~, rx_sym_gray] = min(abs(rxSignal - constellation).^2);
                rx_sym_gray = rx_sym_gray - 1;
                
                rx_sym = find(gray_mapping == rx_sym_gray) - 1;
                rx_bits = de2bi(rx_sym, 4, 'left-msb')';
                
                symbol_errors = (rx_sym_gray ~= tx_sym_gray);
                bit_errors = sum(rx_bits ~= tx_bits);
                
                total_symbol_errors = total_symbol_errors + symbol_errors;
                total_bit_errors = total_bit_errors + bit_errors;
            end
            
            ser(i) = total_symbol_errors / numSymbols;
            ber(i) = total_bit_errors / numBits;
        end
        
        ber_all(run, :) = ber;
        ser_all(run, :) = ser;
    end
    
    %% 4. Average Results over 100 runs
    ber = mean(ber_all, 1);
    ser = mean(ser_all, 1);
    
    %% 5. Theoretical Calculations and Plotting
    theoretical_ser = 2*qfunc(sqrt(2*4*10.^(EbN0_dB/10))*sin(pi/16));
    theoretical_ber = theoretical_ser / 4;
    
    figure;
    semilogy(EbN0_dB, ber, 'bo-', 'DisplayName', 'Simulated BER');
    hold on;
    semilogy(EbN0_dB, theoretical_ber, 'r-', 'DisplayName', 'Theoretical BER');
    grid on; xlabel('Eb/N0 (dB)'); ylabel('BER');
    title('16-PSK Performance with Huffman Coding');
    legend; axis([0 20 1e-6 1]);
    
    figure;
    semilogy(EbN0_dB, ser, 'mo-', 'DisplayName', 'Simulated SER');
    hold on;
    semilogy(EbN0_dB, theoretical_ser, 'c-', 'DisplayName', 'Theoretical SER');
    grid on; xlabel('Eb/N0 (dB)'); ylabel('SER');
    title('16-PSK Symbol Error Rate');
    legend; axis([0 20 1e-6 1]);
    
    %% 6. Save BER Results to CSV
    T = table(EbN0_dB(:), ber(:), theoretical_ber(:), ...
        'VariableNames', {'EbN0_dB', 'SimulatedBER', 'TheoreticalBER'});
    writetable(T, 'ber_huffman_m16psk.csv');
    fprintf('\nBER results saved to "ber_huffman_m16psk.csv"\n');
end
