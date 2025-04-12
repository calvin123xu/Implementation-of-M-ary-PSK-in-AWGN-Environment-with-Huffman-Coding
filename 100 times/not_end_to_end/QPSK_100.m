function awgn_simulation_with_huffman_qpsk()
    %% 1. Generate Bitstream with Single Huffman Dictionary
    filename = 'input.txt'; % Change to your file
    fileID = fopen(filename, 'r'); 
    text_data = fread(fileID, '*char')';
    fclose(fileID);
    
    % Convert entire text to binary first to build dictionary
    ascii_values = uint8(text_data);
    bit_matrix = de2bi(ascii_values, 8, 'left-msb');
    original_bitstream = reshape(bit_matrix', 1, []);
    
    % Create single Huffman dictionary for entire file
    symbols = reshape(original_bitstream, 8, []).';
    symbols = bi2de(symbols, 'left-msb');
    [unique_symbols, ~, idx] = unique(symbols);
    counts = accumarray(idx, 1);
    probabilities = counts / sum(counts);
    dict = huffmandict(unique_symbols, probabilities);
    
    %% 2. Process in Chunks with Fixed Dictionary
    chunk_size = 10000; % Characters per chunk
    num_chars = length(text_data);
    num_chunks = ceil(num_chars/chunk_size);
    encoded_data = [];
    
    for chunk = 1:num_chunks
        start_idx = (chunk-1)*chunk_size + 1;
        end_idx = min(chunk*chunk_size, num_chars);
        chunk_data = text_data(start_idx:end_idx);
        
        % Convert chunk to binary
        ascii_values = uint8(chunk_data);
        bit_matrix = de2bi(ascii_values, 8, 'left-msb');
        chunk_bitstream = reshape(bit_matrix', 1, []);
        
        % Huffman encoding with fixed dictionary
        symbols = reshape(chunk_bitstream, 8, []).';
        symbols = bi2de(symbols, 'left-msb');
        encoded_chunk = huffmanenco(symbols, dict);
        encoded_data = [encoded_data; encoded_chunk(:)]; % Vertical concatenation
    end
    
    %% 3. QPSK Modulation and AWGN Simulation
    EbN0_dB = 0:20;
    numBits = length(encoded_data);
    
    % Ensure even number of bits for QPSK (2 bits per symbol)
    if mod(numBits, 2) ~= 0
        encoded_data = [encoded_data; 0]; % Pad with zero if odd
        numBits = numBits + 1;
    end
    
    numSymbols = numBits/2;
    ber = zeros(size(EbN0_dB));
    process_chunk_size = 1e6; % Process 1 million bits at a time (500k symbols)
    
    % Number of repetitions for averaging
    num_repeats = 100;
    
    % Store total BER for averaging
    total_ber = zeros(size(EbN0_dB));
    
    for repeat = 1:num_repeats
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            
            % For QPSK, Es/N0 = 2*Eb/N0 (2 bits per symbol)
            EsN0 = 2 * EbN0;
            noiseVar = 1 / (2 * EsN0); % Noise variance per dimension
            
            total_errors = 0;
            processed_bits = 0;
            
            % Process in chunks
            for chunk_start = 1:process_chunk_size:numBits
                chunk_end = min(chunk_start + process_chunk_size - 1, numBits);
                current_bits = encoded_data(chunk_start:chunk_end);
                
                % Ensure even number of bits in this chunk
                if mod(length(current_bits), 2) ~= 0
                    current_bits = [current_bits; 0];
                end
                
                % QPSK Modulation
                even_bits = current_bits(1:2:end);
                odd_bits = current_bits(2:2:end);
                
                % Map to QPSK symbols (Gray coding)
                % 00 -> exp(j*pi/4), 01 -> exp(j*3pi/4)
                % 11 -> exp(j*5pi/4), 10 -> exp(j*7pi/4)
                txSymbols = (1/sqrt(2)) * ((1 - 2*even_bits) + 1j*(1 - 2*odd_bits));
                
                % AWGN Channel
                noise = sqrt(noiseVar) * (randn(size(txSymbols)) + 1j*randn(size(txSymbols)));
                rxSymbols = txSymbols + noise;
                
                % QPSK Demodulation
                rx_even_bits = real(rxSymbols) < 0;
                rx_odd_bits = imag(rxSymbols) < 0;
                
                % Combine bits
                rxBits = zeros(2*length(rx_even_bits), 1);
                rxBits(1:2:end) = rx_even_bits;
                rxBits(2:2:end) = rx_odd_bits;
                
                % Calculate errors
                chunk_errors = sum(rxBits ~= current_bits(1:length(rxBits)));
                total_errors = total_errors + chunk_errors;
                processed_bits = processed_bits + length(current_bits);
            end
            
            total_ber(i) = total_ber(i) + total_errors / processed_bits;
        end
    end
    
    % Average the BER after all repetitions
    ber = total_ber / num_repeats;
    
    %% 4. Theoretical BER Calculation for QPSK
    % Theoretical BER for QPSK in AWGN channel (same as BPSK when using Gray coding)
    theoretical_ber = 0.5*erfc(sqrt(10.^(EbN0_dB/10)));
    
    % Display theoretical values
    fprintf('\nTheoretical BER values for QPSK:\n');
    for i = 1:length(EbN0_dB)
        fprintf('Eb/N0 = %2d dB, Theoretical BER = %.2e\n', EbN0_dB(i), theoretical_ber(i));
    end
    
    %% 5. Plot Results
    figure;
    semilogy(EbN0_dB, ber, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'Simulated BER (QPSK)');
    hold on;
    semilogy(EbN0_dB, theoretical_ber, 'r-', 'LineWidth', 2, 'DisplayName', 'Theoretical BER (QPSK)');
    grid on;
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('AWGN Channel Performance with Huffman Coding (QPSK)');
    legend('Location', 'best');
    axis([0 20 1e-6 1]);
    hold off;
    
    %% 6. Additional Analysis
    % Calculate coding gain at BER = 1e-4
    target_ber = 1e-4;
    [~, sim_idx] = min(abs(ber - target_ber));
    [~, theo_idx] = min(abs(theoretical_ber - target_ber));
    
    if ~isempty(sim_idx) && ~isempty(theo_idx)
        coding_gain = EbN0_dB(theo_idx) - EbN0_dB(sim_idx);
        fprintf('\nCoding gain at BER ≈ %0.0e: %.2f dB\n', target_ber, coding_gain);
    end
    %% 7. Save BER Results to CSV
    T = table(EbN0_dB(:), ber(:), theoretical_ber(:), ...
        'VariableNames', {'EbN0_dB', 'SimulatedBER', 'TheoreticalBER'});
    writetable(T, 'ber_huffman_qpsk.csv');
    fprintf('\nBER results saved to "ber_huffman_qpsk.csv"\n');
end
