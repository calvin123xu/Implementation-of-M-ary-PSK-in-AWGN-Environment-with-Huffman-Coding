function awgn_simulation_with_huffman_fixed()
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
    
    %% 3. AWGN Simulation with Memory Management
    EbN0_dB = 0:20;
    numBits = length(encoded_data);
    ber = zeros(size(EbN0_dB));
    process_chunk_size = 1e6; % Bits per simulation chunk
    
    % Initialize arrays to store results from each run
    ber_all = zeros(100, length(EbN0_dB));  % 100 simulations
    
    % Run the simulation 100 times
    for run = 1:100
        for i = 1:length(EbN0_dB)
            EbN0 = 10^(EbN0_dB(i)/10);
            noiseVar = 1 / (2 * EbN0);
            total_errors = 0;
            
            for chunk_start = 1:process_chunk_size:numBits
                chunk_end = min(chunk_start + process_chunk_size - 1, numBits);
                current_chunk = chunk_start:chunk_end;
                
                txSymbols = 2 * encoded_data(current_chunk) - 1;
                noise = sqrt(noiseVar) * randn(size(txSymbols));
                rxSymbols = txSymbols + noise;
                rxBits = rxSymbols > 0;
                
                total_errors = total_errors + sum(rxBits ~= encoded_data(current_chunk));
            end
            
            ber(run, i) = total_errors / numBits;  % Store BER for this run
        end
    end
    
    %% 4. Average BER Calculation
    avg_ber = mean(ber, 1);  % Average BER across all 100 runs
    
    %% 5. Theoretical BER Calculation
    % Theoretical BER for BPSK in AWGN channel
    theoretical_ber = 0.5*erfc(sqrt(10.^(EbN0_dB/10)));
    
    % Display theoretical values
    fprintf('\nTheoretical BER values:\n');
    for i = 1:length(EbN0_dB)
        fprintf('Eb/N0 = %2d dB, Theoretical BER = %.2e\n', EbN0_dB(i), theoretical_ber(i));
    end
    
    %% 6. Plot Results
    figure;
    semilogy(EbN0_dB, avg_ber, 'bo-', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', 'Simulated BER');
    hold on;
    semilogy(EbN0_dB, theoretical_ber, 'r-', 'LineWidth', 2, 'DisplayName', 'Theoretical BER');
    grid on;
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('AWGN Channel Performance with Huffman Coding (BPSK)');
    legend('Location', 'best');
    axis([0 20 1e-6 1]);
    hold off;
    
    %% 7. Additional Analysis
    % Calculate coding gain at BER = 1e-4
    target_ber = 1e-4;
    [~, sim_idx] = min(abs(avg_ber - target_ber));
    [~, theo_idx] = min(abs(theoretical_ber - target_ber));
    
    if ~isempty(sim_idx) && ~isempty(theo_idx)
        coding_gain = EbN0_dB(theo_idx) - EbN0_dB(sim_idx);
        fprintf('\nCoding gain at BER ≈ %0.0e: %.2f dB\n', target_ber, coding_gain);
    end
    
    %% 8. Save BER Results to CSV
    T = table(EbN0_dB(:), avg_ber(:), theoretical_ber(:), ...
        'VariableNames', {'EbN0_dB', 'SimulatedBER', 'TheoreticalBER'});
    writetable(T, 'ber_huffman_bpsk.csv');
    fprintf('\nBER results saved to "ber_huffman_bpsk.csv"\n');
end
