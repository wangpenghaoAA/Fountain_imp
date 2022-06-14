clear;
clc;
clear;

% imdata = imread('Jetplane.jpg');
% msgBin = de2bi(imdata,8,'left-msb');
% len = size(msgBin,1).*size(msgBin,2);
% input_data_save = reshape(double(msgBin).',len,1).';
fid=fopen('poems.txt');  %有12首诗
s=fread(fid,'*char');
fclose(fid);
unis=unique(s,'stable');
SASCII=abs(s)';
msgBin = de2bi(int8(s),8,'left-msb');
len = size(msgBin,1).*size(msgBin,2);
input_data_save = reshape(double(msgBin).',len,1).';




LT_L = 272;

LT_N = 1126;
% input_filename = sprintf('Jetplane.txt');
% [FP] = fopen(input_filename,'r');
% 
% input_data_save = fscanf(FP,'%d');
% fclose(FP);

input_file_bit_size = size(input_data_save,1);
Binary_Data_input = zeros(LT_N,LT_L);
for i=1:LT_N
    for j=1:LT_L
        if((i-1)*LT_L+j<input_file_bit_size + 1)
            Binary_Data_input(i,j)=input_data_save((i-1)*LT_L+j);
        else
            Binary_Data_input(i,j)=randi(2)-1;
        end
            
    end
end
p=1;
while p<=LT_N
    Binary_Data_input_n = zeros(1,LT_L);
    Binary_Data_input_n = Binary_Data_input(p,:);
    LT_n = 34;
    LT_k = 8;
    LT_seed_bit = 10;
    group_n =8;
    LT_c = 0.025;
    LT_delta = 0.001;
    LT_s = LT_c * sqrt(LT_k) * log(LT_k/LT_delta);
    for i=1:LT_k
       for j=1:LT_n
          if((i-1)*LT_n+j<LT_L + 1)
            Binary_Data_input_n_n(i,j)=Binary_Data_input_n((i-1)*LT_n+j);
          end 
       end
    end


    Max_seed_num = power(2,LT_seed_bit)-1;
    LT_generator = zeros(LT_k+1,LT_n); % seed at the very last row. messages that were used in encoding for each columns.
    LFSR_initial_seed = 42;
    RNG_input_seed = 10;

    soliton = zeros(1,LT_k);
    for i=1:LT_k
        if(i==1)
            soliton(i) = 1/LT_k;
        end

        if(i~=1)
            soliton(i) = 1/i/(i-1);
        end
    end

    tau = zeros(1,LT_k);
    KSratio = round(LT_k / LT_s);
    for i=1:LT_k
        if(i == KSratio)
            tau(i) = LT_s * log(LT_s/LT_delta) / LT_k;
        end

        if(i<KSratio)
            tau(i) = LT_s / LT_k / i;
        end
    end

    LT_z = 0;
    for i=1:LT_k
        LT_z = LT_z + soliton(i) + tau(i);
    end

    Robust_soliton = zeros(1,LT_k);
    for i=1:LT_k
        Robust_soliton(i) = (soliton(i) + tau(i)) / LT_z;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sum_test = zeros(1,LT_k);

    for i=1:LT_k
        if(i == 1)
            sum_test(i) = Robust_soliton(1);
        else
            sum_test(i) = sum_test(i-1) + Robust_soliton(i);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    total_degree_sum = 0;

    %%%%%%%%%%%%%%%%%%%%%%RS parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    RS_field_size = 8;
    RS_parity_num = 2;    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    rng(RNG_input_seed);
    %%%%%%%%%%%%%%%%%%%%% encoding trial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    try_encoding = 0;
    success_encoding = 0;
    encoding_max_number = Max_seed_num;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%% homopolymer / GC content%%%%%%%%%%%%%%%%%%%%

    Max_run_length = 3;
    Min_GC_content = 0.45;
    MAX_GC_content = 0.55;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while(try_encoding < encoding_max_number)
        total_degree_sum_test = 0;
        success_GC = 1;
        success_run = 1;

        if(try_encoding == 0)
            temp1_seed = de2bi(LFSR_initial_seed,LT_seed_bit);
            Galois_register = temp1_seed;
        end

        temp2_seed = Galois_register;
        for j=1:LT_seed_bit
            if(j~=4 &&j~=5&&j~=8&&j~=10)
               Galois_register(j) = temp2_seed(j+1); 
            end

            if(j==10)
               Galois_register(j) = temp2_seed(1);
            end

            if(j==4||j==5||j==8)
               Galois_register(j) = mod(temp2_seed(j+1) + temp2_seed(1),2);
            end

        end

        
        %Here, seed is stored at Galois_register. Then, we need to make LT code one symbol.

        Galois_register_dec_temp = bi2de(Galois_register);
        rng(Galois_register_dec_temp);
        degree_prop = rand;
        LT_generator_temp = zeros(LT_k+1,1);


        for j=1:LT_k
            if(j==1)
                if(degree_prop < Robust_soliton(1))
                    symbolselection = datasample([1:LT_k],1,'Replace',false);
                    LT_generator_temp(symbolselection(1),1)=1;
                    total_degree_sum_test = total_degree_sum_test + 1;
                end
            end

            if(j~=1)
                if((sum_test(j-1)<=degree_prop) && (degree_prop < sum_test(j)))
                    symbolselection = datasample([1:LT_k],j,'Replace',false);
                    total_degree_sum_test = total_degree_sum_test + j;
                    for m=1:j
                        LT_generator_temp(symbolselection(m),1)=1;
                    end
                end
            end

            if(degree_prop==1)
                LT_generator_temp(:,1)=1;
                total_degree_sum_test = total_degree_sum_test + LT_k;
            end
        end

        LT_generator_temp(LT_k+1,1) = Galois_register_dec_temp;

        %Here, LT code one symbol is generated.


        Galois_seed_temp = fliplr(Galois_register);
        One_strand_bit_temp(1:LT_seed_bit) = Galois_seed_temp;
        LT_code_one_symbol_encoding_temp = LT_generator_temp(1:LT_k);
        LT_code_one_symbol_temp = mod(transpose(LT_code_one_symbol_encoding_temp)*Binary_Data_input_n_n,2);
        One_strand_bit_temp(LT_seed_bit+1:LT_seed_bit+LT_n) = LT_code_one_symbol_temp;

        % Seed and LT code is done. RS encoding is needed.
        % Changing into ACGT is needed.

        for i=1:(LT_seed_bit+LT_n)/2
            if((One_strand_bit_temp(2*i-1)==0) && (One_strand_bit_temp(2*i)==0))
                One_strand_ACGT_temp(i) = 'A';
            end

            if((One_strand_bit_temp(2*i-1)==0) && (One_strand_bit_temp(2*i)==1))
                One_strand_ACGT_temp(i) = 'C';
            end

            if((One_strand_bit_temp(2*i-1)==1) && (One_strand_bit_temp(2*i)==0))
                One_strand_ACGT_temp(i) = 'G';
            end

            if((One_strand_bit_temp(2*i-1)==1) && (One_strand_bit_temp(2*i)==1))
                One_strand_ACGT_temp(i) = 'T';
            end 
        end

        %One_strand_ACGT_temp is generated. size is (1, x)
        %length is (LT_seed_bit+LT_L+8*RS_parity_num)/2.

        for i=1:(LT_seed_bit+LT_n)/2-Max_run_length
            for j=1:Max_run_length
                if(One_strand_ACGT_temp(i)~=One_strand_ACGT_temp(i+j))
                    break;
                end

                if(j==Max_run_length)
                    success_run = 0;
                end
            end

            if(success_run == 0)
                break;
            end
        end  

        GC_count = 0;
        for i=1:(LT_seed_bit+LT_n)/2
            if((One_strand_ACGT_temp(i)=='G') || (One_strand_ACGT_temp(i)=='C'))
                GC_count = GC_count + 1;
            end
        end

        GC_ratio = GC_count / ((LT_seed_bit+LT_n)/2);

        if((GC_ratio < Min_GC_content) || (GC_ratio > MAX_GC_content))
            success_GC = 0;
        end


        if((success_run == 1) && (success_GC == 1))
            success_encoding = success_encoding + 1;

            Binary_Data_output(p,(LT_seed_bit+LT_n)*(success_encoding-1)+1:success_encoding*(LT_seed_bit+LT_n)) = One_strand_bit_temp;
            ACGT_Data_output(p,((LT_seed_bit+LT_n)/2)*(success_encoding-1)+1:success_encoding*((LT_seed_bit+LT_n)/2)) = One_strand_ACGT_temp
     
        end
        try_encoding = try_encoding + 1;

    end
    L = 0;D=0;
    while L~=LT_n
        random_length = success_encoding;
        rng(D)
        random_num = randperm(random_length,group_n);
        for t = 1:group_n
            random_selection((t-1)*(LT_seed_bit+LT_n)+1:t*(LT_seed_bit+LT_n)) = Binary_Data_output(p,(random_num(t)-1)*(LT_seed_bit+LT_n)+1:random_num(t)*(LT_seed_bit+LT_n));
        end
        random__output_selection = random_selection;
        decoding_flag = zeros(1,LT_k);  
        for k = 1: size(random_selection ,2)/(LT_seed_bit+LT_n)
        temp_vector_candidate = random_selection((LT_seed_bit+LT_n)*(k-1)+1:(LT_seed_bit+LT_n)*(k-1)+10);
        Galois_register = bi2de(temp_vector_candidate, 'left-msb');
        rng(Galois_register);
        degree_prop = rand;
        LT_generator_temp = zeros(LT_k,1); 
        for j=1:LT_k
            if(j==1)
            if(degree_prop < Robust_soliton(1))
                symbolselection = datasample([1:LT_k],1,'Replace',false);
                LT_generator_temp(symbolselection(1),1)=1;
            end
            end

            if(j~=1)
                if((sum_test(j-1)<=degree_prop) && (degree_prop < sum_test(j)))    
                symbolselection = datasample([1:LT_k],j,'Replace',false);
                    for c=1:j
                    LT_generator_temp(symbolselection(c),1)=1;
                    end
                end
            end
        end

        if(degree_prop==1)
        LT_generator_temp(:,1)=1;
        end
        received_xor_state(k,:) = transpose(LT_generator_temp); 
        end
        decoding_flag = zeros(1,LT_k);
        decoding_bitstream =  zeros(LT_k,LT_n);

        for i = 1:size(received_xor_state,1)-1
            for j = i+1:size(received_xor_state,1)    
                 inferred_index = mod(received_xor_state(i,:) + received_xor_state(j,:),2) ; 
                 subtracted_sequence1 = random_selection((LT_seed_bit+LT_n)*(i-1)+11:(LT_seed_bit+LT_n)*i);
                 subtracted_sequence2 = random_selection((LT_seed_bit+LT_n)*(j-1)+11:(LT_seed_bit+LT_n)*j);
                 subtracted_sequence3 = mod(subtracted_sequence1+subtracted_sequence2,2);
                 A=0;B=0;C=0;
                 a=1;b=1;c=1;
                 for a = 1:length(received_xor_state(i,:))
                     A = A + received_xor_state(i,a);
                 end
                 for b =1:length(received_xor_state(j,:))
                     B = B + received_xor_state(j,b);
                 end
                 for c = 1:length(inferred_index)
                     C = C + inferred_index(c);
                 end
                 if C ==1||C == 0
                     if C ==1
                        a = find(inferred_index);
                        received_xor_state(j,:) = inferred_index;
                        decoding_bitstream(a,:) = subtracted_sequence3;
                        random_selection((LT_seed_bit+LT_n)*(j-1)+11:(LT_seed_bit+LT_n)*j) = subtracted_sequence3;
                        decoding_flag(1,a) = 1;
                        for z = 1:size(received_xor_state,1)
                            if received_xor_state(z,a)==1
                                subtracted_sequence4 = random_selection((LT_seed_bit+LT_n)*(z-1)+11:(LT_seed_bit+LT_n)*z);
                                received_xor_state(z,:) = mod(received_xor_state(z,:) + inferred_index ,2) ;
                                random_selection((LT_seed_bit+LT_n)*(z-1)+11:(LT_seed_bit+LT_n)*z) = mod(decoding_bitstream(a,:)+subtracted_sequence4,2);
                            else
                            end
                        end
                        z = 1;
                        while z <= size(received_xor_state,1)
                            number0 = 0;
                            for q = 1:size(received_xor_state,2)
                                number0 = number0 + received_xor_state(z,q);
                            end
                            if number0 == 1
                               a = find(received_xor_state(z,:)) ;
                               inferred_index0 = received_xor_state(z,:);
                               decoding_bitstream(a,:) = random_selection((LT_seed_bit+LT_n)*(z-1)+11:(LT_seed_bit+LT_n)*z);
                               decoding_flag(1,a) = 1;
                               for r = 1:size(received_xor_state,1)
                                   if received_xor_state(r,a)==1
                                      received_xor_state(r,:) = mod(received_xor_state(r,:) + inferred_index0,2);
                                      subtracted_sequence0 = random_selection((LT_seed_bit+LT_n)*(r-1)+11:(LT_seed_bit+LT_n)*r);
                                      random_selection((LT_seed_bit+LT_n)*(r-1)+11:(LT_seed_bit+LT_n)*r) = mod(decoding_bitstream(a,:)+subtracted_sequence0,2);
                                   end
                               end
                            z  = 1;
                            continue;
                            end
                            z = z+1;
                         end
                     else
                      received_xor_state(j,:) = inferred_index;
                     end
                 else
                    if C<B||C<A
                    received_xor_state(j,:) = inferred_index; 
                    random_selection((LT_seed_bit+LT_n)*(j-1)+11:(LT_seed_bit+LT_n)*j) = subtracted_sequence3;
                    end
                 end
            end
            L = 0;
            for i = 1:size(decoding_flag,2)
            L = L + decoding_flag(i);
            end
            if L == LT_k
                break;
            end
        end
        if L == LT_k
            Before_site_selection_bit(p,:) = random__output_selection
            break;
        end
        D = D+1;
    end
    p=p+1;
end
for q= 1:LT_N
    for i=1:length(Before_site_selection_bit (q,:))/(LT_seed_bit+LT_n)
            for j=1:LT_seed_bit+LT_n
                Before_site_selection_bit_temp(j)=Before_site_selection_bit(q,(LT_seed_bit+LT_n)*(i-1)+j);
            end
            Before_site_selection_temp(q,i) = bi2de(Before_site_selection_bit_temp);
    end
end
for i = 1:LT_N
   Index_tem(i) = Before_site_selection_temp(i,1);
   if i>=2
      for j = 1:i-1
          while Index_tem(i)==Index_tem(i-j)
              for m=2:group_n 
                if Index_tem(i) ~= Before_site_selection_temp(i,m)
                   k = Before_site_selection_temp(i,m);
                   Before_site_selection_temp(i,m)=Index_tem(i);
                   Index_tem(i)=k;
                   break;
                end
              end
          end
      end
   end
end
After_site_selection_temp = zeros (LT_N,group_n);
for i = 1:LT_N
    After_site_selection_temp(i,1) = Index_tem(i);
    for j = 2:group_n
       After_site_selection_temp(i,j) = Before_site_selection_temp(i,j);
    end
end
for i = 1:LT_N
    for j = 1:group_n
         After_site_selection_bit_tem(i,(LT_seed_bit+LT_n)*(j-1)+1:(LT_seed_bit+LT_n)*j) = de2bi(After_site_selection_temp(i,j),(LT_seed_bit+LT_n));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For every 8 bits, correspond to GF(256) and encode them.
for m = 1:size( After_site_selection_bit_tem, 1)
    Binary_Data_output_size =length( After_site_selection_bit_tem(m,:));
    RS_input_temp = zeros(1,Binary_Data_output_size/8);
    eightbit_temp = zeros(1,8);
    for i=1:Binary_Data_output_size/8
        for j=1:8
            eightbit_temp(j)= After_site_selection_bit_tem(m,(8*i-j+1));%倒过来
        end
        RS_input_temp(i) = bi2de(eightbit_temp);
    end
    msg_temp = gf(RS_input_temp,8);
    RS_code_temp = rsenc(msg_temp,(Binary_Data_output_size+8*RS_parity_num)/8,Binary_Data_output_size/8);

    RS_code_parity_binary_temp = zeros(1,8*RS_parity_num);
    eightbit_temp2 = zeros(1,8);
    galois_test = RS_code_temp.x;

    for i=1:RS_parity_num
        eightbit_temp2 = de2bi(galois_test(i),8);
        RS_bit(m,(i-1)*8+1:8*i)=fliplr(eightbit_temp2);
    end
    
end
data_RS_bit = [After_site_selection_bit_tem RS_bit];
output_file_name = sprintf('Encoded_oligo_5.9.txt');
[FP_out] = fopen(output_file_name, 'wt');
for i=1:size(data_RS_bit,1)
    for j=1:size(data_RS_bit,2)/2
        if((data_RS_bit(i,(2*j-1))==0) && (data_RS_bit(i,(2*j))==0))
            fprintf(FP_out, 'A');
        elseif((data_RS_bit(i,(2*j-1))==0) && (data_RS_bit(i,(2*j))==1))
            fprintf(FP_out, 'C');
        elseif((data_RS_bit(i,(2*j-1))==1) && (data_RS_bit(i,(2*j))==0))
            fprintf(FP_out, 'G');
        elseif((data_RS_bit(i,(2*j-1))==1) && (data_RS_bit(i,(2*j))==1))
            fprintf(FP_out, 'T');
        end
    end
    fprintf(FP_out, '\n');
end