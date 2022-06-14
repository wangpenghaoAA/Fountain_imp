clc;
clear;

%%%%%%%%%%%%%%%%%%%%parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LT_seed_bit = 10;
RS_parity_num = 2;
LT_k = 8;
LT_n = 34;
LT_c = 0.025;
LT_delta = 0.001;
LT_s = LT_c * sqrt(LT_k) * log(LT_k/LT_delta);
% RNG_input_seed = 10;
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

LT_Z = 0;
for i=1:LT_k
    LT_Z = LT_Z + soliton(i) + tau(i);
end

Robust_soliton = zeros(1,LT_k);
for i=1:LT_k
    Robust_soliton(i) = (soliton(i) + tau(i)) / LT_Z;
end

sum_test = zeros(1,LT_k);

for i=1:LT_k
    if(i == 1)
        sum_test(i) = Robust_soliton(1);
    else
        sum_test(i) = sum_test(i-1) + Robust_soliton(i);
    end
end
random_selection = [0	0	1	1	0	0	0	1	0	1	0	1	0	0	0	0	1	0	1	1	0	0	0	0	1	1	1	0	1	0	0	1	0	1	1	0	1	0	0	1	1	1	1	1	0	1	0	0	0	1	0	1	0	0	1	0	0	0	1	1	0	0	0	0	0	1	0	0	1	1	0	1	1	0	1	0	0	1	0	1	1	0	0	0	1	1	1	0	1	1	1	0	0	1	1	1	1	0	0	1	1	0	0	0	0	0	1	1	0	0	0	1	0	0	0	1	1	0	0	0	0	1	0	1	1	0	1	1	1	1	1	1	1	1	1	0	0	1	1	0	0	0	0	0	1	1	0	0	1	0	0	0	1	1	0	1	1	1	1	1	0	1	0	1	0	1	1	1	1	0	1	0	0	0	0	0	1	1	1	1	0	1	1	1	0	1	1	0	0	0	1	1	1	0	0	0	0	1	0	1	0	0	1	0	0	0	1	1	0	0	1	0	1	1	1	0	1	1	1	0	1	0	1	1	0	1	0	1	1	0	0	0	1	0	0	0	0	1	0	0	1	0	0	0	1	1	1	1	0	1	0	0	0	1	0	0	0	1	1	0	0	0	1	0	0	0	0	0	1	0	0	1	0	1	0	1	0	0	1	1	1	0	1	1	0	0	1	1	1	1	1	0	1	0	1	0	0	1	1	0	1	0	1	0	1	1	1	1	0	0	1	0	1	1	0	1	1	1	1	0	0	0	1	1	0	1	0	0	1	0	1	0	1	1	0	1	1	1	1	1	0	0	1	0	1	1	1	0	1	1	1	0];
for p = 1:size(random_selection,1)
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
        end
    for i = 1:size(decoding_bitstream,1)
    decoding_bitstream_all(p,LT_n*(i-1)+1:LT_n*i) = decoding_bitstream(i,:)
    end
end




   























