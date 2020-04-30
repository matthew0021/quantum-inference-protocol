function output = Cq_generate_concatenated(varargin)
%
% -------------------------------------------------------------------------
% Brief description: 
% 
% Cq_generate_concatenated outputs the inferred quantum statistical
% memory of a process at a chosen L value given N bits.
%
% This code was written as part of the Quantum Inference Project to explore
% how quantum mechanics can be used to study complexity. The code was written
% by Matthew Ho, a PhD candidate at the School of Physical and Mathematical
% Sciences at the Nanyang Technological University (NTU), Singapore & 
% Complexity Institute, NTU, Singapore, under the supervision
% of Dr. Thomas Elliott and Asst. Prof. Mile Gu.
%
% -------------------------------------------------------------------------
% Code description:
% 
% Outputs = Cq_generate_concatenated(L, bitstream)
% 
% Inputs: L = integer, history length -- Caution: Using too long an L will
%                                                 result in long
%                                                 computational time.
%         bits = bitstream (binary, for this function)
%
% Outputs:
% output{1,2} = Cq_concatenated;
% output{2,2} = sqrt_condFuture_probs;
% output{3,2} = stateProbVec;
% output{4,2} = condCounts;
% output{5,2} = stateProbCount;
% output{6,2} = rho;
% output{7,2} = topological_complexity;
% output{8,2} = sqrt(condProb);
% 
% -------------------------------------------------------------------------

% Log:
% Updated by Matthew Ho on 2018-02-18 at 22:40
% Updated by Matthew Ho on 2018-04-25 at 15:11
%   (Added topological complexity)
% Modified by Matthew Ho on 2018-06-21 at 00:30
%   (Replace bi2de with bin2dec, universal across all MATLAB versions)
% Modified by Matthew Ho on 2018-06-21 at 00:51 
%   (Only output Cq_concatenated
% Modified by Matthew Ho on 2018-06-21 at 01:05
%   (Removed genallpastfuturecombs. Put the function into script to save
%   hassle of adding 1 extra .m file to run.)
% Modified by Matthew Ho on 2019-01-07 at 16:51
%   (Edited description)
%   (Cleaned up commented areas of the code)
% Modified by Matthew Ho on 2020-03-02 at 14:13
%   (Included single time-step conditional probabilities as one of the
%   outputs)
% Modified by Matthew Ho on 2020-04-24 at 23:42
%   Added error check for number of arguments inputted. Same as
%   Unitary_generation.
%
% -------------------------------------------------------------------------



if nargin == 1
    % Check if it is L or bitstream that is inputted
    [~, C0] = size(varargin{1});
    if C0 > 1 % Bitstream is inputted
        bits = varargin{1};
        L = floor(log2(length(bits)/1000));
        if L > 8
            L = 10; % Lmax to cap at 10
        end
        if floor(log2(length(bits)/1000)) <= 0
            output = NaN;
            fprintf('Error, check bitstream length \n')
            return
        end
        
    elseif C0 == 1
        L = varargin{1};
        bits = NaN;
        output = NaN;
        fprintf('Error, check input for bitstream \n')
        return
    end
elseif nargin == 2
    % Check if 1st cell array is L, 2nd cell array is bits
    [~, C1] = size(varargin{1});
    [~, C2] = size(varargin{2});
    if C1 == 1 && C2 > 1
        L = varargin{1};
        bits = varargin{2};
        if floor(log2(length(bits)/1000)) <= 0
            output = NaN;
            fprintf('Error, check bitstream length \n')
            return
        end
        
    elseif C1 > 1 && C2 == 1
        L = varargin{2};
        bits = varargin{1};
        if floor(log2(length(bits)/1000)) <= 0
            output = NaN;
            fprintf('Error, check bitstream length \n')
            return
        end
        
    else
        output = NaN;
        fprintf('Error, check the inputs \n')
        return
    end
else
    output = NaN;
    fprintf('Error, check the inputs \n')
    return
end



%%%%%%%%%%%%%
%%% BITSTREAM 
%%%%%%%%%%%%%

% Flip it to a row vector if necessary.
[r, c] = size(bits);
if r > 1 && c == 1
    bits = bits';
end
bitsSize = size(bits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% USER INPUT LENGTH OF BITSTREAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

len = L; %length of bitstream

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TO CONSTRUCT MATRIX WITH INCREASING TIMESTEPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bitsMatrix = zeros(bitsSize(2)-(2*len)+len,(len+1)); % Edited on 24-Mar-2018

for i=1:1:bitsSize(2)-(2*len)+len  % Edited on 24-Mar-2018
    bitsMatrix(i,:) = bits(i:(len+1)+i-1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CONDITIONAL PROBABILITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condCounts = zeros(2^len,2);
for i=1:1:size(bitsMatrix,1)
    r = bin2dec(num2str(bitsMatrix(i,1:len)))+1;
    c = bitsMatrix(i,len+1)+1;
    condCounts(r,c) = condCounts(r,c)+1;
end

condCounts;
condProb = zeros(2^len,2);
for i=1:1:2^len
    for j=1:1:2
        condProb(i,j) = condCounts(i,j)/sum(condCounts(i,:));
    end
end

for i=1:1:2^len
    for j=1:1:2
        if isnan(condProb(i,j)) == 1
            condProb(i,j) = 0;
        end
    end
end
condProb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TO GENERATE MATRIX OF SIZE (2^L, L)
%%% EACH ROW CORRESPONDS TO THE BINARY 
%%% REPRESENTATION OF EACH ROW NUMBER+1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generating the past and future combinations
base = 2;
maxdecimal = base^len;
nums = linspace(1,maxdecimal,maxdecimal);
past = [];
for i=1:1:length(nums)
    str = dec2base(nums(i)-1,base);
    vector = str -'0';
    while length(vector)<len
        vector = [0 vector];
    end
    past = vertcat(past,vector);
end
past;
future = past;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EVALUATING CONDITIONAL PROBABILITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condfut_no_need_to_normalise = zeros(2^len,2^len);
%tic
for i=1:1:2^len
    for j=1:1:2^len
        
        % step (1)
        temp = [past(i,:) future(j,:)];
        
        % step (2)
        mat = zeros(len,len+1);
        for k=1:1:len
            mat(k,:) = temp(k:k+len);
        end
        mat;
        
        % step (3)
        for k=1:1:len
            row_to_access = bin2dec(num2str(mat(k,1:len)))+1;
            col_to_access = mat(k,len+1)+1;
            prob_from_mat(k) = condProb(row_to_access,col_to_access);
        end
        
        condfut_no_need_to_normalise(i,j) = prod(prob_from_mat);
        
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EVALUATING CONDITIONAL PROBABILITIES 2^L BY 2^L DIRECTLY
%%% NOT CONCATENATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bitsSize_2L2L = size(bits);
bitsMatrix_2L2L = zeros(bitsSize_2L2L(2)-2*len+1,2*len);
for i=1:1:length(bits)-2*len+1
    bitsMatrix_2L2L(i,:) = bits(i:2*len+i-1);
end
condFuture_counts_2L2L = zeros(2^len,2^len);
condFuture_probs_2L2L = zeros(2^len,2^len);

for h=1:1:size(bitsMatrix_2L2L,1)
    r = bin2dec(num2str(bitsMatrix_2L2L(h,1:len)))+1;
    c = bin2dec(num2str(bitsMatrix_2L2L(h,len+1:2*len)))+1;
    condFuture_counts_2L2L(r,c) = condFuture_counts_2L2L(r,c)+1;    
end

condFuture_probs_2L2L = zeros(2^len,2^len);
for i=1:1:2^len
    condFuture_probs_2L2L(i,:) = condFuture_counts_2L2L(i,:)/sum(condFuture_counts_2L2L(i,:));
end
condFuture_probs_2L2L(isnan(condFuture_probs_2L2L)) = 0;

condFuture_counts_2L2L;
condFuture_probs_2L2L;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SQRT THE CONDITIONAL PROBABILITIES FOR QUANTUM STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sqrt_condFuture_probs = sqrt(condfut_no_need_to_normalise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FINDING PROBABILITY OF EACH STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bitsMatrix2 = zeros(bitsSize(2)-len+1,len);
for i=1:1:bitsSize(2)-len+1
    bitsMatrix2(i,:) = bits(i:len+i-1);
end
[r6, ~] = size(bitsMatrix2);

stateProbCount = zeros(2^len,1);
for h=1:1:r6
    r = bin2dec(num2str(bitsMatrix2(h,:))) + 1;
    stateProbCount(r) = stateProbCount(r)+1;
end

% Calculating probability vector
stateProbCount;
stateProbVec = stateProbCount/sum(stateProbCount);

%%%%%%%%%%%%%%%%
%%% FINDING \rho
%%%%%%%%%%%%%%%%

rho = zeros(2^len,2^len);
for i=1:1:2^len
    rho = rho + stateProbVec(i) * sqrt_condFuture_probs(i,:)' * sqrt_condFuture_probs(i,:);
end

rho;

% Calculating Cq using initial \rho
[~, eigenvalues] = eig(rho);
log2eigenvalues = real(log2(eigenvalues));  

eigenvalues;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CALCULATING INITIAL CQ
%%%%%%%%%%%%%%%%%%%%%%%%%%

[r4, ~] = size(log2eigenvalues);
for i=1:1:r4
    for j=1:1:r4
        if log2eigenvalues(i,j) == -Inf
            log2eigenvalues(i,j) = 0;
        end
    end
end

Cq_concatenated = abs(-trace(eigenvalues*log2eigenvalues));

count = 0;
for i=1:1:2^len
    if eigenvalues(i,i) ~= 0
        count = count + 1;
    end
end
count;

% For dimensions of the system, use topological complexity
%topological_complexity = log2(count);
topological_complexity = log2(2^len);

%output = Cq_concatenated;

output{1,1} = sprintf('Inferred Cq(L=%d)',L);
output{2,1} = sprintf('Prob amplitudes of 2^%d=%d QMS',L,2^L);
output{3,1} = sprintf('Prob of each QMS');
output{4,1} = sprintf('Counts for prob amplitudes for 2^%d=%d QMS',L,2^L);
output{5,1} = sprintf('Counts for probs of each QMS');
output{6,1} = sprintf('rho: density matrix of QMS 2^%d by 2^%d',L,L);
output{7,1} = sprintf('Topological copmlexity Dq = log2(dim(S))');
output{8,1} = sprintf('Prob amplitudes, single time-step');
output{9,1} = sprintf('Prob amplitudes of 2^%d=%d QMS, not concatenated',L,2^L);
output{10,1} = sprintf('Condcounts 2^%d by 2^%d, not concatenated',L,L);

output{1,2} = Cq_concatenated;
output{2,2} = sqrt_condFuture_probs;
output{3,2} = stateProbVec;
output{4,2} = condCounts;
output{5,2} = stateProbCount;
output{6,2} = rho;
output{7,2} = topological_complexity;
output{8,2} = sqrt(condProb);
output{9,2} = sqrt(condFuture_probs_2L2L);
output{10,2} = condFuture_counts_2L2L;



end

