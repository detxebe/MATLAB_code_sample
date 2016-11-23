function [NPV, NPVs]= Fx_N(x0_N, ro, cwp, cwi, d, NIs, NPs, Nt, T, dt, xd, xu);
%
%      [NPV, NPVs] = function Fx_(x0_N, ro, cwp, cwi, d, NIs, NPs, Nt, T, dt, xd, xu);
%
% Determines discounted NPV (total and as a function of time) for the
% reservoir RES1 (line drive with 4 injectors and 4 producers) for
% production time frame of 3000 days and BHP controls given by x0.
%
% x0_N ..... vector of normalized [0,1] Nt x (NIs + NPs) well controls; 
%            the order is first injectors for first control interval, then 
%            producers for first control interval, then injectors for 
%            second control interval, etc.
% ro ....... oil sale price (USD/bbl)
% cwp ...... water separation cost (USD/bbl)
% cwi ...... water injection cost (USD/bbl)
% d ........ discount factor
% NIs ...... number of injectors
% NPs ...... number of producers
% Nt ....... number of control intervals per well
% T ........ total production time frame (has to be consistent with control.in)
% dt ....... time interval in model response NPV (has to be consistent with control.in)
% xd ....... vector of Nt x (NIs + NPs) lower bound controls; used for denormalization
% xu ....... vector of Nt x (NIs + NPs) upper bound controls; used for denormalization
%
% NPV ...... total discounted NPV in USD million (= sum(NPVs))
% NPVs ..... sequence of partial NPVs (in USD million) with time interval dt
%



% 0. DENORMALIZATION

x0_N = x0_N(:);     % be sure it is a column vector
xd   = xd(:);
xu   = xu(:);

x0   = xd + x0_N.*(xu-xd);


% 1.A. CREATE SIMULATION FOLDER

warning off;

curr_dir = pwd;
mkdir('sim');
cd('sim');


% 1.B COPY BASIC SIMULATION FILES + GPRS SIMULATOR

basic_sim_folder = 'gprs_model';

if (ispc)
   [kk01, kk02, kk03] = copyfile([curr_dir filesep basic_sim_folder filesep 'windows' filesep '*.*'], '.', 'f');
else    
   [kk01, kk02, kk03] = copyfile([curr_dir filesep basic_sim_folder filesep 'linux'   filesep '*.*'], '.', 'f');
end;


% 1.C WRITE INPUT FILES

x0 = reshape(x0, (NIs+NPs), Nt);

ts = 0:round(T/Nt):T;     % different time interval starts
if (ts(end)~=T)
    ts(end) = T;
end    

for i=1:NIs
  
    file_name = sprintf('WellINJ%d.in', i);
    
    fid       = fopen(file_name, 'w');
    for j=1:Nt 
       fprintf(fid, 'TBHP	     %7.2f	     %7.2f	     %7.2f	  2  0  1\n', ...
               ts(j), ts(j+1), x0(i,j)); 
    end    
    fclose(fid);    
    
end

for i=1:NPs
  
    file_name = sprintf('WellPROD%d.in', i);
    
    fid       = fopen(file_name, 'w');
    for j=1:Nt 
       fprintf(fid, 'TBHP	     %7.2f	     %7.2f	     %7.2f\n', ...
               ts(j), ts(j+1), x0(NIs+i,j)); 
    end    
    fclose(fid);    
    
end


% 2. RUN SIMULATION

if (ispc)
    [kk04, kk05] = system('gprs_win');     % recheck particular sintax
else
    [kk04, kk05] = system('gprs_lnx');     % recheck particular sintax
end


% 3. READ OUTPUT FILES

name_res = 'RES1';

time_dts = (dt:dt:T)';
if (time_dts(end)~=T)
   time_dts(end+1) = T;
end   
Ldts     = length(time_dts);

WWIRs    = zeros(Ldts, NIs);
WWPRs    = zeros(Ldts, NPs);
WOPRs    = zeros(Ldts, NPs);

for i=1:NIs
  
    file_name = sprintf([name_res '_INJ%d.out'], i);
    
    fid       = fopen(file_name);
    
    for j=1:5, kk = fgets(fid); end    
        
    kk = fscanf(fid, '%f'); 
    
    fclose(fid);    
    
    if (ispc)     % 4 COLUMNS (first one time, last one water rate)
       
       kk = reshape(kk, 4, length(kk)/4)';     
        
    else          % 5 COLUMNS (first one time, last one water rate)
        
       kk = reshape(kk, 5, length(kk)/5)'; 
       kk = [zeros(1,5);kk]; 
       
    end
    
    % INTERPOLATE EVERY dt DAYS
    
    % ---------
    % AVOID REPEATING ROWS FOR TIME COLUMN
    [kk01, Ikk] = unique(kk(:,1));
    kk          = kk(Ikk, :);
    % ---------
    
    WWIRs(:,i) = interp1(kk(:,1), -kk(:, end), time_dts);
        
end

for i=1:NPs
  
    file_name = sprintf([name_res '_PROD%d.out'], i);
    
    fid       = fopen(file_name);
    
    for j=1:5, kk = fgets(fid); end    
        
    kk = fscanf(fid, '%f'); 
    
    fclose(fid);    
    
    if (ispc)     % 4 COLUMNS (first one time, last but one oil rate, last one water rate)
       
       kk = reshape(kk, 4, length(kk)/4)';     
        
    else          % 5 COLUMNS (first one time, last but one oil rate, last one water rate)
        
       kk = reshape(kk, 5, length(kk)/5)'; 
       kk = [zeros(1,5);kk]; 
       
    end
    
    % INTERPOLATE EVERY dt DAYS

    % ---------
    % AVOID REPEATING ROWS FOR TIME COLUMN
    [kk01, Ikk] = unique(kk(:,1));
    kk          = kk(Ikk, :);
    % ---------
    
    WOPRs(:,i) = interp1(kk(:,1),  kk(:, end-1), time_dts);
    WWPRs(:,i) = interp1(kk(:,1),  kk(:, end  ), time_dts);
        
end


% 4. COMPUTE DISCOUNTED NPV "rate" IN USD MILLION (NPV total is just the sum)
    
FWIRs      = sum(WWIRs, 2); 
FOPRs      = sum(WOPRs, 2); 
FWPRs      = sum(WWPRs, 2); 

NPVs       = diff([0; time_dts]).*(FOPRs * ro - FWPRs * cwp - FWIRs * cwi);
disc_facts = (1/(1+d)).^(time_dts/365.25);       % time expressed in years for discount

NPVs       = NPVs.*disc_facts/1e6;               % discounted NPV in USD million 

NPV        = sum(NPVs);


% C. CLEAN FILES

cd(curr_dir);

rmdir('sim', 's');