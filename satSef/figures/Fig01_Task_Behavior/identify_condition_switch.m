function [ trial_switch ] = identify_condition_switch( behavInfo )
%identify_condition_switch Summary of this function goes here

DEBUG = false;

NUM_SESSION = size(behavInfo,1);

FLAG_ACC = 1;
FLAG_NORM = 2;
FLAG_FAST = 3;

trial_switch = new_struct({'A2F','F2A','A2N','N2F'}, 'dim',[1,NUM_SESSION]);

for kk = 1:NUM_SESSION
  
  num_trials = length(behavInfo.condition{kk});
  condition = double(behavInfo.condition{kk});
  condition(condition == 0) = NaN;
  
  %identify all trials with condition switch
  tmp_F2A = find(diff(condition) == (FLAG_ACC - FLAG_FAST)) + 1 ;
  tmp_A2F = find(diff(condition) == (FLAG_FAST - FLAG_ACC)) + 1 ;
  
  if ismember(behavInfo.session{kk}(1), {'Q','S'})
    tmp_A2N = find((condition(1:end-1) == FLAG_ACC) & (diff(condition) == 1)) + 1;
    tmp_N2F = find((condition(1:end-1) == FLAG_NORM) & (diff(condition) == 1)) + 1;
  end
  
  %remove those trials closest to session end
  tmp_F2A(tmp_F2A > (num_trials-4)) = [];
  tmp_A2F(tmp_A2F > (num_trials-4)) = [];
  
  if (DEBUG)
    figure(); hold on
    plot(condition, 'ko-')
  end%if:DEBUG
  
  trial_switch(kk).F2A = tmp_F2A;
  trial_switch(kk).A2F = tmp_A2F;
  
  if ismember(behavInfo.session{kk}(1), {'Q','S'})
    tmp_N2F(tmp_N2F > (num_trials-9)) = [];
    tmp_A2N(tmp_A2N > (num_trials-9)) = [];
    trial_switch(kk).A2N = tmp_A2N;
    trial_switch(kk).N2F = tmp_N2F;
  end
  
end%for:sessions(kk)

trial_switch = struct2table(trial_switch);

end%function:identify_condition_switch()

