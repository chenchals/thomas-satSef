function [ ] = anova_TwoWay_Between_SAT( Y , Condition , Efficiency , varargin )
%anova_TwoWay_Between_SAT Summary of this function goes here
%   Input
%     'display' -- {'on','off'} -- Show ANOVA table as output?
%     'model' -- {'linear','interaction','full'}
%     'sstype' -- {'1','2','3'} -- Form of computation of sum of squares
% 


% Modification:
%   21/04/2020 - removed getopt and using inputParser to parse argument
%                list

%args = getopt(varargin, {{'display=','on'}, {'model=','full'}, {'sstype=',3}});

argParser = inputParser();
argParser.addParameter('display','on');
argParser.addParameter('model','full');
argParser.addParameter('sstype',3);
argParser.parse(varargin{:});
args = argParser.Results;

[~,tbl] = anovan(Y, {Condition Efficiency}, ...
  'display',args.display, ...
  'model',args.model, ...
  'sstype',args.sstype, ...
  'varnames',{'Condition','Efficiency'});

end % fxn : anova_TwoWay_Between_SAT()

