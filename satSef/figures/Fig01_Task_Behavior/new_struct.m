function [ S ] = new_struct( fields , varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

args = getopt(varargin, {{'dim=',1}});

%check inputs
if ~isnumeric(args.dim)
  error('Input "dim" must be of type numeric.')
end
if isscalar(args.dim)
  dim = [1,args.dim];
else
  dim = args.dim;
end

%initialize new struct
S(dim(1),dim(2)).(fields{1}) = [];

%create empty fields
N_fields = length(fields);
for f = 2:N_fields
  S(1).(fields{f}) = [];
end

end
