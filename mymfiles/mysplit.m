function [lines] = mysplit(in)
% MYSPLIT: split a string by newlines.

cells = {};
idx = 1;
while idx < length(in)
  tok = strtok(in(idx:end), char(10));
  cells{length(cells) + 1} = tok;
  idx = idx + length(tok) + 1;
end

lines = char(cells);

