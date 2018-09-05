
function x = mnrnd (n, p, s)

%   % Check arguments
%   if (nargin == 3)
%     if (~ isscalar (n) || n < 0 || round (n) ~= n)
%       error ('mnrnd: n must be a non-negative integer');
%     end
%     if (~ isvector (p) || any (p < 0 | p > 1))
%       error ('mnrnd: p must be a vector of probabilities');
%     end
%     if (~ isscalar (s) || s < 0 || round (s) ~= s)
%       error ('mnrnd: s must be a non-negative integer');
%     end
%   elseif (nargin == 2)
%     if (isvector (p) && size (p, 1) > 1)
%       p = p';
%     end
%     if (~ isvector (n) || any (n < 0 | round (n) ~= n) || size (n, 2) > 1)
%       error ('mnrnd: n must be a non-negative integer column vector');
%     end
%     if (~ ismatrix (p) | isempty (p) | any (p < 0 | p > 1))
%       error ('mnrnd: p must be a non-empty matrix with rows of probabilities');
%     end
%     if (~ isscalar (n) && size (p, 1) > 1 && length (n) ~= size (p, 1))
%       error ('mnrnd: the length of n must match the number of rows of p');
%     end
%   else
%     error ('print usage');
%   end

  % Adjust input sizes
  if (nargin == 3)
    n = n * ones (s, 1);
    p = repmat (p(:)', s, 1);
  elseif (nargin == 2)
    if (isscalar (n) && size (p, 1) > 1)
      n = n * ones (size (p, 1), 1);
    elseif (size (p, 1) == 1)
      p = repmat (p, length (n), 1);
    end
  end
  sz = size (p);

  % Upper bounds of categories
  ub = cumsum (p, 2);
  % Make sure that the greatest upper bound is 1
  gub = ub(:, end);
  ub(:, end) = 1;
  % Lower bounds of categories
  lb = [zeros(sz(1), 1) ub(:, 1:(end-1))];

  % Draw multinomial samples
  x = zeros (sz);
  for i = 1:sz(1)
    % Draw uniform random numbers
    r = repmat (rand (n(i), 1), 1, sz(2));
    % Compare the random numbers of r to the cumulated probabilities of p and
    % count the number of samples for each category
    x(i, :) =  sum (r <= repmat (ub(i, :), n(i), 1) & r > repmat (lb(i, :), n(i), 1), 1);
  end
  % Set invalid rows to NaN
  k = (abs (gub - 1) > 1e-6);
  x(k, :) = NaN;

end

%~test
%~ n = 10;
%~ p = [0.2, 0.5, 0.3];
%~ x = mnrnd (n, p);
%~ assert (size (x), size (p));
%~ assert (all (x >= 0));
%~ assert (all (round (x) == x));
%~ assert (sum (x) == n);

%~test
%~ n = 10 * ones (3, 1);
%~ p = [0.2, 0.5, 0.3];
%~ x = mnrnd (n, p);
%~ assert (size (x), [length(n), length(p)]);
%~ assert (all (x >= 0));
%~ assert (all (round (x) == x));
%~ assert (all (sum (x, 2) == n));

%~test
%~ n = (1:2)';
%~ p = [0.2, 0.5, 0.3; 0.1, 0.1, 0.8];
%~ x = mnrnd (n, p);
%~ assert (size (x), size (p));
%~ assert (all (x >= 0));
%~ assert (all (round (x) == x));
%~ assert (all (sum (x, 2) == n));
