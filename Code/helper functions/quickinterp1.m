function Vout = quickinterp1(X,V,Xq)
%INTERP1 1-D interpolation (table lookup)
%
%   Vq = INTERP1(X,V,Xq) interpolates to find Vq, the values of the
%   underlying function V=F(X) at the query points Xq. 
%
%   X must be a vector. The length of X is equal to N.
%   If V is a vector, V must have length N, and Vq is the same size as Xq.
%   If V is an array of size [N,D1,D2,...,Dk], then the interpolation is
%   performed for each D1-by-D2-by-...-Dk value in V(i,:,:,...,:). If Xq
%   is a vector of length M, then Vq has size [M,D1,D2,...,Dk]. If Xq is 
%   an array of size [M1,M2,...,Mj], then Vq is of size
%   [M1,M2,...,Mj,D1,D2,...,Dk].
%
%   Vq = INTERP1(V,Xq) assumes X = 1:N, where N is LENGTH(V)
%   for vector V or SIZE(V,1) for array V.
%
%   Interpolation is the same operation as "table lookup".  Described in
%   "table lookup" terms, the "table" is [X,V] and INTERP1 "looks-up"
%   the elements of Xq in X, and, based upon their location, returns
%   values Vq interpolated within the elements of V.
%
%   Vq = INTERP1(X,V,Xq,METHOD) specifies the interpolation method.
%   The available methods are:
%
%     'linear'   - (default) linear interpolation
%     'nearest'  - nearest neighbor interpolation
%     'next'     - next neighbor interpolation
%     'previous' - previous neighbor interpolation
%     'spline'   - piecewise cubic spline interpolation (SPLINE)
%     'pchip'    - shape-preserving piecewise cubic interpolation
%     'cubic'    - same as 'pchip'
%     'v5cubic'  - the cubic interpolation from MATLAB 5, which does not
%                  extrapolate and uses 'spline' if X is not equally
%                  spaced.
%     'makima'   - modified Akima cubic interpolation
%
%   Vq = INTERP1(X,V,Xq,METHOD,'extrap') uses the interpolation algorithm
%   specified by METHOD to perform extrapolation for elements of Xq outside
%   the interval spanned by X.
%
%   Vq = INTERP1(X,V,Xq,METHOD,EXTRAPVAL) replaces the values outside of
%   the interval spanned by X with EXTRAPVAL.  NaN and 0 are often used for
%   EXTRAPVAL.  The default extrapolation behavior with four input
%   arguments is 'extrap' for 'spline', 'pchip' and 'makima', and
%   EXTRAPVAL = NaN (NaN+NaN*1i for complex values) for the other methods.
%
%   PP = INTERP1(X,V,METHOD,'pp') is not recommended. Use
%   griddedInterpolant instead.
%   PP = INTERP1(X,V,METHOD,'pp') uses the interpolation algorithm
%   specified by METHOD to generate the ppform (piecewise polynomial form)
%   of V. The method may be any of the above METHOD except for 'v5cubic'
%   and 'makima'. PP may then be evaluated via PPVAL. PPVAL(PP,Xq) is the
%   same as INTERP1(X,V,Xq,METHOD,'extrap').
%
%   For example, generate a coarse sine curve and interpolate over a
%   finer abscissa:
%       X = 0:10; V = sin(X); Xq = 0:.25:10;
%       Vq = interp1(X,V,Xq); plot(X,V,'o',Xq,Vq,':.')
%
%   For a multi-dimensional example, we construct a table of functional
%   values:
%       X = [1:10]'; V = [ X.^2, X.^3, X.^4 ];
%       Xq = [ 1.5, 1.75; 7.5, 7.75]; Vq = interp1(X,V,Xq);
%
%   creates 2-by-2 matrices of interpolated function values, one matrix for
%   each of the 3 functions. Vq will be of size 2-by-2-by-3.
%
%   Class support for inputs X, V, Xq, EXTRAPVAL:
%      float: double, single
%
%   See also INTERPFT, SPLINE, PCHIP, INTERP2, INTERP3, INTERPN, PPVAL,
%            griddedInterpolant, scatteredInterpolant.

%   Copyright 1984-2018 The MathWorks, Inc.


% Quick method: input needs to be sorted and all(X(1) < Xq < X(end))

% INTERP1(X,V,Xq)
F = griddedInterpolant(X,V,'linear');
Vout = F(Xq);

end % INTERP1

