function rho = SplineInterp(varargin)
switch numel(varargin)
    case 4
        rho = SplineInterp3D(varargin{:});
    case 3
        rho = SplineInterp2D(varargin{:});
    otherwise
        rho = NaN;
        error('invalid number of inputs');
end