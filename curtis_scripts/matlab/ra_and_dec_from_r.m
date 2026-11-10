% ALGORITHM 4.1: OBTAIN THE RIGHT ASCENSION AND DECLINATION
% FROM THE POSITION VECTOR

function [ra dec] = ra_and_dec_from_r(r)
%{
    This function calculates the right ascension and the
    declination from the geocentric equatorial position vector.

    r       - position vector
    l, m, n - direction cosines of r
    ra      - right ascension (degrees)
    dec     - declination (degrees)
%}
% ------------------------------------------------
if numel(r) ~= 3 || any(~isfinite(r(:))) || norm(r) == 0
    error('Position must be a finite, nonzero three-vector.')
end
dec = atan2d(r(3),hypot(r(1),r(2)));
ra = mod(atan2d(r(2),r(1)),360);
end