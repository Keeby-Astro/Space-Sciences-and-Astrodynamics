function r = quat_rotate(q,v)
%{
    quat_rotate rotates a vector by a unit quaternion.
    r = quat_rotate(q,v) calculates the rotated vector r for a
        quaternion q and a vector v.
    q   is a 1-by-4 matrix whose norm must be 1. q(1) is the scalar part
        of the quaternion.
    v   is a 1-by-3 matrix.
    r   is a 1-by-3 matrix.

    The 3-vector v is made into a pure quaternion 4-vector V = [0 v]. r is
    produced by the quaternion product R = q*V*qinv. r = [R(2) R(3) R(4)].
    
    MATLAB M-functions used: quatmultiply, quatinv.
%}
% ---------------------------------------------------
qinv = quatinv(q);
r    = quatmultiply(quatmultiply(q,[0 v]),qinv);
r    = r(2:4);
end %quat_rotate