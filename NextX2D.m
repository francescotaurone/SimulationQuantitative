function [x1kp1, x2kp1] = NextX2D(x1k, x2k, u1, u2, Size)
    x1kp1 = NextX(x1k, u1, Size);
    x2kp1 = NextX(x2k, u2, Size);
end

