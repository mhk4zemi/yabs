function [Ix, Iy] = calculate_area_moment_of_inertia(Dout,Din)
    Ix  = pi .* (Dout.^4.0 - Din.^4.0) ./ 64.0;
    Iy = Ix;
end
