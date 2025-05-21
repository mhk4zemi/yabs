function results = make_timoshenko_beam(comb_st, spacing, rho, E, nu)

    Zint= linspace(min(comb_st.zTop), ... 
                   max(comb_st.zBot),...
                   floor((max(comb_st.zBot)-min(comb_st.zTop))/spacing));
    initialization_array = zeros(1,length(Zint));

    Zarr = initialization_array;
    Marr = initialization_array;
    Aarr = initialization_array;
    Ixarr = initialization_array;
    Iyarr = initialization_array;
    rxarr = initialization_array;
    ryarr = initialization_array;
    Garr = initialization_array;
    Jarr = initialization_array;
    HcogArr = initialization_array;
    MassArr = initialization_array;
    flag = false;

    for i =1:length(Zint)
        hTop = comb_st.zTop;
        hBot = comb_st.zBot;
        idx1 = find(Zint(i)>=hTop,1,'last');
        data_can.zTop(1) = hTop(idx1);
        data_can.zBot(1) = hBot(idx1);
        data_can.dBotm(1) = comb_st.dBotm(idx1);
        data_can.dTop(1) =  comb_st.dTop(idx1);
        data_can.Thick(1) =  comb_st.Thick(idx1);

        [A,Din,Dout] = area_calculation(data_can,1,Zint(i));
         mL = A *rho;
        [Ix, Iy] = calculate_area_moment_of_inertia(Dout,Din);
         G = E / (2 * (1 + nu));
         rix = (Ix ./ A) .^ 0.5;
         riy = (Iy ./ A) .^ 0.5;  

        Zarr(i) = Zint(i);
        Marr(i) = mL;
        Aarr(i)  = A;
        Ixarr(i) = Ix;
        Iyarr(i) = Iy;
        rxarr(i) = rix;
        ryarr(i) = riy;
    end

    zero_array  = repmat(0, length(Zarr), 1);

    results = table(Zarr', ...
                    Marr', ...
                    zero_array,zero_array, ...
                    rxarr',ryarr', zero_array,zero_array, ...
                    repmat(E, length(Zarr), 1), ... 
                    repmat(G, length(Zarr), 1), ... 
                    Ixarr', Iyarr', ...
                    2*Ixarr', ...
                    repmat(2*(1+nu)/(4+3*nu), length(Zarr), 1), repmat(2*(1+nu)/(4+3*nu), length(Zarr), 1), ...
                    Aarr', ...
                    zero_array, zero_array , zero_array , ...
                    'VariableNames', {'H', 'm', 'x_cg', 'y_cg', 'ri_x', 'ri_y', 'x_sh', 'y_sh', ...
                                      'E', 'G', 'I_x', 'I_y', 'J', 'Kx', 'Ky', 'A', 'pitch', 'xe', 'ye'});
end
