function visualize_structure(comb_st, b_TopMass,L_soil,Dout, b_Clamp)
    figure();
    hold on
    
    % Draw the main structure
    for i=1:size(comb_st.dBotm,1)
            drawSymmetricRectangleFromBottom(comb_st.dBotm(i),comb_st.dTop(i), ... 
                (comb_st.zBot(i)-comb_st.zTop(i)), ... 
                comb_st.zTop(i)); hold on; hold on
    end
    xlim([-5 5])
    max_height =  max([comb_st.zTop;comb_st.zBot]);
    ylim([0 max_height*1.2])

    % Soil
    if b_Clamp
        % Nothing
        % ToDo visually show clamped state.
    else
        sss = linspace(0,L_soil,10);
        hold on
        for i=1:length(sss)
            add_spring(1, sss(i),Dout/2);
        end
        xlim([-5 10])
        ylim([-5 max(comb_st.zBot)*1.2])
    end
    
    % Top mass
    if b_TopMass
        [x_corners,y_corners] = draw_top_mass(max(comb_st.dBotm)*3, ... 
                                            max(comb_st.dBotm)*2, ...
                                            0, ... 
                                            max_height+max(comb_st.dBotm));
        patch(x_corners,y_corners, [0.5 0.5 0.5], 'EdgeColor', 'k', 'LineWidth', 1)
    end
    axis equal; 
end
