function [x_corners,y_corners] =draw_top_mass(width,height,center_x,center_y);
    % Calculate the vertices of the box
    half_width = width / 2;
    half_height = height / 2;

    % Define the corner points of the box
    x_corners = [center_x - half_width, center_x + half_width, center_x + half_width, center_x - half_width, center_x - half_width];
    y_corners = [center_y - half_height, center_y - half_height, center_y + half_height, center_y + half_height, center_y - half_height];

end
