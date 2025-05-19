function drawSymmetricRectangleFromBottom(width2,width, height, bottomY)
    % Function to draw a rectangle symmetric around x = 0
    % Inputs:
    % - width: width of the rectangle
    % - height: height of the rectangle
    % - bottomY: y-coordinate of the bottom side of the rectangle
    
    % Calculate rectangle vertices
    x = [-width/2, width/2, width2/2, -width2/2]; % Symmetric around x = 0
    y = [bottomY, bottomY, bottomY + height, bottomY + height];
    
    % Plot the rectangle
    fill(x, y, 'cyan'); % Cyan-filled rectangle
    hold on;
    axis equal;
    
    % Add grid and labels
    grid on;
    xlabel('x-axis');
    ylabel('y-axis');
    title('Rectangle Symmetric Around x = 0');
    
    % Add symmetry line for reference
%     xline(0, 'r--', 'LineWidth', 1.5); % Dashed red line at x = 0
    
%     hold off;
end
