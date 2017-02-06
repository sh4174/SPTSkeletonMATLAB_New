function [ tangents, normals ] = calc_curve_tangent_normal( pts )
    dx_pts = [];    
    dy_pts = [];
    
    for j = 1:size( pts, 1 )
        
        if pts( j, 1 ) == 0 && pts( j, 2 ) == -30
            disp( 'idx' );
            disp( j );
        end
        
        if j == 1
            p_x_n1 = pts( size( pts, 1 ), 1 );
            p_y_n1 = pts( size( pts, 1 ), 2 );
            
            p_x_p1 = pts( 2, 1 );
            p_y_p1 = pts( 2, 2 );
        elseif j == size( pts, 1 )
            p_x_n1 = pts( j - 1, 1 );
            p_y_n1 = pts( j - 1, 2 );
            
            p_x_p1 = pts( 1, 1 );
            p_y_p1 = pts( 1, 2 );    
        else
            p_x_n1 = pts( j - 1, 1 );
            p_y_n1 = pts( j - 1, 2 );
            
            p_x_p1 = pts( j + 1, 1 );
            p_y_p1 = pts( j + 1, 2 );                   
        end
        
        dx = p_x_p1 - p_x_n1;
        dy = p_y_p1 - p_y_n1;
        
        dx_n = dx / sqrt( dx * dx + dy * dy );
        dy_n = dy / sqrt( dx * dx + dy * dy );
        
        dx_pts = [ dx_pts; dx_n ];
        dy_pts = [ dy_pts; dy_n ];
    end
    
    tangents = [ dx_pts, dy_pts ];
    normals = [ -dy_pts, dx_pts ];    
end

