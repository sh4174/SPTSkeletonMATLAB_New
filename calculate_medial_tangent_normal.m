function [ gradRList, tanList, normalList ] = calculate_medial_tangent_normal( skelIdx, radiiList )
    % Gradient of Radii on Skeleton
    gradRList = zeros( size( skelIdx ) );
    tanList = zeros( size( skelIdx ) );
    normalList = zeros( size( skelIdx ) );

    % Gauss Filter
    sigma = 0.5;
    sizeG = 5;
    hSize = floor( sizeG / 2 );
    Gauss_x = linspace( -hSize, hSize , sizeG );
    gaussF = exp( -Gauss_x .^ 2 / ( 2 * sigma ^2 ) );
    gaussF = gaussF ./ sum( gaussF(:) );

    for j = 1:size( skelIdx, 1 )
        if j == 1 || j == 2
            x_parts = skelIdx( j:j+2, 2 );
            y_parts = skelIdx( j:j+2, 1 );
            r_parts = radiiList( j:j+2 );

            x_parts_f = filter( gaussF( hSize+1:end ), 1, x_parts );
            y_parts_f = filter( gaussF( hSize+1:end ), 1, y_parts );
            r_parts_f = filter( gaussF( hSize+1:end ), 1, r_parts );

            radius_p = r_parts( 2 );
            radius_n = r_parts( 1 );

            x_p = x_parts_f( 2 );
            y_p = y_parts_f( 2 );

            x_n = x_parts_f( 1 );
            y_n = y_parts_f( 1 );
            dt = 1;
        elseif j == size( skelIdx, 1 ) || j == size( skelIdx, 1 ) - 1 
            x_parts = skelIdx( j-2:j, 2 );
            y_parts = skelIdx( j-2:j, 1 );
            r_parts = radiiList( j-2:j );

            x_parts_f = filter( gaussF( 1:hSize+1 ), 1, x_parts );
            y_parts_f = filter( gaussF( 1:hSize+1 ), 1, y_parts );
            r_parts_f = filter( gaussF( 1:hSize+1 ), 1, r_parts );

            radius_p = r_parts( hSize + 1 );
            radius_n = r_parts( hSize );

            x_p = x_parts_f( hSize + 1 );
            y_p = y_parts_f( hSize + 1 );

            x_n = x_parts_f( hSize );
            y_n = y_parts_f( hSize );

            dt = 1;

        else
            x_parts = skelIdx( j-2:j+2, 2 );
            y_parts = skelIdx( j-2:j+2, 1 );
            r_parts = radiiList( j-2:j+2 );

            x_parts_f = filter( gaussF, 1, x_parts );
            y_parts_f = filter( gaussF, 1, y_parts );
            r_parts_f = filter( gaussF, 1, r_parts );

            radius_p = r_parts( hSize + 2 );
            radius_n = r_parts( hSize );

            x_p = x_parts_f( hSize + 2 );
            y_p = y_parts_f( hSize + 2 );

            x_n = x_parts_f( hSize );
            y_n = y_parts_f( hSize );

            dt = 2;
        end

        % x', y', r'
        dxdt = ( x_p - x_n ) / dt;
        dydt = ( y_p - y_n ) / dt;
        drdt = ( radius_p - radius_n ) / dt;

        % Unit Tangent
        dxdy_mag = sqrt( dxdt^2 + dydt^2 );
        u_t = [ dxdt, dydt ] ./ dxdy_mag;
        tanList( j, : ) = u_t;
        normalList( j, : ) = [ -u_t( 2 ), u_t(1) ];

        grad_r = ( drdt / dxdy_mag ) .* u_t;
        gradRList( j, : ) = grad_r;
    end
end

