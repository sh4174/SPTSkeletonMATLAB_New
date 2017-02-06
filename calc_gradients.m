function [ gradList, grad2List ] = calc_gradients( skelIdx )
    % Gradient of Skeleton
    gradList = zeros( size( skelIdx ) );
    grad2List = zeros( size( skelIdx ) );

    % Gauss Filter
    sigma = 1.0;
    sizeG = 3;
    hSize = floor( sizeG / 2 );
    Gauss_x = linspace( -hSize, hSize , sizeG );
    gaussF = exp( -Gauss_x .^ 2 / ( 2 * sigma ^2 ) );
    gaussF = gaussF ./ sum( gaussF(:) );

    pad_v = 6;
    padmat = ones( pad_v, 1 );
    skelIdx_r_pad = [ padmat * skelIdx( 1, 1 ); skelIdx( :, 1 ); padmat * skelIdx( end, 1 ) ];
    skelIdx_c_pad = [ padmat * skelIdx( 1, 2 ); skelIdx( :, 2 ); padmat * skelIdx( end, 2 ) ];
    
    skelIdx_r_f = filter( gaussF, 1, skelIdx_r_pad );
    skelIdx_c_f = filter( gaussF, 1, skelIdx_c_pad );
    
    skelIdx_r = skelIdx_r_f( pad_v:size( skelIdx_r_f, 1) - pad_v + 1, :);
    skelIdx_c = skelIdx_c_f( pad_v:size( skelIdx_c_f, 1) - pad_v + 1, :);
    
    for j = 1:size( skelIdx, 1 )
        if j == 1
            dt = 1;
            x_p = skelIdx_c( 1, 1 );
            x_n = skelIdx_c( 2, 1 );            
            y_p = skelIdx_r( 1, 1 );
            y_n = skelIdx_r( 2, 1 );            
        elseif j == size( skelIdx, 1 )
            dt = 1;
            x_p = skelIdx_c( end-1, 1 );
            x_n = skelIdx_c( end, 1 );            
            y_p = skelIdx_r( end-1, 1 );
            y_n = skelIdx_r( end, 1 );            
        else
            dt = 2;
            x_p = skelIdx_c( j, 1 );
            x_n = skelIdx_c( j + 2, 1 );
            y_p = skelIdx_r( j, 1 );
            y_n = skelIdx_r( j + 2, 1 );
        end
        
        % x', y', r'
        dxdt = ( x_p - x_n ) / dt;
        dydt = ( y_p - y_n ) / dt;

        % Unit Tangent
        dxdy_mag = sqrt( dxdt^2 + dydt^2 );
        u_t = [ dxdt, dydt ];
        gradList( j, : ) = u_t;
    end
    
     
    for i = 1:size( gradList, 1 )
        if i == 1
            x_u_p = gradList( 1, 1 );
            y_u_p = gradList( 1, 2 );
            
            x_u_n = gradList( 2, 1 );
            y_u_n = gradList( 2, 2 );
            dt = 1;
        elseif i == size( gradList, 1 )
            x_u_p = gradList( size( gradList, 1 ) - 1, 1 );
            y_u_p = gradList( size( gradList, 1 ) - 1, 2 );
            
            x_u_n = gradList( size( gradList, 1 ), 1 );
            y_u_n = gradList( size( gradList, 1 ), 2 );
            dt = 1;
        else
            x_u_p = gradList( i - 1, 1 );
            y_u_p = gradList( i - 1, 2 );
            
            x_u_n = gradList( i + 1, 1 );
            y_u_n = gradList( i + 1, 2 );
            dt = 2;
        end
        
        ddxdt = ( x_u_p - x_u_n ) / dt;
        ddydt = ( y_u_p - y_u_n ) / dt;
        ddxdy_mag = sqrt( ddxdt^2 + ddydt^2 );
        du_t = [ ddxdt, ddydt ];
        grad2List( j, : ) = du_t;
    end
end
