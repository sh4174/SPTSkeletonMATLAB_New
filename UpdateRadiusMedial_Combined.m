clear all;
addpath( 'ellipse' );

%% Read Data
for i = 0:5
    filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];
    fileFilledPath = [ 'ellipse_flow__t_' num2str( i ) '_filled.mat' ];

    filePtsPath = [ 'ellipse_flow__t_' num2str( i ) '_Pts.mat' ];
    fileSkelDPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo.mat' ];

    % Read Image
    bwImg = load( filePath );
    bw = bwImg.ellipseImg_p;


    bwFilled = load( fileFilledPath );
    bwFilled = bwFilled.ellipseF;


    % Read Surface points (Parametrized Curve)
    pts = load( filePtsPath );
    pts = pts.SourcePts_s;

    % Read Estimated Skeleton
    % skel = load( fileSkelPath );
    % skel = skel.skel;

    % Read Estimated Skeleton - Single Line
    skelD = load( fileSkelDPath );
    skelD = skelD.skelD;

    % % Display Object Surface and Skeleton
    % surf_skel = imfuse( bw, skel, 'blend', 'Scaling', 'joint' );
    % figure;
    % imshow( surf_skel );

    %% Calculate Initial Maximal Inscribed Circle - Closest Distance to Surface

    % Calculate Distance Map
    distSurf = bwdist( bw );

    % Skeleton Position - Line
    skelIdx = [];
    for r = 1:size( skelD, 1 )
    for c = 1:size( skelD, 2 )
        if skelD( r, c ) > 0
            skelIdx = [ skelIdx; r, c; ];
        end    
    end
    end

    %% Skeleton Parametrization - Single Curve
    skel_endPts = bwmorph( skelD, 'endpoints' );

    [ end_r, end_c ] = find( skel_endPts );

    if end_c( 1 ) < end_c( 2 )
        b_c = end_c( 1 );
        b_r = end_r( 1 );
        e_c = end_c( 2 );
        e_r = end_r( 2 );
    else
        b_c = end_c( 2 );
        b_r = end_r( 2 );
        e_c = end_c( 1 );
        e_r = end_r( 2 );
    end

    % Calculate Distance of Skeleton Points from Start Point
    distList = zeros( [ size( skelIdx, 1 ), 1 ] );

    for j = 1:size( skelIdx, 1 ) 
        skel_r = skelIdx( j, 1 );
        skel_c = skelIdx( j, 2 );

        dist = sqrt( ( b_r - skel_r ) * ( b_r - skel_r ) + ( b_c - skel_c ) * ( b_c - skel_c ) );
        distList( j ) = dist;
    end

    if issorted( distList ) == 0
        [ distList, sortIdx] = sort( distList );
        skelIdx = skelIdx( sortIdx, : );    
    end

    % Radius List of Skeleton to Object Surface
    radiiList = zeros( [ size( skelIdx, 1 ), 1 ] );

    for j = 1:size( skelIdx, 1 )
        skel_r = skelIdx( j, 1 );
        skel_c = skelIdx( j, 2 );

        radius = distSurf( skel_r, skel_c );
        radiiList( j ) = radius;
    end

    % Calculate Medial curve's radii gradient, curve tangent, curve normal
    [ gradRList, tanList, normalList ] = calculate_medial_tangent_normal( skelIdx, radiiList );
    disp( gradRList );

%     figure;
%     quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), gradRList( 1:10:end, 1 ), gradRList( 1:10:end, 2 ) );
%     set(gca,'Ydir','reverse')
%     title( 'Gradients' );

    % Calculate Corresponding Object Surface Normals
    [ U_p, U_n ] = calculate_surface_normal_from_medial( skelIdx, gradRList, normalList );

    %% Update Medial Curves 
    % Energy Function
    object_area = sum( bwFilled(:) );

    % Reconstruct Object Surface by Skeleton + Radii Function
    estimated_surface = zeros( size( bw ) );

    for j = 1:20:size( skelIdx, 1 )
        skel_r = skelIdx( j, 1 );
        skel_c = skelIdx( j, 2 );

        radius = distSurf( skel_r, skel_c );
        estimated_surface = MidpointCircle( estimated_surface, radius, skel_r, skel_c, 1 );    
    end

    overlapped_surface = bitand( logical( estimated_surface ), logical( bwFilled ) );
    overlapped_area = sum( overlapped_surface(:) );
    estimated_area = sum( estimated_surface(:) );
    energy_estimated = 1 / ( 2 * overlapped_area / ( object_area + estimated_area ) )


    %% Reconstruction 
    % Reconstruct Object Surface by Skeleton + Radii Function - Before Update
    skelIdx_t = skelIdx( 1:20:end, : );
    radiiList_t = radiiList( 1:20:end, : );

    recon = zeros( size( bw ) );

    for j = 1:size( skelIdx_t, 1 )
        skel_r = skelIdx_t( j, 1 );
        skel_c = skelIdx_t( j, 2 );

        radius = radiiList_t( j );
        recon = MidpointCircle( recon, radius, skel_r, skel_c, 1 );    
    end

    reconRGB = double( cat( 3, recon, recon, recon ) );

    for r = 1:size( skelD, 1 )
    for c = 1:size( skelD, 2 )
        if skelD( r, c ) > 0
            reconRGB( r, c, 1 ) = 0;
            reconRGB( r, c, 2 ) = 1;
            reconRGB( r, c, 3 ) = 0;
        end

        if bw( r, c ) > 0
            reconRGB( r, c, 1 ) = 1;
            reconRGB( r, c, 2 ) = 0;
            reconRGB( r, c, 3 ) = 0;
        end

    end
    end

    % imwrite( reconRGB, [ 'ellipse/SingleSkel/ellipse_flow__t_' num2str( i ) '_singleSkel.png' ] );

    % % Display Reconstructed Object and Surface Normal Plot
    % figure;
    % subplot( 1, 2, 1 );
    % imagesc( reconRGB );
    % hold on 
    % quiver( skelIdx_t( 1:end, 2 ), skelIdx_t( 1:end, 1 ), radiiList( 1:20:end ) .* U_p( 1:20:end, 1 ), radiiList( 1:20:end ) .* U_p( 1:20:end, 2 ) );
    % hold on
    % quiver( skelIdx_t( 1:end, 2 ), skelIdx_t( 1:end, 1 ), radiiList( 1:20:end ) .* U_n( 1:20:end, 1 ), radiiList( 1:20:end ) .* U_n( 1:20:end, 2 ) );
    % hold off
    % title( [ 'reconstructed : Img ' num2str( i ) ] );

    % Update medial curve positions by gradient descent with Gateaux derivative
    [ gradList, grad2List ] = calculate_gradients( skelIdx_t );

    gradNorm = gradList .* gradList;
    gradNorm = sum( gradNorm, 2 );
    gradNorm_prev = mean( sqrt( gradNorm ) );

    grad2Norm = grad2List .* grad2List;
    grad2Norm_prev = mean( sqrt( sum( grad2Norm, 2 ) ) );

    gd_m_step = 1;
    gd_r_step = 1;


    gd_alpha = 2000;
    gd_r_alpha = 2000;

    smooth_eps = 0.5;
    grad2_eps = 1.0;

    energy_prev = energy_estimated + smooth_eps * ( gradNorm_prev + grad2_eps * grad2Norm_prev );
    energy_0 = energy_estimated + smooth_eps * ( gradNorm_prev + grad2_eps * grad2Norm_prev );

    energy_prev = energy_prev / energy_0;

    disp( 'Initial Energy' );
    disp( energy_prev );

    % Interpolation Ratio
    interp_ratio = 1;
    
    
    for t = 1:10
        % Update Medial Curve
        for j = 1:size( skelIdx_t, 1 )
        for k = 1:2
            skelIdx_temp = skelIdx_t;
            skelIdx_temp( j, k ) = skelIdx_t( j, k ) + gd_m_step;        
            estimated_surface = zeros( size( bw ) );

            % Interpolate
            skelIdx_temp_int_r = round( interp( skelIdx_temp( :, 1 ), interp_ratio ) );
            skelIdx_temp_int_c = round( interp( skelIdx_temp( :, 2 ), interp_ratio ) );
            
            radiiList_t_int = round( interp( radiiList_t, interp_ratio ) );

            skelIdx_temp_int = [ skelIdx_temp_int_r( 1:(end-interp_ratio+1), 1 ), skelIdx_temp_int_c( 1:(end-interp_ratio+1), 1 ) ];           
            radiiList_t_int = radiiList_t_int( 1:(end-interp_ratio+1), 1 );
                       
            
            for j2 = 1:size( skelIdx_temp_int, 1 )
                skel_r = skelIdx_temp_int( j2, 1 );
                skel_c = skelIdx_temp_int( j2, 2 );

                radius = radiiList_t_int( j2 );
                try
                    estimated_surface = MidpointCircle( estimated_surface, radius, skel_r, skel_c, 1 );    
                catch
                    dkdkslslsl = 0;
                end
            end

            [ gradList, grad2List ] = calculate_gradients( skelIdx_temp );

            gradNorm_temp = gradList .* gradList;
            gradNorm_temp = mean( sqrt( sum( gradNorm_temp, 2 ) ) );                

            grad2Norm_temp = grad2List .* grad2List;
            grad2Norm_temp = mean( sqrt( sum( grad2Norm_temp, 2 ) ) );                

            overlapped_surface = bitand( logical( estimated_surface ), logical( bwFilled ) );
            overlapped_area = sum( overlapped_surface(:) );
            estimated_area = sum( estimated_surface(:) );

            energy_i_t = 1 / ( 2 * overlapped_area / ( object_area + estimated_area ) );
            energy_s_t = gradNorm_temp + grad2_eps * grad2Norm_temp;

            energy_estimated = ( energy_i_t + smooth_eps * energy_s_t ) / energy_0;

            gateaux_dev = ( energy_estimated - energy_prev ) / gd_m_step;

            skelIdx_t( j, k ) = skelIdx_t( j, k ) - gd_alpha * gateaux_dev;      
            energy_prev = energy_estimated;
        end
        end

        disp( 'Medial Curve Update' );
        disp( [ 'Estimated Energy at Iter - ' num2str( t ) ' : ' num2str( energy_estimated ) ] );
        disp( [ 'Energy Overlapped - ' num2str( energy_i_t ) ] );
        disp( [ 'Dice Coeff - ' num2str( 1.0 / energy_i_t ) ] );
        disp( [ 'Energy Curve Smoothness - ' num2str( energy_s_t ) ] );
        disp( [ 'Gateaux Derivative : ' num2str( gateaux_dev ) ] );

        % Update Radial Function
        % Medial Curve Regularization does not change but included for energy
        % normalization
        [ gradList_r, grad2List_r ] = calculate_gradients( skelIdx_t );

        gradNorm_r = gradList_r .* gradList_r;
        gradNorm_r = mean( sqrt( sum( gradNorm_r, 2 ) ) );                

        grad2Norm_r = grad2List_r .* grad2List_r;
        grad2Norm_r = mean( sqrt( sum( grad2Norm_r, 2 ) ) );                
        energy_s_t_r = gradNorm_r + grad2_eps * grad2Norm_r;

        for j = 1:size( radiiList_t, 1 )
            radiiList_temp = radiiList_t;
            radiiList_temp( j ) = radiiList_t( j ) + gd_r_step;
            estimated_surface = zeros( size( bw ) );

            skelIdx_t_int_r = round( interp( skelIdx_t( :, 1 ), interp_ratio ) );
            skelIdx_t_int_c = round( interp( skelIdx_t( :, 2 ), interp_ratio ) );
            
            radiiList_temp_int = round( interp( radiiList_temp, interp_ratio ) );
      
            skelIdx_t_int = [ skelIdx_t_int_r( 1:(end-interp_ratio+1), 1 ), skelIdx_t_int_c( 1:(end-interp_ratio+1), 1 ) ];           
            radiiList_temp_int = radiiList_temp_int( 1:(end-interp_ratio+1), 1 );
            
            for j2 = 1:size( skelIdx_t_int, 1 )
                skel_r = skelIdx_t_int( j2, 1 );
                skel_c = skelIdx_t_int( j2, 2 );

                radius = radiiList_temp_int( j2 );
                try
                    estimated_surface = MidpointCircle( estimated_surface, radius, skel_r, skel_c, 1 );    
                catch
                    dkdkslslsl = 0;
                end
            end

            overlapped_surface = bitand( logical( estimated_surface ), logical( bwFilled ) );
            overlapped_area = sum( overlapped_surface(:) );
            estimated_area = sum( estimated_surface(:) );

            energy_i_t = 1 / ( 2 * overlapped_area / ( object_area + estimated_area ) );

            energy_estimated = ( energy_i_t + smooth_eps * energy_s_t_r ) / energy_0;

            gateaux_dev = ( energy_estimated - energy_prev ) / gd_r_step;

            radiiList_t( j ) = radiiList_t( j ) - gd_r_alpha * gateaux_dev;      
            energy_prev = energy_estimated;        
        end    


        if abs( gateaux_dev ) < 0.0000001
            disp( 'Optimization Terminated : dE/dm < threshold' );
            break;        
        end

        disp( 'Radius Update' );
        disp( [ 'Estimated Energy at Iter - ' num2str( t ) ' : ' num2str( energy_estimated ) ] );
        disp( [ 'Energy Overlapped - ' num2str( energy_i_t ) ] );
        disp( [ 'Dice Coeff - ' num2str( 1.0 / energy_i_t ) ] );
        disp( [ 'Energy Curve Smoothness - ' num2str( energy_s_t ) ] );
        disp( [ 'Gateaux Derivative : ' num2str( gateaux_dev ) ] );
    end


%     %% Display
%     % Display Surface Normal Estimated from Medial Axis
%     figure;
%     quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), U_p( 1:10:end, 1 ), U_p( 1:10:end, 2 ) );
%     hold on
%     quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), U_n( 1:10:end, 1 ), U_n( 1:10:end, 2 ) );
%     set(gca,'Ydir','reverse')
%     hold off
%     title( 'Surface Normals' );
%     
%     figure;
%     quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), normalList( 1:10:end, 1 ), normalList( 1:10:end, 2 ) );
%     hold on
%     quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), -normalList( 1:10:end, 1 ), -normalList( 1:10:end, 2 ) );
%     set(gca,'Ydir','reverse')
%     hold off
%     title( 'Normals' );

    %% Reconstruction 
    % Reconstruct Object Surface by Skeleton + Radii Function
    visInterp_ratio = 50;
    
    recon = zeros( size( bw ) );
    skelIdx_t_int_r = round( interp( skelIdx_t( :, 1 ), visInterp_ratio ) );
    skelIdx_t_int_c = round( interp( skelIdx_t( :, 2 ), visInterp_ratio ) );
    
    radiiList_t_int = round( interp( radiiList_t, visInterp_ratio ) );
    
    skelIdx_t_int = [ skelIdx_t_int_r( 1:(end-visInterp_ratio+1), 1 ), skelIdx_t_int_c( 1:(end-visInterp_ratio+1), 1 ) ];           
    radiiList_t_int = radiiList_t_int( 1:(end-visInterp_ratio+1), 1 );

    
    for j = 1:size( skelIdx_t_int, 1 )
        skel_r = skelIdx_t_int( j, 1 );
        skel_c = skelIdx_t_int( j, 2 );

        radius = radiiList_t_int( j );
        recon = MidpointCircle( recon, radius, skel_r, skel_c, 1 );    
    end

    reconRGB = double( cat( 3, recon, recon, recon ) );

    for j = 1:size( skelIdx_t_int, 1 )
        r = skelIdx_t_int( j, 1 );
        c = skelIdx_t_int( j, 2 );
        
        for dr = -2:2
            for dc = -2:2
                reconRGB( r + dr, c + dc, 1 ) = 0;
                reconRGB( r + dr, c + dc, 2 ) = 1;
                reconRGB( r + dr, c + dc, 3 ) = 0;
            end
        end
    end
    
    for r = 1:size( skelD, 1 )
    for c = 1:size( skelD, 2 )
        if bw( r, c ) > 0
            for dr = -2:2
                for dc = -2:2
                    reconRGB( r + dr, c + dc, 1 ) = 1;
                    reconRGB( r + dr, c + dc, 2 ) = 0;
                    reconRGB( r + dr, c + dc, 3 ) = 0;
                end
            end
        end
    end
    end

    imwrite( reconRGB, [ 'ellipse/SingleSkel3/ellipse_flow__t_' num2str( i ) '_singleSkel.png' ] );

%     % Calculate Medial curve's radii gradient, curve tangent, curve normal
%     [ gradRList_updated, tanList_updated, normalList_updated ] = calculate_medial_tangent_normal( skelIdx_t, radiiList_t );
% 
%     % Calculate Corresponding Object Surface Normals
%     [ U_p_updated, U_n_updated ] = calculate_surface_normal_from_medial( skelIdx_t, gradRList_updated, normalList_updated );


%     Display Reconstructed Object and Surface Normal Plot
%     subplot( 1, 2, 2 );
%     imagesc( reconRGB );
%     hold on 
%     quiver( skelIdx_t( 1:end, 2 ), skelIdx_t( 1:end, 1 ), radiiList_t( 1:end ) .* U_p_updated( 1:end, 1 ), radiiList_t( 1:end ) .* U_p_updated( 1:end, 2 ) );
%     hold on
%     quiver( skelIdx_t( 1:end, 2 ), skelIdx_t( 1:end, 1 ), radiiList_t( 1:end ) .* U_n_updated( 1:end, 1 ), radiiList_t( 1:end ) .* U_n_updated( 1:end, 2 ) );
%     hold on 
%     plot( skelIdx_t( :, 2 ), skelIdx_t( :, 1 ), 'r', 'LineWidth', 2 );
%     hold off
%     title( [ 'reconstructed : Img ' num2str( i ) ] );
end