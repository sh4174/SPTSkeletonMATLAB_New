clear all;
addpath( 'ellipse' );

%% Read Data
sampling_interval = 20;

for i = 1:1
    filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];
    fileFilledPath = [ 'ellipse_flow__t_' num2str( i ) '_filled.mat' ];
    
    if i == 0
        fileSkelDPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo.mat' ];
    else
        fileSkelDPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo.mat' ];
        fileSkelInitPath = [ 'ellipse_flow__t_' num2str( i - 1 ) '_skel_0topo_Updated.mat' ];
        fileRadiusInitPath = [ 'ellipse_flow__t_' num2str( i - 1 ) '_radius_Updated.mat' ];
    end

    % Read Image
    % Target Image (Contour)
    bwImg = load( filePath );
    bwImg = bwImg.ellipseImg_p;
    % Target Image (Filled)
    bwFilled = load( fileFilledPath );
    bwFilled = bwFilled.ellipseF;

    % Read Reference Skeleton - Not Used for Processing : Skeleton from
    % Morphological Op
    skelD = load( fileSkelDPath );
    skelD = skelD.skelD;

    % Read Initial Skeleton from i - 1 th data and radii 
    skelIdx = load( fileSkelInitPath );
    skelIdx = skelIdx.skelIdx_t;
    radiiList = load( fileRadiusInitPath );
    radiiList = radiiList.radiiList_t;
    
    skelIdx_t = round( skelIdx );
    radiiList_t = radiiList;

    % Target Area
    object_area = sum( bwFilled(:) );
    % Reconstruct Object Surface by Skeleton + Radii Function
    estimated_surface = zeros( size( bwImg ) );

    for j = 1:size( skelIdx, 1 )
        skel_r = round( skelIdx( j, 1 ) );
        skel_c = round( skelIdx( j, 2 ) );

        radius = radiiList( j );
        estimated_surface = MidpointCircle( estimated_surface, radius, skel_r, skel_c, 1 );    
    end       
    overlapped_surface = bitand( logical( estimated_surface ), logical( bwFilled ) );
    overlapped_area = sum( overlapped_surface(:) );
    estimated_area = sum( estimated_surface(:) );

    % Image Matching Energy - Inverse Dice Coefficient
    energy_estimated = 1 / ( 2 * overlapped_area / ( object_area + estimated_area ) )

    % Medial Curve Smoothness/Elasticity Energy (Smoothness Regularization)
    [ gradListX, gradList ] = gradient( skelIdx_t );    
    [ grad2ListX, grad2List ] = gradient( gradList );

    gradNorm = gradList .* gradList;
    gradNorm = sum( gradNorm, 2 );
    gradNorm_prev = mean( sqrt( gradNorm ) );

    grad2Norm = grad2List .* grad2List;
    grad2Norm_prev = mean( sqrt( sum( grad2Norm, 2 ) ) );

    gd_m_step = 0.5;
    gd_r_step = 0.5;

    gd_alpha = 1000;
    gd_r_alpha = 1000;

    smooth_eps = 0.00005;
    grad2_eps = 50;

    energy_smooth_prev = gradNorm_prev + grad2_eps * grad2Norm_prev;
    energy_0 = energy_estimated + smooth_eps * ( energy_smooth_prev );
    energy_prev = energy_0;

    disp( 'Initial Energy' );
    disp( energy_prev );

    % Interpolation Ratio
    interp_ratio = 1;
        
    iter_t = 5;
    iter_m_t = 50;
    iter_r_t = 50;
    
    EnergyArr_m = zeros( iter_t, iter_m_t );
    EnergyArr_r = zeros( iter_t, iter_r_t );
    
    DiceArr_m = zeros( iter_t, iter_m_t );
    DiceArr_r = zeros( iter_t, iter_r_t );
     
    visInterp_ratio = 50;        
    [ recon_org, reconRGB_Org ] = ReconOutputSkel( skelIdx_t, radiiList_t, bwImg, skelD, visInterp_ratio );
    imwrite( reconRGB_Org, [ 'ellipse/SingleSkel_Result/ellipse_flow__t_' num2str( i ) '_singleSkel_Org.png' ] );
        
    for t = 1:iter_t
        % Update Medial Curve
        for t_m = 1:iter_m_t
            gateaux_dev_mc_max = 0;
            gateux_dev_mc_sum = 0;
            for j = 1:size( skelIdx_t, 1 )
            for k = 1:2
                skelIdx_temp = skelIdx_t;
                skelIdx_temp( j, k ) = skelIdx_t( j, k ) + gd_m_step;        
                estimated_surface = zeros( size( bwImg ) );

                % Interpolate
                skelIdx_temp_int_r = interp( skelIdx_temp( :, 1 ), interp_ratio );
                skelIdx_temp_int_c = interp( skelIdx_temp( :, 2 ), interp_ratio );

                radiiList_t_int = interp( radiiList_t, interp_ratio );

                skelIdx_temp_int = [ skelIdx_temp_int_r( 1:(end-interp_ratio+1), 1 ), skelIdx_temp_int_c( 1:(end-interp_ratio+1), 1 ) ];           
                radiiList_t_int = radiiList_t_int( 1:(end-interp_ratio+1), 1 );

                for j2 = 1:size( skelIdx_temp_int, 1 )
                    skel_r = round( skelIdx_temp_int( j2, 1 ) );
                    skel_c = round( skelIdx_temp_int( j2, 2 ) );

                    radius = round( radiiList_t_int( j2 ) );
                    try
                        estimated_surface = MidpointCircle( estimated_surface, radius, skel_r, skel_c, 1 );    
                    catch
                        dkdkslslsl = 0;
                    end
                end
                
                [ gradListX, gradList ] = gradient( skelIdx_temp );    
                [ grad2ListX, grad2List ] = gradient( gradList );
        
                gradNorm_temp = gradList .* gradList;
                gradNorm_temp = mean( sqrt( sum( gradNorm_temp, 2 ) ) );                

                grad2Norm_temp = grad2List .* grad2List;
                grad2Norm_temp = mean( sqrt( sum( grad2Norm_temp, 2 ) ) );                

                overlapped_surface = bitand( logical( estimated_surface ), logical( bwFilled ) );
                overlapped_area = sum( overlapped_surface(:) );
                estimated_area = sum( estimated_surface(:) );

                energy_i_t = 1 / ( 2 * overlapped_area / ( object_area + estimated_area ) );
                energy_s_t = gradNorm_temp + grad2_eps * grad2Norm_temp;

                energy_estimated = ( energy_i_t + smooth_eps * energy_s_t ); % / energy_0;

                gateaux_dev = ( energy_estimated - energy_prev ) / gd_m_step;

                gateux_dev_mc_sum = gateux_dev_mc_sum + abs( gateaux_dev );

                if abs( gateaux_dev ) > gateaux_dev_mc_max
                    gateaux_dev_mc_max = abs( gateaux_dev );
                end
                skelIdx_t( j, k ) = skelIdx_t( j, k ) - ( gd_alpha * gateaux_dev );      
                energy_prev = energy_estimated;
            end
            end
            
            EnergyArr_m( t, t_m ) = energy_estimated;            
            DiceArr_m( t, t_m ) =  1.0 / energy_i_t ;
            
            disp( 'Medial Curve Update : i ' );
            disp( [ 'Estimated Energy at Iter - ' num2str( t_m ) ' : ' num2str( energy_estimated ) ] );
            disp( [ 'Energy Overlapped - ' num2str( energy_i_t ) ] );
            disp( [ 'Dice Coeff - ' num2str( 1.0 / energy_i_t ) ] );
            disp( [ 'Energy Curve Smoothness - ' num2str( energy_s_t ) ] );
            disp( [ 'Mean Gateaux Derivative : ' num2str( gateux_dev_mc_sum / ( 2 * size( skelIdx_t, 1 ) ) ) ] );
            disp( [ 'Max Gateaux Derivative : ' num2str( gateaux_dev_mc_max ) ] );
            
            if  gateaux_dev_mc_max < 0.000005
                disp( 'Optimization Terminated : dE/dm < threshold' );
                break;        
            end
        end
        
        skelIdx_t( :, 1 ) = smooth( skelIdx_t( :, 1 ), 3 );
        skelIdx_t( :, 2 ) = smooth( skelIdx_t( :, 2 ), 3 );
        
        dev_m = gateux_dev_mc_sum / ( 2 * size( skelIdx_t, 1 ) );
        disp( '====================================================================' );
        disp( '====================================================================' );
        disp( 'Medial Curve Update' );
        disp( [ 'Estimated Energy at Iter - ' num2str( t ) ' : ' num2str( energy_estimated ) ] );
        disp( [ 'Energy Overlapped - ' num2str( energy_i_t ) ] );
        disp( [ 'Dice Coeff - ' num2str( 1.0 / energy_i_t ) ] );
        disp( [ 'Energy Curve Smoothness - ' num2str( energy_s_t ) ] );
        disp( [ 'Mean Gateaux Derivative : ' num2str( gateux_dev_mc_sum / ( 2 * size( skelIdx_t, 1 ) ) ) ] );
        disp( [ 'Max Gateaux Derivative : ' num2str( gateaux_dev_mc_max ) ] );
        disp( '====================================================================' );
        disp( '====================================================================' );
        
        visInterp_ratio = 50;        
        [ recon_mc, reconRGB_mc ] = ReconOutputSkel( skelIdx_t, radiiList_t, bwImg, skelD, visInterp_ratio );
        imwrite( reconRGB_mc, [ 'ellipse/SingleSkel_Result/ellipse_flow__t_' num2str( i ) '_singleSkel_MC' num2str( t ) '.png' ] );
        
        % Update Radial Function
        % Medial Curve Regularization does not change but included for energy
        % normalization
        [ gradListX, gradList_r ] = gradient( skelIdx_t );    
        [ grad2ListX, grad2List_r ] = gradient( gradList_r );
        
        gradNorm_r = gradList_r .* gradList_r;
        gradNorm_r = mean( sqrt( sum( gradNorm_r, 2 ) ) );                

        grad2Norm_r = grad2List_r .* grad2List_r;
        grad2Norm_r = mean( sqrt( sum( grad2Norm_r, 2 ) ) );                
        energy_s_t_r = gradNorm_r + grad2_eps * grad2Norm_r;

        energy_s_t = energy_s_t_r;
        
        for t_r = 1:iter_r_t
            gateux_dev_r_sum = 0;
            gateaux_dev_r_max = 0;

            for j = 1:size( radiiList_t, 1 )
                radiiList_temp = radiiList_t;
                radiiList_temp( j ) = radiiList_t( j ) + gd_r_step;
                estimated_surface = zeros( size( bwImg ) );

                skelIdx_t_int_r = interp( skelIdx_t( :, 1 ), interp_ratio );
                skelIdx_t_int_c = interp( skelIdx_t( :, 2 ), interp_ratio );

                radiiList_temp_int = interp( radiiList_temp, interp_ratio );

                skelIdx_t_int = [ skelIdx_t_int_r( 1:(end-interp_ratio+1), 1 ), skelIdx_t_int_c( 1:(end-interp_ratio+1), 1 ) ];           
                radiiList_temp_int = radiiList_temp_int( 1:(end-interp_ratio+1), 1 );

                for j2 = 1:size( skelIdx_t_int, 1 )
                    skel_r = round( skelIdx_t_int( j2, 1 ) );
                    skel_c = round( skelIdx_t_int( j2, 2 ) );

                    radius = round( radiiList_temp_int( j2 ) );
                    estimated_surface = MidpointCircle( estimated_surface, radius, skel_r, skel_c, 1 );    
                end

                overlapped_surface = bitand( logical( estimated_surface ), logical( bwFilled ) );
                overlapped_area = sum( overlapped_surface(:) );
                estimated_area = sum( estimated_surface(:) );

                energy_i_t = 1 / ( 2 * overlapped_area / ( object_area + estimated_area ) );

                energy_estimated = ( energy_i_t + smooth_eps * energy_s_t_r );

                gateaux_dev = ( energy_estimated - energy_prev ) / gd_r_step;

                gateux_dev_r_sum = gateux_dev_r_sum + abs( gateaux_dev );

                if abs( gateaux_dev ) > gateaux_dev_r_max
                    gateaux_dev_r_max = abs( gateaux_dev );
                end

                radiiList_t( j ) = radiiList_t( j ) - ( gd_r_alpha * gateaux_dev );      
                energy_prev = energy_estimated;                  
            end    

            dev_r = gateux_dev_r_sum / size( radiiList_t, 1 );

            disp( 'Radius Update: i' );
            disp( [ 'Estimated Energy at Iter - ' num2str( t_r ) ' : ' num2str( energy_estimated ) ] );
            disp( [ 'Energy Overlapped - ' num2str( energy_i_t ) ] );
            disp( [ 'Dice Coeff - ' num2str( 1.0 / energy_i_t ) ] );
            disp( [ 'Energy Curve Smoothness - ' num2str( energy_s_t_r ) ] );
            disp( [ 'Mean Gateaux Derivative : ' num2str( dev_r ) ] );
            disp( [ 'Max Gateaux Derivative : ' num2str( gateaux_dev_r_max ) ] );   
            
            EnergyArr_r( t, t_r ) = energy_estimated;            
            DiceArr_r( t, t_r ) =  1.0 / energy_i_t ;
            
            if  gateaux_dev_r_max < 0.000005
                disp( 'Optimization Terminated : dE/dr < threshold' );
                break;        
            end

        end
        
        disp( '=============================================================================' );
        disp( '=============================================================================' );
        disp( 'Radius Update' );
        disp( [ 'Estimated Energy at Iter - ' num2str( t ) ' : ' num2str( energy_estimated ) ] );
        disp( [ 'Energy Overlapped - ' num2str( energy_i_t ) ] );
        disp( [ 'Dice Coeff - ' num2str( 1.0 / energy_i_t ) ] );
        disp( [ 'Energy Curve Smoothness - ' num2str( energy_s_t ) ] );
        disp( [ 'Mean Gateaux Derivative : ' num2str( dev_r ) ] );
        disp( [ 'Max Gateaux Derivative : ' num2str( gateaux_dev_r_max ) ] );   
        disp( '====================================================================' );
        disp( '====================================================================' );
        
        visInterp_ratio = 50;        
        [ recon_r, reconRGB_r ] = ReconOutputSkel( skelIdx_t, radiiList_t, bwImg, skelD, visInterp_ratio );
        imwrite( reconRGB_r, [ 'ellipse/SingleSkel_Result/ellipse_flow__t_' num2str( i ) '_singleSkel_R' num2str( t ) '.png' ] );        
    end

    %% Save Updated Skel and Radius
    fileSkelResultPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo_Updated.mat' ];
    fileRadiusResultPath = [ 'ellipse_flow__t_' num2str( i ) '_radius_Updated.mat' ];

    save( fileSkelResultPath, 'skelIdx_t' );
    save( fileRadiusResultPath, 'radiiList_t' );
    
    %% Reconstruction 
    % Reconstruct Object Surface by Skeleton + Radii Function
    visInterp_ratio = 50;        
    [ recon, reconRGB ] = ReconOutputSkel( skelIdx_t, radiiList_t, bwImg, skelD, visInterp_ratio );
    imwrite( reconRGB, [ 'ellipse/SingleSkel_Result/ellipse_flow__t_' num2str( i ) '_singleSkel.png' ] );
%}
end