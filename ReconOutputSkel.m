function [ recon, reconRGB ] = ReconOutputSkel( skelIdx_t, radiiList_t, bwImg, skelD, visInterp_ratio )
    recon = zeros( size( bwImg ) );
    skelIdx_t_int_r = interp( skelIdx_t( :, 1 ), visInterp_ratio );
    skelIdx_t_int_c = interp( skelIdx_t( :, 2 ), visInterp_ratio );
    
    radiiList_t_int = interp( radiiList_t, visInterp_ratio );
    
    skelIdx_t_int = [ skelIdx_t_int_r( 1:(end-visInterp_ratio+1), 1 ), skelIdx_t_int_c( 1:(end-visInterp_ratio+1), 1 ) ];           
    radiiList_t_int = radiiList_t_int( 1:(end-visInterp_ratio+1), 1 );

    for j = 1:size( skelIdx_t_int, 1 )
        skel_r = round( skelIdx_t_int( j, 1 ) );
        skel_c = round( skelIdx_t_int( j, 2 ) );

        radius = round( radiiList_t_int( j ) );
        recon = MidpointCircle( recon, radius, skel_r, skel_c, 1 );    
    end

    reconRGB = double( cat( 3, recon, recon, recon ) );
    lineWidth = 1;
    
    
    for j = 1:size( skelIdx_t_int, 1 )
        r = round( skelIdx_t_int( j, 1 ) );
        c = round( skelIdx_t_int( j, 2 ) );
        
        for dr = -lineWidth:lineWidth
            for dc = -lineWidth:lineWidth
                reconRGB( r + dr, c + dc, 1 ) = 0;
                reconRGB( r + dr, c + dc, 2 ) = 1;
                reconRGB( r + dr, c + dc, 3 ) = 0;
            end
        end
    end
    
    for r = 1:size( skelD, 1 )
    for c = 1:size( skelD, 2 )
        if bwImg( r, c ) > 0
            for dr = -lineWidth:lineWidth
                for dc = -lineWidth:lineWidth
                    reconRGB( r + dr, c + dc, 1 ) = 1;
                    reconRGB( r + dr, c + dc, 2 ) = 0;
                    reconRGB( r + dr, c + dc, 3 ) = 0;
                end
            end
        end
    end
    end
end

