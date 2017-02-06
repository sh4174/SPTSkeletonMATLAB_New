close all;
clear all;
addpath( 'ellipse' );

%% Read Data

for i = 0:19
    filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];
    filePtsPath = [ 'ellipse_flow__t_' num2str( i ) '_Pts.mat' ];
    fileSkelDPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo.mat' ];

    % Read Image
    bwImg = load( filePath );
    bw = bwImg.ellipseImg_p;

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

    % Reconstruct Object Surface by Skeleton + Radii Function
    recon = zeros( size( bw ) );

    for j = 1:size( skelIdx, 1 )
        skel_r = skelIdx( j, 1 );
        skel_c = skelIdx( j, 2 );

        radius = distSurf( skel_r, skel_c );
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

    imwrite( reconRGB, [ 'ellipse/SingleSkel/ellipse_flow__t_' num2str( i ) '_singleSkel.png' ] );
    
%     figure;
%     imshow( reconRGB );
%     title( [ 'reconstructed : Img ' num2str( i ) ] );
end


