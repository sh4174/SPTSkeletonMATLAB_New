close all;
clear all;
addpath( 'ellipse' );

%% Read Data
i = 19;
filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];
filePtsPath = [ 'ellipse_flow__t_' num2str( i ) '_Pts.mat' ];
fileSkelDPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo.mat' ];
fileSkelPath = [ 'ellipse_flow__t_' num2str( i ) '_skel.mat' ];

% Read Image
bwImg = load( filePath );
bw = bwImg.ellipseImg_p;

% Read Surface points (Parametrized Curve)
pts = load( filePtsPath );
pts = pts.SourcePts_s;

% Read Estimated Skeleton
skel = load( fileSkelPath );
skel = skel.skel;

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
    skel_x = skelIdx( j, 2 );
    skel_y = skelIdx( j, 1 );
    radius = radiiList( j );
    
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

disp( gradRList );

figure;
subplot( 1, 3, 1 );
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), gradRList( 1:10:end, 1 ), gradRList( 1:10:end, 2 ) );
set(gca,'Ydir','reverse')
title( 'Gradients' );

% Calculate Corresponding Object Surface Normals
U_p = zeros( size( skelIdx ) );
U_n = zeros( size( skelIdx ) );

for j = 1:size( skelIdx, 1 )
    gradR_j = gradRList( j, : );
    normal_j = normalList( j, : );
    
    U_p_j = -gradR_j + sqrt( 1 - ( gradR_j( 1 ) * gradR_j( 1 ) + gradR_j( 2 ) * gradR_j( 2 ) ) ) .* normal_j;
    U_n_j = -gradR_j - sqrt( 1 - ( gradR_j( 1 ) * gradR_j( 1 ) + gradR_j( 2 ) * gradR_j( 2 ) ) ) .* normal_j;
    
    U_p( j, : ) = U_p_j;
    U_n( j, : ) = U_n_j;
end

subplot( 1, 3, 2 );
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), U_p( 1:10:end, 1 ), U_p( 1:10:end, 2 ) );
hold on
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), U_n( 1:10:end, 1 ), U_n( 1:10:end, 2 ) );
set(gca,'Ydir','reverse')
hold off
title( 'Surface Normals' );

subplot( 1, 3, 3 );
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), normalList( 1:10:end, 1 ), normalList( 1:10:end, 2 ) );
hold on
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), -normalList( 1:10:end, 1 ), -normalList( 1:10:end, 2 ) );
set(gca,'Ydir','reverse')
hold off
title( 'Normals' );

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
    
%     if skel( r, c ) > 0
%         reconRGB( r, c, 1 ) = 0;
%         reconRGB( r, c, 2 ) = 0;
%         reconRGB( r, c, 3 ) = 1;
%     end
    
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

% imwrite( reconRGB, [ 'ellipse/SingleSkel/ellipse_flow__t_' num2str( i ) '_Skels.png' ] );

figure;
imagesc( reconRGB );
hold on 
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), radiiList( 1:10:end ) .* U_p( 1:10:end, 1 ), radiiList( 1:10:end ) .* U_p( 1:10:end, 2 ) );
hold on
quiver( skelIdx( 1:10:end, 2 ), skelIdx( 1:10:end, 1 ), radiiList( 1:10:end ) .* U_n( 1:10:end, 1 ), radiiList( 1:10:end ) .* U_n( 1:10:end, 2 ) );
hold off
title( [ 'reconstructed : Img ' num2str( i ) ] );

