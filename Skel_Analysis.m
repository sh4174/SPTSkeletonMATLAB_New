close all;
clear all;


for i = 0:0
    filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];

    bwImg = load( filePath );
    bw = bwImg.ellipseImg_p;
% 
%     figure;
%     imshow( bw );

    bwF = imfill( bw );
%     figure;
%     imshow( bwF );

    se = strel('disk', 40 );
    bwF_s = imclose( bwF, se );
    
    bwG = imgaussfilt( bwF_s, 5 );

    [ Gmag, Gdir ] = imgradient( bwG, 'intermediate' );
    
    Gdir_mask = Gdir;
    Gmag_mask = Gmag;
    
    Gdir_mask( bw == 0 ) = 0;
    Gmag_mask( bw == 0 ) = 0;
    
    figure;
    imshowpair( Gmag_mask, Gdir_mask, 'montage' );
    
    Gdir_s = imgaussfilt( Gdir, 0.3 );
    
    Gdir_mask = Gdir_s;
    Gmag_mask = Gmag;
    
    Gdir_mask( bw == 0 ) = 0;
    Gmag_mask( bw == 0 ) = 0;
    
    [ x, y ] = meshgrid( 1:size( Gdir_mask, 2 ), 1:size( Gdir_mask, 1 ) );
    
    u = zeros( size( x ) );
    v = zeros( size( x ) );
    
    for xx = 1:size( x, 1 )
    for yy = 1:size( x, 2 )
        x_u = x( xx, yy );
        y_u = y( xx, yy );
        
        theta = Gdir_mask( y_u, x_u );
        theta = theta / 180 * pi;
        
        u( xx, yy ) = cos( theta );
        v( xx, yy ) = sin( theta );
        
        if theta == 0
            u( xx, yy ) = 0;
            v( xx, yy ) = 0;
        end
    end
    end
    
    
    figure;
    quiver( x, y, u, v );
%     figure;
%     imshow( bwF_s );

    skel = bwmorph( bwF_s, 'skel', Inf );
 
%     figure;
%     imshow( skel );

    B = bwmorph(skel, 'branchpoints');
    E = bwmorph(skel, 'endpoints');
    [y,x] = find(E);
    B_loc = find(B);
    Dmask = false(size(skel));
    for k = 1:numel(x)
        D = bwdistgeodesic(skel,x(k),y(k));
        distanceToBranchPt = min(D(B_loc));
        Dmask(D < distanceToBranchPt) =true;
    end
    skelD = skel - Dmask;

%     figure;
%     imshow(skelD);
%     hold all;
%     [y,x] = find(B); plot(x,y,'ro')
end

