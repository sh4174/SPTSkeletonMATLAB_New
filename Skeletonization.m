close all;
clear all;


for i = 0:19
    filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];

    bwImg = load( filePath );
    bw = bwImg.ellipseImg_p;

    figure;
    imshow( bw );

    bwF = imfill( bw );
    figure;
    imshow( bwF );

    se = strel('disk', 40 );
    bwF_s = imclose( bwF, se );

    figure;
    imshow( bwF_s );

    skel = bwmorph( bwF_s, 'skel', Inf );

    figure;
    imshow( skel );

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

    figure;
    imshow(skelD);
    hold all;
    [y,x] = find(B); plot(x,y,'ro')
    
    fileSkelPath = [ 'ellipse_flow__t_' num2str( i ) '_skel.mat' ];
    fileSkelSPath = [ 'ellipse_flow__t_' num2str( i ) '_skel_0topo.mat' ];
    
    save( fileSkelPath, 'skel' );
    save( fileSkelSPath, 'skelD' );
end

