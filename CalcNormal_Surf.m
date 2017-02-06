close all;
clear all;

addpath( 'ellipse' );

for i = 0:0
    filePath = [ 'ellipse_flow__t_' num2str( i ) '.mat' ];
    filePtsPath =  [ 'ellipse_flow__t_' num2str( i ) '_Pts.mat' ];
    
    bwImg = load( filePath );
    bw = bwImg.ellipseImg_p;
    
    pts = load( filePtsPath );
    pts = pts.SourcePts_s;

    [ tangents, normals ] = calc_curve_tangent_normal( pts );
   
    figure;
    quiver( pts( :, 1 ),  pts( :, 2 ), normals( :, 1 ), normals( :, 2 ) );    
end
