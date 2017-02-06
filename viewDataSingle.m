close all
clear all

addpath '../../utilities/matlab/'

% load source data
% SourcePts contains the coordinates of the vertices, SourceEdges contains a connectivity matrix for the definition of the edges



for i = 0:0
    filename = [ 'ellipse_flow__t_' num2str( i ) ];
    [SourcePts SourceEdges] = VTKPolyDataReader([ filename '.vtk']);
    
    ellipseImg_p = zeros( [ 700, 1200 ] );
    
    SourcePts_s = [];
    
    % Change VTK to BW Image
    for k = 1:size( SourcePts, 1 )
        pt_x = SourcePts( k, 1 );
        pt_y = SourcePts( k, 2 );
        
        pt_x_o = pt_x + 80;
        pt_y_o = pt_y + 30;
        
        pt_x_s = pt_x_o * 5;
        pt_y_s = pt_y_o * 5;
        
        pt_x_i = round( pt_x_s ) + 200;
        pt_y_i = round( pt_y_s ) + 200;
        
        ellipseImg_p( pt_y_i, pt_x_i ) = 1;
        
        SourcePts_s = [ SourcePts_s; pt_x, pt_y ];
    end
    
    figure;
    imshow( ellipseImg_p );
    
    % Interpolate between points
    for k = 1:size( SourceEdges, 1 )
        line_x = ( SourcePts( SourceEdges( k, : ), 1 ) + 80 ) .* 5 + 200;
        line_y = ( SourcePts( SourceEdges( k, : ), 2 ) + 30 ) .* 5 + 200;
        
        if ( sqrt( ( line_x( 2 ) - line_x( 1 ) )^2 + ( line_y( 2 ) - line_y( 1 ) )^2 ) > 1.0 )
            if( line_x( 2 ) > line_x( 1 ) )
                xq = line_x(1):0.01:line_x(2);
                yq = interp1( line_x, line_y, xq );
            else
                xq = line_x(1):-0.01:line_x(2);
                yq = interp1( line_x, line_y, xq );
            end
            
            for t = 1:size( xq(:), 1 )
                ellipseImg_p( round( yq( t ) ), round( xq( t ) ) ) = 1;
            end
        end
    end
        
    figure;
    imshow( ellipseImg_p );
    title( 'Image Interp' );
    
    
    ellipseF = imfill( ellipseImg_p );
    figure;
    imshow( ellipseF );
    title( 'Filled' );
% 
%     save( [ filename '.mat' ], 'ellipseImg_p' );
%     save( [ filename '_filled.mat' ], 'ellipseF' );
%     save( [ filename '_Pts.mat' ], 'SourcePts_s' );
%     save( [ filename '_Edges.mat' ], 'SourceEdges' );

end
