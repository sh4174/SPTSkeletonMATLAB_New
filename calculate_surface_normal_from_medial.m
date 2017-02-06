function [ U_p, U_n ] = calculate_surface_normal_from_medial( skelIdx, gradRList, normalList )

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

end

