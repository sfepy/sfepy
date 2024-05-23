function [vals, vecs] = matlab_eig(mtx_filename, eigs_filename)

    data = load(mtx_filename);

    if data.verbose
        disp(data)
        disp(data.eigs_options)
    end

    if and(isequal(data.method, 'eig'), isequal(data.n_eigs, 'None'))
        if data.verbose
            disp('using eig()')
        end
        if isequal(data.B, 'None')
            if data.eigenvectors
                [vecs, vals] = eig(data.A, data.balance, data.algorithm);
                vals = diag(vals);
            else
                vals = eig(data.A, data.balance);
                vecs = 'None';
            end
        else
            if data.eigenvectors
                [vecs, vals] = eig(data.A, data.B, data.algorithm);
                vals = diag(vals);
            else
                vals = eig(data.A, data.B);
                vecs = 'None';
            end
        end
    else
        if data.verbose
            disp('using eigs()')
        end
        if isequal(data.n_eigs, 'None')
            data.n_eigs = size(data.A, 1);
        end
        if isequal(data.B, 'None')
            if data.eigenvectors
                [vecs, vals] = eigs(data.A, data.n_eigs, data.which, ...
                                    data.eigs_options);
                vals = diag(vals);
            else
                vals = eigs(data.A, data.n_eigs, data.which, ...
                            data.eigs_options);
                vecs = 'None';
            end
        else
            if data.eigenvectors
                [vecs, vals] = eigs(data.A, data.B, data.n_eigs, data.which, ...
                                    data.eigs_options);
                vals = diag(vals);
            else
                vals = eigs(data.A, data.B, data.n_eigs, data.which, ...
                            data.eigs_options);
                vecs = 'None';
            end
        end
    end

    save(eigs_filename, 'vals', 'vecs', '-v6');
