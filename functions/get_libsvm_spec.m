function libsvmstring = get_libsvm_spec(cost,gamma,epsilon)
    error('this is outmoded. delete this')
    svr_type = 3; % hard coded epsilon-SVR (Zhang et al., 2014)
    kernel_type = 2; % hard coded kernel type -- RBF (Zhang et al., 2014)
    % for libsvm, standardize is done beforehand... not in the code flag.
    % note -e epsilon is not the loss-function epsilon, which is -p...
    libsvmstring = sprintf('-s %.4f -t %.4f -c %.4f -g %.4f -p %.4f -q',svr_type,kernel_type,cost,gamma,epsilon);
