files = dir('*.csv');

for i=1:length(files)
    load(files(i).name, '-ascii');
end

%%

Exn = zeros(3, 6, 6, 10000);
Eyn = zeros(3, 6, 6, 10000);

Exn(1, 1, :, :) = x2n_beta0_nu0_001;
Exn(1, 2, :, :) = x2n_beta0_nu0_0025;
Exn(1, 3, :, :) = x2n_beta0_nu0_005;
Exn(1, 4, :, :) = x2n_beta0_nu0_01;
Exn(1, 5, :, :) = x2n_beta0_nu0_02;
Exn(1, 6, :, :) = x2n_beta0_nu0_03;
Exn(2, 1, :, :) = x2n_beta0_002_nu0_001;
Exn(2, 2, :, :) = x2n_beta0_002_nu0_0025;
Exn(2, 3, :, :) = x2n_beta0_002_nu0_005;
Exn(2, 4, :, :) = x2n_beta0_002_nu0_01;
Exn(2, 5, :, :) = x2n_beta0_002_nu0_02;
Exn(2, 6, :, :) = x2n_beta0_002_nu0_03;
Exn(3, 1, :, :) = x2n_beta0_02_nu0_001;
Exn(3, 2, :, :) = x2n_beta0_02_nu0_0025;
Exn(3, 3, :, :) = x2n_beta0_02_nu0_005;
Exn(3, 4, :, :) = x2n_beta0_02_nu0_01;
Exn(3, 5, :, :) = x2n_beta0_02_nu0_02;
Exn(3, 6, :, :) = x2n_beta0_02_nu0_03;

Eyn(1, 1, :, :) = y2n_beta0_nu0_001;
Eyn(1, 2, :, :) = y2n_beta0_nu0_0025;
Eyn(1, 3, :, :) = y2n_beta0_nu0_005;
Eyn(1, 4, :, :) = y2n_beta0_nu0_01;
Eyn(1, 5, :, :) = y2n_beta0_nu0_02;
Eyn(1, 6, :, :) = y2n_beta0_nu0_03;
Eyn(2, 1, :, :) = y2n_beta0_002_nu0_001;
Eyn(2, 2, :, :) = y2n_beta0_002_nu0_0025;
Eyn(2, 3, :, :) = y2n_beta0_002_nu0_005;
Eyn(2, 4, :, :) = y2n_beta0_002_nu0_01;
Eyn(2, 5, :, :) = y2n_beta0_002_nu0_02;
Eyn(2, 6, :, :) = y2n_beta0_002_nu0_03;
Eyn(3, 1, :, :) = y2n_beta0_02_nu0_001;
Eyn(3, 2, :, :) = y2n_beta0_02_nu0_0025;
Eyn(3, 3, :, :) = y2n_beta0_02_nu0_005;
Eyn(3, 4, :, :) = y2n_beta0_02_nu0_01;
Eyn(3, 5, :, :) = y2n_beta0_02_nu0_02;
Eyn(3, 6, :, :) = y2n_beta0_02_nu0_03;