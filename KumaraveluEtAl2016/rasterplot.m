function [] = rasterplot()

  load('simulate_network_model')

  strings = {'Ctx'; 'GPe'; 'GPi'; 'STN'; 'D1'; 'D2'; 'TH'};
  aps = [Cor_APs, GPe_APs, GPi_APs, STN_APs, Striat_APs_dr, Striat_APs_indr, TH_APs];
  len = [20, 10, 10, 10, 10, 10, 10];
  n = 7;
  k = 0;

  for j = 1:n
    ax = subplot(n, 1, j);
    times = [];
    dots = [];
    nspk = 0;
    for i = 1:len(j)
      k = k + 1;
      ap = aps(k);
      times = [times; ap.times];
      dots = [dots; ones(length(ap.times), 1) .* i];
      nspk = nspk + length(ap.times);
    end
    plot(ax, times, dots, 'LineStyle', 'none', 'Marker', '.');
    xlim(ax, [0,2])
    ylim(ax, [0,len(j)])
    ylabel(ax, strings{j})
    disp(strings{j}); disp(nspk / 2. / len(j))
  end
end