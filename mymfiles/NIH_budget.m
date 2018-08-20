function NIH_budget(sals, calmonths)

far = 0.69;
fbr = [0.28 0.28 0.21];  % modify this as needed

disp(sprintf('person\t months\t Inst sal\t sal req\t fringe\t total\n'));
for i = 1:length(sals)
    salreq = (calmonths(i)/12*sals(i));
    fb = (salreq*fbr(i));
    tot = (salreq+fb);
    line(i,:) = [calmonths(i) round(sals(i)) round(salreq)   round(fb) round(tot)]; 
    disp(sprintf('guy%i\t %0.2f\t %i\t\t %i\t\t %i\t %i\n', i, line(i,:)));
end

subtot = sum(line);
disp(sprintf('tot:\t%0.2f\t %i\t\t %i\t\t %i\t %i\n', subtot));

direct = subtot(end);
indirect = direct*far;
total = direct+indirect;

disp(sprintf('direct:  %i\n', round(direct)));
disp(sprintf('indir:  %i\n', round(indirect)));
disp(sprintf('TOTAL:  %i\n', round(total)));
