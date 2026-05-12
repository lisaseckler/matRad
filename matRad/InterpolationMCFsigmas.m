%% Chat GPT interpolation MCF sigmas

E1 = arrayfun(@(x) x.energy, machine.data);  % 121x1
E2 = arrayfun(@(x) x.energy, MCF.data);  % 166x1

[E1, order] = sort(E1);
machine.data = machine.data(order);

n2 = numel(MCF.data);
idxLow  = zeros(n2,1);
idxHigh = zeros(n2,1);
w       = zeros(n2,1);

for j = 1:n2
    if E2(j) <= E1(1)
        idxLow(j) = 1;
        idxHigh(j) = 1;
        w(j) = 0;
    elseif E2(j) >= E1(end)
        idxLow(j) = numel(E1);
        idxHigh(j) = numel(E1);
        w(j) = 0;
    else
        idxHigh(j) = find(E1 >= E2(j), 1, 'first');
        idxLow(j)  = idxHigh(j) - 1;
        w(j) = (E2(j) - E1(idxLow(j))) / ...
               (E1(idxHigh(j)) - E1(idxLow(j)));
    end
end



%%
for j = 1:n2

    d2 = MCF.data(j).depths;    % 150x1

    % untere Energie
    d1a = machine.data(idxLow(j)).depths;
    s1a = machine.data(idxLow(j)).sigma;

    s2a = interp1(d1a, s1a, d2, 'linear', 'extrap');

    % obere Energie
    if idxLow(j) ~= idxHigh(j)
        d1b = machine.data(idxHigh(j)).depths;
        s1b = machine.data(idxHigh(j)).sigma;

        s2b = interp1(d1b, s1b, d2, 'linear', 'extrap');

        % Energie-Interpolation
        s2 = (1-w(j))*s2a + w(j)*s2b;
    else
        s2 = s2a;
    end

    MCF.data(j).sigma = s2;
end

%% New try
% 1) Energien extrahieren


E1 = arrayfun(@(x) x.energy, machine.data);
E2 = arrayfun(@(x) x.energy, MCF.data);

% Sortieren (SEHR wichtig)
[E1, order] = sort(E1);
machine.data = machine.data(order);

n2 = numel(MCF.data);

%% =========================
% 2) Gewichte definieren
%    (hier erstmal konstant – später optimierbar)
%% =========================

weights = ones(n2,1);   % später anpassbar

%% =========================
% 3) SOBP berechnen
% =========================

for j = 1:n2
    
    Ej = E2(j);
    d2 = MCF.data(j).depths(:);   % 150x1
    
    % ---- Energie-Interpolation vorbereiten ----
    
    if Ej <= E1(1)
        idxLow = 1;
        idxHigh = 1;
        wE = 0;
    elseif Ej >= E1(end)
        idxLow = numel(E1);
        idxHigh = numel(E1);
        wE = 0;
    else
        idxHigh = find(E1 >= Ej, 1, 'first');
        idxLow  = idxHigh - 1;
        
        wE = (Ej - E1(idxLow)) / ...
             (E1(idxHigh) - E1(idxLow));
    end
    
    % ---- Untere Energie ----
    
    d1a = machine.data(idxLow).depths(:);
    s1a = machine.data(idxLow).sigma(:);
    
    % Tiefe muss sortiert sein!
    [d1a, sortIdx] = sort(d1a);
    s1a = s1a(sortIdx);
    
    s2a = interp1(d1a, s1a, d2, 'pchip', 'extrap');
    
    % ---- Obere Energie (falls nötig) ----
    
    if idxLow ~= idxHigh
        
        d1b = machine.data(idxHigh).depths(:);
        s1b = machine.data(idxHigh).sigma(:);
        
        [d1b, sortIdx] = sort(d1b);
        s1b = s1b(sortIdx);
        
        s2b = interp1(d1b, s1b, d2, 'pchip', 'extrap');
        
        % Energie-Interpolation
        sInterp = (1 - wE)*s2a + wE*s2b;
        
    else
        sInterp = s2a;
    end
    
    % ---- Untere Energie-Extrapolation dämpfen ----
    
    if Ej < E1(1)
        scale = Ej / E1(1);
        scale = max(scale,0);
        sInterp = scale * sInterp;
    end
    
    % ---- Gewicht anwenden ----
    
    %Tiefe2.data(j).sigma = weights(j) * sInterp;
    MCF.data(j).sigma = sInterp;
    
end

%% =========================
% 4) SOBP aufsummieren
%% =========================

depthCommon = Tiefe2.data(1).depth;
SOBP = zeros(size(depthCommon));

for j = 1:n2
    SOBP = SOBP + Tiefe2.data(j).sigma;
end

%% =========================
% 5) Normieren
%% =========================

SOBP = SOBP / max(SOBP);

%% =========================
% 6) Plot
%% =========================

figure
plot(depthCommon, SOBP, 'LineWidth', 2)
xlabel('Tiefe')
ylabel('Normierte Dosis')
title('Kohlenstoff SOBP')
grid on
