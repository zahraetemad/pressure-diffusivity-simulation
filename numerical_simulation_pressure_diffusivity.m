% Pressure diffusivity in a 2D porous medium.
%
% assumptions:
%   - top-left boundary is treated as a known pressure inlet.
%   - bottom-right boundary is treated as a known pressure outlet.
%   - other boundaries are no-flow boundaries.
%


clear; close all; clc;

params = defaultParameters();
params = readUserInputs(params);
params = prepareNumerics(params);

fprintf('\nPressure diffusivity coefficient alpha = %.6g m^2/s\n', params.alpha);
fprintf('Using dx = %.6g m, dy = %.6g m, dt = %.6g s, final time = %.6g s\n', ...
    params.dx, params.dy, params.dt, params.tFinal);
fprintf('Explicit stability number alpha*dt*(1/dx^2 + 1/dy^2) = %.6g\n', ...
    params.stabilityNumber);

resultsDir = makeResultsDirectory();

explicitResult = solveExplicit(params);
implicitResult = solveImplicit(params);

fprintf('\nExplicit method elapsed time: %.6f s\n', explicitResult.elapsedTime);
fprintf('Implicit method elapsed time: %.6f s\n', implicitResult.elapsedTime);

profileDifference = rmsDifference(explicitResult.profile.pressure, ...
    implicitResult.profile.pressure) / max(abs(params.pin - params.pout), eps);
fprintf('Relative RMS difference on A-A line at final time: %.6g\n', profileDifference);

plotContourSnapshots(explicitResult, 'Explicit FTCS', ...
    fullfile(resultsDir, 'explicit_pressure_snapshots.png'));
plotContourSnapshots(implicitResult, 'Implicit backward Euler', ...
    fullfile(resultsDir, 'implicit_pressure_snapshots.png'));
plotLineComparison(explicitResult, implicitResult, ...
    fullfile(resultsDir, 'aa_profile_explicit_vs_implicit.png'));

if params.runConvergenceStudy
    convergence = runGridConvergenceStudy(params);
    plotConvergenceStudy(convergence, fullfile(resultsDir, 'grid_convergence.png'));
end

fprintf('\nPlots saved in: %s\n', resultsDir);

function params = defaultParameters()
    params.kPerm = 1e-12;          % Permeability, m^2
    params.phi = 0.25;             % Porosity, -
    params.mu = 1e-3;              % Viscosity, Pa s
    params.c = 1e-9;               % Compressibility, Pa^-1

    params.Lx = 100;               % Domain length, m
    params.Ly = 50;                % Domain width, m
    params.nxCells = 20;           % Number of cells in x
    params.nyCells = 10;           % Number of cells in y
    params.tFinal = 10;            % Final time, s

    params.p0 = 0;                 % Initial pressure, Pa
    params.pin = 1e5;              % Inlet pressure, Pa
    params.pout = 0;               % Outlet pressure, Pa
    params.runConvergenceStudy = true;
end

function params = readUserInputs(params)
    params.Lx = promptPositiveScalar('Domain length L in x direction, m', params.Lx);
    params.Ly = promptPositiveScalar('Domain width W in y direction, m', params.Ly);
    params.nxCells = promptPositiveInteger('Number of grid cells in x direction', params.nxCells);
    params.nyCells = promptPositiveInteger('Number of grid cells in y direction', params.nyCells);

    if mod(params.nxCells, 2) ~= 0
        params.nxCells = params.nxCells + 1;
        fprintf('Adjusted x grid cells to %d so the L/2 boundary split lies on a grid node.\n', ...
            params.nxCells);
    end

    alpha = params.kPerm / (params.phi * params.mu * params.c);
    dx = params.Lx / params.nxCells;
    dy = params.Ly / params.nyCells;
    maxStableDt = 0.5 / (alpha * (1 / dx^2 + 1 / dy^2));
    defaultDt = 0.8 * maxStableDt;

    params.dtRequested = promptPositiveScalar('Time step dt, s', defaultDt);
    params.tFinal = promptPositiveScalar('Final simulation time, s', params.tFinal);
    params.runConvergenceStudy = promptLogical('Run grid convergence study', params.runConvergenceStudy);
end

function params = prepareNumerics(params)
    params.alpha = params.kPerm / (params.phi * params.mu * params.c);
    params.nx = params.nxCells + 1;
    params.ny = params.nyCells + 1;
    params.dx = params.Lx / params.nxCells;
    params.dy = params.Ly / params.nyCells;

    maxStableDt = 0.5 / (params.alpha * (1 / params.dx^2 + 1 / params.dy^2));
    if params.dtRequested > maxStableDt
        fprintf(['Requested dt = %.6g s is unstable for the explicit method. ', ...
            'Using %.6g s instead.\n'], params.dtRequested, 0.8 * maxStableDt);
        params.dtRequested = 0.8 * maxStableDt;
    end

    params.nSteps = max(1, ceil(params.tFinal / params.dtRequested));
    params.dt = params.tFinal / params.nSteps;
    params.rx = params.alpha * params.dt / params.dx^2;
    params.ry = params.alpha * params.dt / params.dy^2;
    params.stabilityNumber = params.rx + params.ry;
    params.x = linspace(0, params.Lx, params.nx);
    params.y = linspace(0, params.Ly, params.ny);
    params.midX = params.Lx / 2;
end

function result = solveExplicit(params)
    tic;
    P = initialPressure(params);
    snapshotSteps = chooseSnapshotSteps(params.nSteps);
    snapshots = cell(numel(snapshotSteps), 1);
    snapshotTimes = zeros(numel(snapshotSteps), 1);
    snapshotCount = 0;

    for step = 1:params.nSteps
        Pold = P;
        Pnew = Pold;

        for iy = 2:params.ny-1
            for ix = 2:params.nx-1
                Pnew(ix, iy) = Pold(ix, iy) ...
                    + params.rx * (Pold(ix+1, iy) - 2 * Pold(ix, iy) + Pold(ix-1, iy)) ...
                    + params.ry * (Pold(ix, iy+1) - 2 * Pold(ix, iy) + Pold(ix, iy-1));
            end
        end

        P = applyBoundaryConditions(Pnew, params);

        if any(step == snapshotSteps)
            snapshotCount = snapshotCount + 1;
            snapshots{snapshotCount} = P;
            snapshotTimes(snapshotCount) = step * params.dt;
        end
    end

    result.method = 'explicit';
    result.P = P;
    result.snapshots = snapshots;
    result.snapshotTimes = snapshotTimes;
    result.elapsedTime = toc;
    result.params = params;
    result.profile = extractAALine(P, params);
end

function result = solveImplicit(params)
    tic;
    P = initialPressure(params);
    [A, dirichletSource] = buildImplicitSystem(params);

    snapshotSteps = chooseSnapshotSteps(params.nSteps);
    snapshots = cell(numel(snapshotSteps), 1);
    snapshotTimes = zeros(numel(snapshotSteps), 1);
    snapshotCount = 0;

    for step = 1:params.nSteps
        rhs = interiorVector(P, params) + dirichletSource;
        interiorPressure = A \ rhs;
        P = insertInteriorVector(P, interiorPressure, params);
        P = applyBoundaryConditions(P, params);

        if any(step == snapshotSteps)
            snapshotCount = snapshotCount + 1;
            snapshots{snapshotCount} = P;
            snapshotTimes(snapshotCount) = step * params.dt;
        end
    end

    result.method = 'implicit';
    result.P = P;
    result.snapshots = snapshots;
    result.snapshotTimes = snapshotTimes;
    result.elapsedTime = toc;
    result.params = params;
    result.profile = extractAALine(P, params);
end

function P = initialPressure(params)
    P = params.p0 * ones(params.nx, params.ny);
    P = applyBoundaryConditions(P, params);
end

function P = applyBoundaryConditions(P, params)
    tol = 100 * eps(params.Lx);

    % No-flow side boundaries.
    P(1, 2:params.ny-1) = P(2, 2:params.ny-1);
    P(params.nx, 2:params.ny-1) = P(params.nx-1, 2:params.ny-1);

    % No-flow top-right and bottom-left boundaries.
    topNoFlow = params.x > params.midX + tol;
    bottomNoFlow = params.x < params.midX - tol;
    P(topNoFlow, params.ny) = P(topNoFlow, params.ny-1);
    P(bottomNoFlow, 1) = P(bottomNoFlow, 2);

    % No-flow corners with no pressure boundary.
    P(1, 1) = P(2, 2);
    P(params.nx, params.ny) = P(params.nx-1, params.ny-1);

    % Known-pressure top-left inlet and bottom-right outlet.
    topPressure = params.x <= params.midX + tol;
    bottomPressure = params.x >= params.midX - tol;
    P(topPressure, params.ny) = params.pin;
    P(bottomPressure, 1) = params.pout;
end

function [A, dirichletSource] = buildImplicitSystem(params)
    nInteriorX = params.nx - 2;
    nInteriorY = params.ny - 2;
    nUnknown = nInteriorX * nInteriorY;

    A = spalloc(nUnknown, nUnknown, 5 * nUnknown);
    dirichletSource = zeros(nUnknown, 1);

    for iy = 2:params.ny-1
        for ix = 2:params.nx-1
            row = interiorIndex(ix, iy, params);
            center = 1 + 2 * params.rx + 2 * params.ry;

            if ix + 1 <= params.nx - 1
                A(row, interiorIndex(ix+1, iy, params)) = -params.rx;
            else
                center = center - params.rx;
            end

            if ix - 1 >= 2
                A(row, interiorIndex(ix-1, iy, params)) = -params.rx;
            else
                center = center - params.rx;
            end

            if iy + 1 <= params.ny - 1
                A(row, interiorIndex(ix, iy+1, params)) = -params.ry;
            elseif isTopPressure(params.x(ix), params)
                dirichletSource(row) = dirichletSource(row) + params.ry * params.pin;
            else
                center = center - params.ry;
            end

            if iy - 1 >= 2
                A(row, interiorIndex(ix, iy-1, params)) = -params.ry;
            elseif isBottomPressure(params.x(ix), params)
                dirichletSource(row) = dirichletSource(row) + params.ry * params.pout;
            else
                center = center - params.ry;
            end

            A(row, row) = center;
        end
    end
end

function value = isTopPressure(x, params)
    value = x <= params.midX + 100 * eps(params.Lx);
end

function value = isBottomPressure(x, params)
    value = x >= params.midX - 100 * eps(params.Lx);
end

function index = interiorIndex(ix, iy, params)
    index = (iy - 2) * (params.nx - 2) + (ix - 1);
end

function values = interiorVector(P, params)
    values = zeros((params.nx - 2) * (params.ny - 2), 1);
    for iy = 2:params.ny-1
        for ix = 2:params.nx-1
            values(interiorIndex(ix, iy, params)) = P(ix, iy);
        end
    end
end

function P = insertInteriorVector(P, values, params)
    for iy = 2:params.ny-1
        for ix = 2:params.nx-1
            P(ix, iy) = values(interiorIndex(ix, iy, params));
        end
    end
end

function profile = extractAALine(P, params)
    [~, ixMid] = min(abs(params.x - params.midX));
    profile.xIndex = ixMid;
    profile.xCoordinate = params.x(ixMid);
    profile.distance = params.y(:);
    profile.pressure = P(ixMid, :).';
end

function snapshotSteps = chooseSnapshotSteps(nSteps)
    snapshotSteps = unique(max(1, round([0.25, 0.50, 0.75, 1.00] * nSteps)));
end

function convergence = runGridConvergenceStudy(baseParams)
    fprintf('\nRunning implicit grid convergence study...\n');
    refinementFactors = [1, 2, 4];
    profiles = cell(numel(refinementFactors), 1);
    dxValues = zeros(numel(refinementFactors), 1);
    errors = zeros(numel(refinementFactors), 1);
    labels = cell(numel(refinementFactors), 1);

    for i = 1:numel(refinementFactors)
        factor = refinementFactors(i);
        params = baseParams;
        params.nxCells = makeEven(baseParams.nxCells * factor);
        params.nyCells = max(2, baseParams.nyCells * factor);
        params.dtRequested = baseParams.dtRequested / factor^2;
        params = prepareNumerics(params);

        result = solveImplicit(params);
        profiles{i} = result.profile;
        dxValues(i) = params.dx;
        labels{i} = sprintf('%d x %d cells', params.nxCells, params.nyCells);
        fprintf('  %s: dx = %.6g m, dt = %.6g s, elapsed = %.4f s\n', ...
            labels{i}, params.dx, params.dt, result.elapsedTime);
    end

    fineProfile = profiles{end};
    for i = 1:numel(refinementFactors)-1
        coarseOnFine = interp1(profiles{i}.distance, profiles{i}.pressure, ...
            fineProfile.distance, 'linear');
        errors(i) = rmsDifference(coarseOnFine, fineProfile.pressure) ...
            / max(abs(baseParams.pin - baseParams.pout), eps);
    end
    errors(end) = NaN;

    convergence.profiles = profiles;
    convergence.dxValues = dxValues;
    convergence.errors = errors;
    convergence.labels = labels;
end

function plotContourSnapshots(result, plotTitle, outputPath)
    fig = figure('Name', plotTitle, 'Color', 'w');
    nPlots = numel(result.snapshots);

    for i = 1:nPlots
        subplot(2, 2, i);
        contourf(result.params.x, result.params.y, result.snapshots{i}.', 30, ...
            'LineColor', 'none');
        colorbar;
        axis equal tight;
        xlabel('x (m)');
        ylabel('y (m)');
        title(sprintf('t = %.3g s', result.snapshotTimes(i)));
    end

    addSuperTitle(plotTitle);
    saveFigure(fig, outputPath);
end

function plotLineComparison(explicitResult, implicitResult, outputPath)
    fig = figure('Name', 'A-A pressure profile', 'Color', 'w');
    plot(explicitResult.profile.distance, explicitResult.profile.pressure, ...
        'LineWidth', 1.5);
    hold on;
    plot(implicitResult.profile.distance, implicitResult.profile.pressure, ...
        '--', 'LineWidth', 1.5);
    hold off;
    grid on;
    xlabel('Distance along A-A, y (m)');
    ylabel('Pressure (Pa)');
    title(sprintf('A-A pressure profile at t = %.3g s', explicitResult.params.tFinal));
    legend('Explicit FTCS', 'Implicit backward Euler', 'Location', 'best');
    saveFigure(fig, outputPath);
end

function plotConvergenceStudy(convergence, outputPath)
    fig = figure('Name', 'Grid convergence study', 'Color', 'w');

    subplot(1, 2, 1);
    hold on;
    for i = 1:numel(convergence.profiles)
        plot(convergence.profiles{i}.distance, convergence.profiles{i}.pressure, ...
            'LineWidth', 1.5);
    end
    hold off;
    grid on;
    xlabel('Distance along A-A, y (m)');
    ylabel('Pressure (Pa)');
    title('Final A-A profiles');
    legend(convergence.labels, 'Location', 'best');

    subplot(1, 2, 2);
    valid = ~isnan(convergence.errors);
    loglog(convergence.dxValues(valid), convergence.errors(valid), 'o-', ...
        'LineWidth', 1.5);
    grid on;
    set(gca, 'XDir', 'reverse');
    xlabel('Grid spacing dx (m)');
    ylabel('Relative RMS difference vs finest grid');
    title('Grid convergence metric');

    saveFigure(fig, outputPath);
end

function addSuperTitle(textValue)
    if exist('sgtitle', 'file') == 2
        sgtitle(textValue);
    end
end

function saveFigure(fig, outputPath)
    if exist('exportgraphics', 'file') == 2
        exportgraphics(fig, outputPath, 'Resolution', 200);
    else
        saveas(fig, outputPath);
    end
end

function resultsDir = makeResultsDirectory()
    scriptPath = mfilename('fullpath');
    if isempty(scriptPath)
        projectRoot = pwd;
    else
        projectRoot = fileparts(scriptPath);
    end
    resultsDir = fullfile(projectRoot, 'results');
    if ~exist(resultsDir, 'dir')
        mkdir(resultsDir);
    end
end

function value = promptPositiveScalar(promptText, defaultValue)
    raw = input(sprintf('%s [%.6g]: ', promptText, defaultValue), 's');
    if isempty(strtrim(raw))
        value = defaultValue;
    else
        value = str2double(raw);
    end

    while ~isfinite(value) || value <= 0
        raw = input(sprintf('Please enter a positive number for %s: ', promptText), 's');
        value = str2double(raw);
    end
end

function value = promptPositiveInteger(promptText, defaultValue)
    value = round(promptPositiveScalar(promptText, defaultValue));
    while value < 2
        raw = input(sprintf('Please enter an integer >= 2 for %s: ', promptText), 's');
        value = round(str2double(raw));
    end
end

function value = promptLogical(promptText, defaultValue)
    if defaultValue
        defaultText = 'y';
    else
        defaultText = 'n';
    end

    raw = lower(strtrim(input(sprintf('%s? y/n [%s]: ', promptText, defaultText), 's')));
    if isempty(raw)
        value = defaultValue;
    else
        value = strcmp(raw(1), 'y') || strcmp(raw(1), '1');
    end
end

function value = makeEven(value)
    if mod(value, 2) ~= 0
        value = value + 1;
    end
end

function value = rmsDifference(a, b)
    diff = a(:) - b(:);
    value = sqrt(mean(diff.^2));
end
