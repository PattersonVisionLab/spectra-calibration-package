function results = runColorTestSuite(fName, opts)

    arguments
        fName               (1,1)       string      = ""
        opts.CoverageReport (1,1)       logical     = true
        opts.KeepFiles      (1,1)       logical     = false
    end

    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.codecoverage.CoverageReport

    % Run the suite in this function's directory ('test')
    userDirectory = pwd();
    testFolder = fullfile(fileparts(mfilename('fullpath')), "test");
    cd(testFolder);

    if fName == ""
        fName = pwd();
    end

    % Create the test suite and runner
    suite = testsuite(fName);
    runner = testrunner("textoutput");

    % Code coverage plugin
    if opts.CoverageReport
        reportFolder = fullfile(testFolder, "ColorCoverageReport");
        if ~exist(reportFolder, 'dir')
            mkdir(reportFolder);
        end
        p = CodeCoveragePlugin.forPackage("pattersonlab.core.color",...
            'IncludingSubpackages', true,...
            'Producing', CoverageReport(reportFolder));
        runner.addPlugin(p);
    end

    % Run the test suite
    results = runner.run(suite);

    % Open the code coverage results, if requested
    if opts.CoverageReport
        open(fullfile(reportFolder, 'index.html'));
    end

    % Clean up files produced by tests
    if ~opts.KeepFiles
        % TODO: Cleanup function
    end

    % Return to user's previous working directory
    cd(userDirectory);
