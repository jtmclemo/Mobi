function results = Mobi_run_all_tests()
    results = runtests('tests');
    disp(table(results));
end
