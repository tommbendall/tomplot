When running the tomplot test suite, command line arguments can be passed to
show and save the plots produced.

* Add --show_plots to display the plots as the test suite is run.
* Add --save_plots to save the plots in tests/tmp_plots. This allows the plots
  to be compared with the known good plots in tests/plots. The plots in this
  directory should not be added to the git repository.
* Add --overwrite_plots to save these plots as new known good plots. Only do
  this when happy with the plots!

This text file also acts a placeholder file for this directory to ensure that it
is created and stored in the repository.