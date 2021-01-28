# clusterPower 0.7.0

- Major update with empirical (simulation-based) power calculation for many new study designs
- Now two shiny apps included; see ?runExample

# clusterPower 0.6.200

- Adding functions for difference in difference models
- Adding shiny app with runExample for GUI access to the package

# clusterPower 0.6.1

- Adding functions for closed-form solutions to normal, count (Poisson), and dichotomous outcomes that mirror power.t.test, etc.  
- Under the hood, documentation is migrated to Roxygen2, and all functions outside of base are explicitly called, via, e.g, stats::.
- Also, Nick Reich graciously transferring maintainership to Ken Kleinman.