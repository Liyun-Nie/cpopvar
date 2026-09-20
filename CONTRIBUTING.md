# Contributing

Please open an issue before proposing a substantial change. Include a minimal
reproducible example, the R version, operating system, configuration fragment,
and the smallest input data that demonstrates the problem.

Contributions should:

- keep analysis controllers configuration-driven;
- avoid hard-coded column names for auxiliary user files;
- obtain plotting colors through the package color helper;
- document exported functions with roxygen2;
- add or update tests for changed behavior;
- avoid raw genomes, private data, credentials, absolute machine paths, and
  generated analysis results.

Before submitting a change, run:

```bash
R CMD build .
R CMD check --as-cran cpopvar_*.tar.gz
```

By contributing, you agree that your contribution is distributed under the
MIT License.
