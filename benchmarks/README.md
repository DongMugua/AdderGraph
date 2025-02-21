# Benchmarks

Benchmarks are decomposed into two main steps, first the filter generation (see `filtergeneration\` for the scripts) which outputs a first set of plots in `results\plots\` and calls that are meant to be executed using [FloPoCo](https://flopoco.gforge.inria.fr/). Second, scripts in `synthesis\` launch the calls and extract relevant results. Finally from these results we output new plots.
