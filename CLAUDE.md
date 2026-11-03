# Repository rules

For the DANCE app, read `DANCE/CLAUDE.md` before changing anything.

Its first rule applies to every change: the exported R script must replicate
every analysis without the app. Any change that can alter a result goes into
the script in the same commit, with a round-trip test showing that the
script's numbers equal the app's.
