Followed the guidelines from : https://r-pkgs.org/check.html

## Output of devtools::check()

There were no ERRORs, WARNINGs or NOTEs

## Output of devtools::check_win_devel()

checking CRAN incoming feasibility ... [14s] NOTE
Maintainer: 'Alexandre Jacquemain <aljacquemain@gmail.com>'

New submission

## Output of usethis::use_github_action_check_standard()

R-CMD-check / macos-latest (release)
Succeeded in 5 minutes and 44 seconds 

R-CMD-check / windows-latest (release)
Succeeded in 7 minutes and 45 seconds 

R-CMD-check / ubuntu-latest (release)
Succeeded in 4 minutes and 56 seconds 

## Output of devtools::check_rhub()

The two following notes appeared on Windows Server 2022, R-devel, 64 bit. I am unable to reproduce the first note on any other platform and I assume it is a bug. Concerning the second, as noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this is most likely a bug and can be ignored.

Checking HTML version of manual ... NOTE
Skipping checking math rendering: package 'V8' unavailable

checking for detritus in the temp directory ... NOTE
Found the following files/directories: 'lastMiKTeXException'


