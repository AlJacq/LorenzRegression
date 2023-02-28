Followed the guidelines from : https://r-pkgs.org/check.html

## Output of devtools::check()

There were no ERRORs, WARNINGs or NOTEs

## Output of devtools::check_win_devel()

Status: 1 NOTE
R Under development (unstable) (2023-02-27 r83911 ucrt)

## Output of usethis::use_github_action_check_standard()

Success of

R-CMD-check / macos-latest (release)
R-CMD-check / windows-latest (release)
R-CMD-check / ubuntu-latest (release)

## Output of devtools::check_rhub()

The two following notes appeared on Windows Server 2022, R-devel, 64 bit. I am unable to reproduce the first note on any other platform and I assume it is a bug. Concerning the second, as noted in [R-hub issue #503](https://github.com/r-hub/rhub/issues/503), this is most likely a bug and can be ignored.

Checking HTML version of manual ... NOTE
Skipping checking math rendering: package 'V8' unavailable

checking for detritus in the temp directory ... NOTE
Found the following files/directories: 'lastMiKTeXException'


