# glorenz

<details>

* Version: 0.1.0
* GitHub: NA
* Source code: https://github.com/cran/glorenz
* Date/Publication: 2025-06-04 12:00:05 UTC
* Number of recursive dependencies: 62

Run `revdepcheck::revdep_details(, "glorenz")` for more info

</details>

## Newly broken

*   checking examples ... ERROR
    ```
    Running examples in ‘glorenz-Ex.R’ failed
    The error most likely occurred in:
    
    > ### Name: tlorenz
    > ### Title: Compute relevant probabilities and estimates for selecting
    > ###   performance criteria
    > ### Aliases: tlorenz
    > 
    > ### ** Examples
    > 
    > df_samp <- data.frame(x1 = rnorm(500, mean = 5, sd = 2),newwts = rep(1, 500))
    > df_samp2 <- data.frame(x1 = rnorm(500, mean = 4.5, sd = 2),newwts = rep(1, 500))
    > p_vals <- seq(0, 1, length.out = 100)
    > lc_vals <- tlorenz(p_vals, d1 = df_samp, group = "x1", d2 = df_samp2)
    Error in Lorenz.curve(y = d1[[group]], graph = FALSE, na.rm = TRUE, ties.method = c("mean",  : 
      unused argument (graph = FALSE)
    Calls: tlorenz
    Execution halted
    ```

*   checking whether package ‘glorenz’ can be installed ... WARNING
    ```
    Found the following significant warnings:
      Note: possible error in 'Lorenz.curve(y = d1[[group]], ': unused argument (graph = FALSE) 
      Note: possible error in 'Lorenz.curve(y = d2[[group]], ': unused argument (graph = FALSE) 
    See ‘/Users/Jacquemain/Library/CloudStorage/OneDrive-UCL/LorenzRegression/revdep/checks.noindex/glorenz/new/glorenz.Rcheck/00install.out’ for details.
    Information on the location(s) of code generating the ‘Note’s can be
    obtained by re-running with environment variable R_KEEP_PKG_SOURCE set
    to ‘yes’.
    ```

*   checking R code for possible problems ... NOTE
    ```
    tlorenz: possible error in Lorenz.curve(y = d1[[group]], graph = FALSE,
      na.rm = TRUE, ties.method = c("mean", "random"), seed = NULL, weights
      = d1[[newwts]]): unused argument (graph = FALSE)
    tlorenz: possible error in Lorenz.curve(y = d2[[group]], graph = FALSE,
      na.rm = TRUE, ties.method = c("mean", "random"), seed = NULL, weights
      = d2[[newwts]]): unused argument (graph = FALSE)
    ```

