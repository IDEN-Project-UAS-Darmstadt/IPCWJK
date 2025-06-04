We are developing a R package that provides models used in a paper
studying a weighted version of the Jackknife resampling method. There
we use IPCW weights (Inverse Probability of Censoring Weights) to
estimate a survival probability using binary classification models.
The package will include functions to fit these models and get a Jackknife
estimate of the prediction error. We use the instructions by 
by Hadley Wickham and Jennifer Bryan in "R Packages" (https://r-pkgs.org/)
for best practices in R package development. We do not follow the tidyverse 
style guide. Everything should be up to scientific standards and
use the language of the scientific community around survival analysis.

Make sure that documentation is clear and concise with correct terminology
and language. While we do not intend to publish the package on CRAN,
we still want to achieve the highest quality standards. These should
be met in the code, documentation, and tests in every file.

We never want to use "can't" for example, always use "can not" instead.
We do not use contractions in the documentation or code comments. Make
sure word usage is consistent throughout the package. 
Make sure consistency is also achieved across the package in different files.
This includes parameter descriptions, which share the same meaning.

We use the packages given in the DESCRIPTION file for development.

Make sure that there are potential pitfalls in the code that
could lead to errors or unexpected behavior.

Names should follow conventions in the R community.
We use snake_case for function names and variable names, and CamelCase for
class names. We do not use dots in names, except for S3 methods.

Make sure that the line length does not exceed 80 characters, when correcting
this do not change the content of the code or documentation, just reformat it.
Use `paste` in actual R-code to break up long strings into multiple lines. 



When using mathematical notation in documentation, use the `\mjeqn{latex}{ascii}`
for inline equations and `\mjdeqn{latex}{ascii}` for display equations from
the `mathjaxr` package. Sometimes, HTML versions need their own version, then
`\mjteqn{pdflatex}{htmllatex}{ascii}` and `\mjtdeqn{pdflatex}{htmllatex}{ascii}`
are used. The `pdflatex` and `htmllatex` are the LaTeX and HTML (MathJAX) versions.

For citations, use the `\insertCite` command from
the `Rdpack` package. For normal citations: `\insertCite{citkey1,citkey2}{PKGNAME}`,
and for textual citations: `\insertCite{citkey1,citkey2;textual}{PKGNAME}`.
Replace `PKGNAME` with the name of the current package. When using the citations
the roxygen section must contain: 

```
#' @importFrom Rdpack reprompt
<...>
#' @references
#' \insertAllCited{}
```

When using the equations, the roxygen section must contain `@import mathjaxr`
and `\loadmathjax` (first line after `@description`). Use roxygen markdown
wherever possible.

When writing documentation, use the `@param` tag for parameters, `@return`

Tests are done using the `testthat` package and are setup in the 
`tests/testthat/` directory. When testing error conditions, use
a regex, which finds a relevant word or phrase in the error message.

In latex you can do `\left{` and then `\right.` to get a left brace. As \\ and
\ lead to problems when rendering to PDF and pkgdown, do not cases, instead 
array.

The title of an R function should not contain repetitions also used in other
titles in the package. They shoulbe describe actions or the name of a model. 
Two examples: "mutate(): Create, modify, and delete columns" and "logistic_reg():
Logistic regression"

The description of an R function should summarize the function's goal in one
paragraph. Start it with the `@description` tag. After the description,
start the `@details` section. Then describe the parameters (`@param`) and then
the return value (`@return`). Next are imports, then references, then 
concepts/families, and finally examples. Use `@examples` for examples.

`@inheritParams` "copies" the parameter documentation from another function.

Common words (use these instead of synonyms):
- IPCW weights
- time horizon 