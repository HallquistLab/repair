# Package website

The pkgdown site at <https://hallquistlab.github.io/repair/> documents the current
`main` development version. It includes a short orientation/installation page,
a complete executable getting-started guide, all five articles, task-oriented
function reference, search, and NEWS. No generated HTML is committed.

## Sources and editing workflow

- `pkgdown/index.md`: concise landing page, example, function roles, guide chooser,
  and interpretation boundaries. The first-model example follows the executable
  README/workflow examples; keep them consistent when changing the interface.
- `_pkgdown.yml`: canonical URL, Bootstrap 5 theme, navigation, article grouping,
  and reference sections. Every documented topic and vignette is listed.
- `R/` roxygen comments and generated `man/*.Rd`: reference content. Update the
  roxygen source and regenerate help when changing function documentation.
- `vignettes/*.Rmd`: executable articles shared with the package documentation.
- `NEWS.md`: changelog. `DESCRIPTION` supplies authors, license, website,
  repository, and issue links.

The site leads with the public workflow and stable table extractors. Advanced
matrix interfaces remain available in their own reference section. Posterior
reliability, scoring-error covariance, and structural coefficient uncertainty are
distinguished explicitly. There is no claim that conditional standard errors
propagate uncertainty from measurement-model fitting.

Audit findings addressed here: useful guides were only discoverable as repository
source; reference topics lacked a task-based website index; there was no site
build/deployment check; and the main workflow guide retained an obsolete sentence
rejecting multiple measurement models. That sentence now points to the supported
explicit cross-model covariance policies.

## Local build and verification

From the package root, with the package build toolchain and Pandoc installed:

```r
install.packages("pak")  # once, if needed
pak::pkg_install(c("deps::.", "any::pkgdown"), dependencies = TRUE)
pkgdown::check_pkgdown()
pkgdown::build_site_github_pages(install = TRUE, new_process = TRUE)
```

This installs the current package in a temporary library before executing
reference examples and every article. Check the result:

```sh
Rscript .github/scripts/check-site.R
python3 -m http.server 8000 --directory docs
```

For search results to stay on localhost during a local preview, build with a
local URL override (do not change the committed canonical URL):

```r
pkgdown::build_site_github_pages(install = TRUE, new_process = TRUE,
  override = list(url = "http://localhost:8000/"))
```

The production build and downloaded CI artifacts intentionally use the public
URL for search results. Open <http://localhost:8000> and inspect the home page, a reference page, and an
article at desktop and narrow widths. Check search and the article navigation.
The verification script requires every reference/article page, search data, the
404 page, and `.nojekyll`; it also checks local hyperlinks, HTML anchors, and
image/script/stylesheet files throughout the generated site, plus search-result
pages and anchors. External URLs are
not fetched by this check, so external outages do not fail site validation.

`docs/`, `pkgdown/`, and `_pkgdown.yml` are excluded from R source-package builds;
`docs/` is also Git-ignored. pkgdown is a website build tool, not a package runtime
or test dependency. The existing package-check workflow remains separate.

## Continuous integration and publishing

`.github/workflows/pkgdown.yaml` builds on pull requests, pushes to `main`, and
manual dispatch. CI installs the local package through setup-r-dependencies,
then builds against that installed checkout with examples/articles enabled.
It runs the link/content check and uploads `pkgdown-site` as a downloadable HTML
preview for 14 days. Download/extract the artifact and serve it locally to review
a PR; PRs do not publish over the public site.

Only a successful build on `main` uploads a Pages artifact and runs the deployment
job. A dispatch on another branch builds a preview only. The deployment job uses
GitHub's official Pages action and the `github-pages` environment, with
`pages: write` and `id-token: write` confined to that job. The build has read-only
repository permissions. Actions are pinned to verified upstream revisions.
No personal publishing token, generated-content branch, or custom domain is used.

### One-time repository setup

A repository administrator/maintainer or a user with permission to manage Pages
must select **Settings → Pages → Build and deployment → Source → GitHub Actions**.
The equivalent API setup is `POST /repos/HallquistLab/repair/pages` with
`{"build_type":"workflow"}` when no site exists. Review any `github-pages`
environment rules so that deployments from `main` are allowed. Enabling this
setting does not publish a site by itself: the workflow must first reach `main`
and complete successfully. The public URL is
<https://hallquistlab.github.io/repair/>.

After merging, verify both **Build documentation** and **Deploy documentation**
in the pkgdown workflow, then open the public URL. If deployment fails before
publishing, check the Pages source setting, environment rules, and job logs.
Content failures should be resolved in sources and rechecked; do not bypass the
build/link check to deploy. To rebuild without source changes, dispatch the
pkgdown workflow on `main`.

The workflow follows [pkgdown's site builder](https://pkgdown.r-lib.org/reference/build_site.html)
and [GitHub's Pages workflow contract](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).
