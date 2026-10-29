# Validate the built site without depending on external websites being available.
# Run from the repository root after pkgdown::build_site_github_pages().
site <- normalizePath("docs", mustWork = TRUE)
required <- c(
  "index.html", "404.html", "reference/index.html", "articles/index.html",
  "news/index.html", "search.json", ".nojekyll",
  paste0("reference/", sub("\\.Rd$", ".html", list.files("man", pattern = "\\.Rd$"))),
  paste0("articles/", sub("\\.Rmd$", ".html", list.files("vignettes", pattern = "\\.Rmd$")))
)
missing <- required[!file.exists(file.path(site, required))]
if (length(missing)) stop("Missing website pages/assets: ", paste(missing, collapse = ", "))

files <- list.files(site, pattern = "\\.html$", recursive = TRUE, full.names = TRUE)
pages <- setNames(lapply(files, xml2::read_html), files)
anchors <- lapply(pages, function(page) {
  unique(c(xml2::xml_attr(xml2::xml_find_all(page, "//*[@id]"), "id"),
           xml2::xml_attr(xml2::xml_find_all(page, "//a[@name]"), "name")))
})
problems <- character()
for (file in files) {
  links <- xml2::xml_find_all(pages[[file]], "//a/@href | //img/@src | //script/@src | //link/@href")
  links <- unique(xml2::xml_text(links))
  # External URLs, data URIs, and empty fragment controls are not local links.
  links <- links[nzchar(links) & links != "#" &
                   !grepl("^([A-Za-z][A-Za-z0-9+.-]*:|//)", links)]
  for (link in links) {
    path <- URLdecode(sub("[?#].*$", "", link))
    fragment <- if (grepl("#", link, fixed = TRUE)) URLdecode(sub("^[^#]*#", "", link)) else ""
    target <- if (!nzchar(path)) file else if (startsWith(path, "/repair/")) {
      file.path(site, substring(path, nchar("/repair/") + 1L))
    } else if (startsWith(path, "/")) {
      file.path(site, substring(path, 2L))
    } else file.path(dirname(file), path)
    if (dir.exists(target)) target <- file.path(target, "index.html")
    label <- paste(substring(file, nchar(site) + 2L), "->", link)
    if (!file.exists(target)) {
      problems <- c(problems, paste("Missing target:", label))
      next
    }
    target <- normalizePath(target, mustWork = TRUE)
    if (nzchar(fragment) && endsWith(target, ".html") &&
        !fragment %in% anchors[[target]]) {
      problems <- c(problems, paste("Missing anchor:", label))
    }
  }
}
# Search results use the site's canonical URL; validate their local destinations
# too, so publication cannot introduce search links to missing pages/anchors.
search <- jsonlite::fromJSON(file.path(site, "search.json"), simplifyVector = FALSE)
# pkgdown can emit empty {"path": []} placeholders for unindexed sections.
search <- Filter(function(entry) !(identical(names(entry), "path") &&
                                   length(entry$path) == 0L), search)
if (!length(search) || any(!vapply(search, function(entry) {
  is.character(entry$path) && length(entry$path) == 1L && nzchar(entry$path)
}, logical(1)))) stop("Website search index is missing or malformed.")
search <- data.frame(
  path = vapply(search, function(entry) entry$path, character(1)),
  id = vapply(search, function(entry) if (is.null(entry$id)) NA_character_ else entry$id, character(1))
)
base_url <- sub("/reference/?$", "/", yaml::read_yaml(file.path(site, "pkgdown.yml"))$urls$reference)
for (i in seq_len(nrow(search))) {
  if (!startsWith(search$path[i], base_url)) {
    problems <- c(problems, paste("Search URL outside site:", search$path[i]))
    next
  }
  path <- substring(search$path[i], nchar(base_url) + 1L)
  target <- file.path(site, URLdecode(path))
  if (!file.exists(target)) {
    problems <- c(problems, paste("Missing search target:", path))
  } else if (!is.na(search$id[i]) && nzchar(search$id[i]) &&
             !search$id[i] %in% anchors[[normalizePath(target)]]) {
    problems <- c(problems, paste("Missing search anchor:", path, search$id[i]))
  }
}
if (length(problems)) stop(paste(c("Broken website links:", problems), collapse = "\n"))
cat("Validated ", length(files), " HTML pages, required content, internal links/assets, and search targets.\n", sep = "")
