test_that("scrolling defers rendering and keeps model columns out of the queue grid", {
  skip_if_not_installed("DT")
  df_Data <- data.frame(ID = seq_len(5100), Date = as.Date("2020-01-01") + seq_len(5100))
  for (i in seq_len(20)) df_Data[[paste0("measure", i)]] <- seq_len(nrow(df_Data))
  result <- DetectDataAnomalies(df_Data, names(df_Data), id_var = "ID", methods = "isolation")
  before <- result
  queue <- CreateAnomalyReviewQueue(result, include_ranked = TRUE, top_n = 200)
  tables <- ExplorerTables(ExploreDataAnomalies(queue, IncludePlots = FALSE))
  grid <- Filter(function(x) "IssueID" %in% x$x$options$scidrMeta$columns, tables)[[1]]
  expect_equal(nrow(grid$x$data), nrow(queue$Issues))
  expect_gt(nrow(grid$x$data), 50)
  expect_false("measure1" %in% grid$x$options$scidrMeta$columns)
  expect_true(grid$x$options$deferRender)
  expect_equal(grid$x$options$scrollY, "500px")
  expect_true(grid$x$options$scrollX)
  expect_false(grepl("p", grid$x$options$dom, fixed = TRUE))
  expect_true("Scroller" %in% grid$x$extensions)
  expect_identical(result, before)
  expect_s3_class(ExploreDataAnomalies(queue, TableHeight = 400, IncludePlots = FALSE), "shiny.tag.list")
  expect_error(ExploreDataAnomalies(queue, TableHeight = 0), "positive")
})

test_that("history plotting remains optional", {
  skip_if_not_installed("DT")
  df_Data <- expand.grid(ID = 1:12, Time = 0:3)
  set.seed(43)
  df_Data$Value <- rnorm(nrow(df_Data))
  result <- DetectLongitudinalAnomalies(df_Data, "Value", "ID", "Time", methods = "change")
  testthat::local_mocked_bindings(PlotDataAnomalyTrajectories = function(...) stop("Plots invoked"))
  expect_s3_class(ExploreDataAnomalies(result), "shiny.tag.list")
  expect_error(ExploreDataAnomalies(result, IncludePlots = NA), "IncludePlots")
  if (!nzchar(system.file("htmlwidgets/lib/plotlyjs/plotly-latest.min.js", package = "plotly"))) skip("Plotly JavaScript unavailable")
  expect_s3_class(ExploreDataAnomalies(result, IncludePlots = TRUE), "shiny.tag.list")
})

test_that("the detail renderer escapes values and labels", {
  skip_if_not_installed("DT")
  if (!nzchar(Sys.which("node"))) skip("Node unavailable")
  # reviewValues were only ever attached to the removed technical record grid; the
  # escaping guarantee still has to hold for every value the detail panel prints.
  meta <- list(columns = c("Label", "Observed", "WhyFlagged"), labels = list(text = "<b>Label</b>"),
    methodNames = list(), evidence = list(), multipleAnalyses = FALSE)
  values <- list("<b>Label</b>", "<script>alert('x')</script>", "unusual")
  fn_Input <- tempfile(fileext = ".json"); fn_Script <- tempfile(fileext = ".js")
  on.exit(unlink(c(fn_Input, fn_Script)), add = TRUE)
  jsonlite::write_json(list(renderer = ScidrAnomalyDetailJS(), meta = meta, values = values), fn_Input, auto_unbox = TRUE)
  writeLines(c("const x = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));",
    "const render = eval('(' + x.renderer + ')');const html = render(x.values, x.meta);",
    "if (html.includes('<script>')) throw Error('Unescaped script tag');",
    "if (!html.includes('&lt;script&gt;')) throw Error('Missing escaped value');",
    "if (!html.includes('&lt;b&gt;Label&lt;/b&gt;')) throw Error('Missing escaped label');",
    "console.log('PASS');"), fn_Script)
  output <- system2(Sys.which("node"), c(shQuote(fn_Script), shQuote(fn_Input)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(output, "status")); expect_true("PASS" %in% output)
})
