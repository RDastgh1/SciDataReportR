# Shared temporal declarations for intake and anomaly modeling.
ScidrTimeOverrides <- function(value, variables, argument) {
  if (is.null(value)) return(stats::setNames(character(), character()))
  if (!is.character(value) || is.null(names(value)) || anyNA(value) || anyNA(names(value)) ||
      any(!nzchar(trimws(value))) || any(!nzchar(names(value))) || anyDuplicated(names(value)) || any(!names(value) %in% variables))
    stop(argument, " must be a named character vector keyed by existing columns.")
  value
}

ScidrTimeUnit <- function(unit) {
  aliases <- list(secs = c("secs", "sec", "second", "seconds", "s"),
    mins = c("mins", "min", "minute", "minutes"), hours = c("hours", "hour", "hr", "hrs", "h"),
    days = c("days", "day", "d"), weeks = c("weeks", "week", "w"))
  unit <- tolower(trimws(unit))
  found <- names(aliases)[vapply(aliases, function(x) unit %in% x, logical(1))]
  if (length(found) != 1L) stop("Unknown TimeUnit: ", unit, ". Use secs, mins, hours, days or weeks.")
  found
}

ScidrTimeDeclarations <- function(data, codebook = NULL, time_types = NULL, time_units = NULL, time_formats = NULL, detection = FALSE) {
  types <- ScidrTimeOverrides(time_types, names(data), "time_types")
  units <- ScidrTimeOverrides(time_units, names(data), "time_units")
  formats <- ScidrTimeOverrides(time_formats, names(data), "time_formats")
  types[] <- tolower(trimws(types))
  if (any(!types %in% c("clock", "duration"))) stop("time_types values must be clock or duration.")
  if (length(units)) units[] <- vapply(units, ScidrTimeUnit, character(1))
  if (length(formats)) {
    formats[] <- toupper(trimws(formats))
    if (any(!formats %in% c("HMS", "HM", "MS"))) stop("time_formats must use HMS, HM or MS.")
  }
  Get <- function(variable, field) {
    if (is.null(codebook) || !all(c("Variable", field) %in% names(codebook))) return(NA_character_)
    value <- as.character(codebook[[field]][match(variable, codebook$Variable)])
    if (!length(value) || is.na(value) || !nzchar(trimws(value))) NA_character_ else trimws(value)
  }
  pieces <- list()
  ambiguous <- character()
  for (variable in names(data)) {
    x <- data[[variable]]
    explicit <- types[variable]
    book <- tolower(Get(variable, "Type"))
    attached <- attr(x, "scidr_type", exact = TRUE)
    kind <- NA_character_
    origin <- NA_character_
    if (!is.na(explicit)) { kind <- unname(explicit); origin <- "override" }
    else if (!is.na(book) && book %in% c("clock", "duration")) { kind <- book; origin <- "codebook" }
    else if (length(attached) == 1L && attached %in% c("clock", "duration")) { kind <- attached; origin <- "attached" }
    else if (detection && inherits(x, "difftime") && !inherits(x, "hms")) { kind <- "duration"; origin <- "difftime class" }
    else if (detection && inherits(x, "hms")) ambiguous <- c(ambiguous, variable)
    if (is.na(kind)) next
    unit <- if (!is.na(units[variable])) unname(units[variable]) else Get(variable, "TimeUnit")
    format <- if (!is.na(formats[variable])) unname(formats[variable]) else Get(variable, "TimeFormat")
    if (!is.na(unit)) unit <- ScidrTimeUnit(unit)
    if (!is.na(format)) {
      format <- toupper(format)
      if (!format %in% c("HMS", "HM", "MS")) stop("Invalid TimeFormat for ", variable, ": use HMS, HM or MS.")
    }
    if (kind == "clock" && identical(format, "MS")) stop("Clock times require HMS or HM: ", variable)
    pieces[[length(pieces) + 1L]] <- tibble::tibble(Variable = variable, TimeType = kind, TimeUnit = unit,
      TimeFormat = format, Source = origin, InputClass = paste(class(x), collapse = "/"),
      InputUnit = if (inherits(x, "difftime")) as.character(attr(x, "units")) else unit)
  }
  # An hms column may hold a time of day or an elapsed duration, and the two are
  # analysed differently. Rather than stopping, leave the column undeclared and
  # say so: it is then compared as plain seconds, which is correct for a duration
  # and merely loses the wrap-around handling for a clock. The names travel on
  # the result so the caller can perform that conversion -- without it the
  # column reaches a modelling backend as hms and is rejected outright, which
  # stops the whole method rather than costing it some accuracy.
  if (length(ambiguous)) warning("Ambiguous hms columns compared as plain seconds without clock-aware handling: ",
    paste(ambiguous, collapse = ", "),
    ". Declare Type = clock/duration in the codebook, run RevalueData(), or supply time_types = c(column = 'clock', other = 'duration').",
    call. = FALSE)
  output <- dplyr::bind_rows(c(list(tibble::tibble(Variable = character(), TimeType = character(), TimeUnit = character(),
    TimeFormat = character(), Source = character(), InputClass = character(), InputUnit = character())), pieces))
  unused <- setdiff(unique(c(names(units), names(formats))), output$Variable)
  if (length(unused)) stop("Temporal units/formats need a clock or duration declaration: ", paste(unused, collapse = ", "))
  attr(output, "scidr_ambiguous") <- ambiguous
  output
}

ScidrTimeSeconds <- function(x, declaration, strict = TRUE) {
  if (inherits(x, "difftime")) return(as.numeric(x, units = "secs"))
  if (inherits(x, c("Date", "POSIXt")) || is.factor(x) || is.logical(x) || !is.null(dim(x)))
    stop("Temporal inputs must be hms/difftime, numeric, or text; calendar dates/timestamps are not implicitly stripped.")
  if (is.numeric(x)) {
    if (is.na(declaration$TimeUnit)) stop("Numeric temporal values require TimeUnit or time_units.")
    multiplier <- c(secs = 1, mins = 60, hours = 3600, days = 86400, weeks = 604800)[[declaration$TimeUnit]]
    return(as.numeric(x) * multiplier)
  }
  if (!is.character(x)) stop("Unsupported temporal storage.")
  format <- declaration$TimeFormat
  if (is.na(format)) stop("Text temporal values require TimeFormat or time_formats: HMS, HM or MS.")
  result <- rep(NA_real_, length(x))
  bad <- integer()
  for (i in which(!is.na(x))) {
    text <- trimws(x[i])
    pattern <- if (format == "HM") "^[+-]?[0-9]+:[0-9]{1,2}$" else if (format == "MS")
      "^[+-]?[0-9]+:[0-9]{1,2}(\\.[0-9]+)?$" else "^[+-]?[0-9]+:[0-9]{1,2}:[0-9]{1,2}(\\.[0-9]+)?$"
    if (!grepl(pattern, text)) {bad <- c(bad, i); next}
    sign <- if (startsWith(text, "-")) -1 else 1
    parts <- as.numeric(strsplit(sub("^[+-]", "", text), ":", fixed = TRUE)[[1]])
    valid <- if (format == "HMS") parts[2] < 60 && parts[3] < 60 else parts[2] < 60
    if (!valid || (declaration$TimeType == "clock" && grepl("^[+-]", text))) {bad <- c(bad, i); next}
    result[i] <- sign * sum(parts * switch(format, HMS = c(3600, 60, 1), HM = c(3600, 60), MS = c(60, 1)))
  }
  if (length(bad) && strict) stop("Unparseable temporal values at rows: ", paste(bad, collapse = ", "))
  if (length(bad)) attr(result, "parse_failures") <- bad
  result
}

ScidrTimeFormat <- function(seconds, kind) {
  if (kind == "duration") return(vapply(seconds, function(x) if (is.na(x)) NA_character_ else paste0(format(signif(x, 8), trim = TRUE, scientific = FALSE), " seconds"), character(1), USE.NAMES = FALSE))
  vapply(seconds, function(x) {
    if (is.na(x) || !is.finite(x)) return(NA_character_)
    x <- round(x, 6)
    hour <- floor(x / 3600); minute <- floor((x - hour * 3600) / 60); second <- x - hour * 3600 - minute * 60
    text <- sprintf("%02d:%02d:%09.6f", as.integer(hour), as.integer(minute), second)
    sub("\\.?0+$", "", text)
  }, character(1), USE.NAMES = FALSE)
}

ScidrTimeMissing <- function(x, tokens) {
  missing <- !is.na(x) & trimws(as.character(x)) %in% tokens
  if (inherits(x, "difftime")) {
    codes <- suppressWarnings(as.numeric(tokens)); codes <- codes[!is.na(codes)]
    missing <- missing | (!is.na(x) & as.numeric(x) %in% codes)
  } else if (is.numeric(x)) missing <- missing | (!is.na(x) & x %in% suppressWarnings(as.numeric(tokens)))
  x[missing] <- NA
  x
}

ScidrTimeRevalue <- function(x, declaration, tokens) {
  original <- x
  x <- ScidrTimeMissing(x, tokens)
  seconds <- ScidrTimeSeconds(x, declaration)
  bad <- which(!is.na(seconds) & (!is.finite(seconds) | (declaration$TimeType == "clock" & (seconds < 0 | seconds >= 86400))))
  if (length(bad)) stop("Invalid ", declaration$TimeType, " values at rows: ", paste(bad, collapse = ", "),
    if (declaration$TimeType == "clock") ". Clock values must be in [00:00, 24:00); values were not wrapped." else ". Durations must be finite.")
  if (!inherits(x, "difftime")) {
    if (declaration$TimeType == "clock") {
      if (!requireNamespace("hms", quietly = TRUE)) stop("Install hms to construct clock times.")
      x <- hms::as_hms(seconds)
    } else x <- as.difftime(seconds, units = "secs")
    label <- attr(original, "label", exact = TRUE)
    if (!is.null(label)) attr(x, "label") <- label
  }
  attr(x, "scidr_type") <- declaration$TimeType
  attr(x, "scidr_time_source") <- declaration$Source
  x
}

# Seconds, whatever the column wore. as.numeric() on a difftime returns its
# stored units -- minutes for units = "mins" -- so the unit is named explicitly
# rather than trusted, which would otherwise rescale a column silently. An hms
# is already stored as seconds.
ScidrTimeModelNumeric <- function(x) {
  if (inherits(x, "hms")) return(as.numeric(x))
  if (inherits(x, "difftime")) return(as.numeric(x, units = "secs"))
  x
}

ScidrTimeEncode <- function(df_Data, feature_map) {
  output <- list()
  for (variable in names(feature_map)) {
    features <- feature_map[[variable]]
    # A non-clock temporal column passes through unencoded, but it must not
    # reach a modelling backend still wearing its hms or difftime class: those
    # are rejected as unsupported column types and abort the whole fit. Seconds
    # are the numeric the rest of the pipeline already uses for these columns.
    if (length(features) == 1L && identical(features, variable)) output[[variable]] <-
      ScidrTimeModelNumeric(df_Data[[variable]])
    else {
      angle <- df_Data[[variable]] * (2 * pi / 86400)
      if (any(grepl("_sin$", features))) output[[features[grepl("_sin$", features)]]] <- sin(angle)
      if (any(grepl("_cos$", features))) output[[features[grepl("_cos$", features)]]] <- cos(angle)
    }
  }
  as.data.frame(output, optional = TRUE)
}

ScidrTimeFeatureMap <- function(variables, declarations) {
  clocks <- declarations$Variable[declarations$TimeType == "clock"]
  output <- stats::setNames(lapply(variables, identity), variables)
  occupied <- variables
  for (variable in intersect(variables, clocks)) {
    prefix <- paste0(".scidr_clock_", match(variable, declarations$Variable))
    while (any(paste0(prefix, c("_sin", "_cos")) %in% occupied)) prefix <- paste0(prefix, "_")
    output[[variable]] <- paste0(prefix, c("_sin", "_cos"))
    occupied <- c(occupied, output[[variable]])
  }
  output
}

ScidrTimePredictData <- function(df_Data, adapter = NULL) {
  if (is.null(adapter)) return(df_Data)
  ScidrTimeEncode(df_Data, adapter$FeatureMap)[adapter$ModelFeatures]
}

ScidrTimeGower <- function(df_Data, temporal) {
  clocks <- intersect(names(df_Data), temporal$Variable[temporal$TimeType == "clock"])
  ordinary <- setdiff(names(df_Data), clocks)
  n <- nrow(df_Data)
  total <- matrix(0, n, n)
  if (length(ordinary)) total <- as.matrix(cluster::daisy(df_Data[ordinary], metric = "gower")) * length(ordinary)
  for (variable in clocks) {
    delta <- abs(outer(df_Data[[variable]], df_Data[[variable]], "-"))
    total <- total + pmin(delta, 86400 - delta) / 43200
  }
  stats::as.dist(total / ncol(df_Data))
}

ScidrTimeRegion <- function(component, comparison, values) {
  values <- as.numeric(unlist(values))
  bounds <- values[is.finite(values) & values >= -1 & values <= 1]
  roots <- unlist(lapply(bounds, function(value) {
    angle <- if (component == "sin") c(asin(value), pi - asin(value)) else c(acos(value), 2 * pi - acos(value))
    (angle %% (2 * pi)) * 86400 / (2 * pi)
  }))
  if (comparison %in% c("=", "==", "in", "!=", "<>", "not in")) {
    points <- sort(unique(roots))
    negative <- comparison %in% c("!=", "<>", "not in")
    text <- if (length(points)) paste(ScidrTimeFormat(points, "clock"), collapse = " or ") else "no clock positions"
    return(list(Text = if (negative) paste0("all clock times except ", text) else text,
      Intervals = if (negative) NULL else cbind(Start = points, End = points), ExcludedPoints = if (negative) points else numeric()))
  }
  edges <- sort(unique(c(0, roots, 86400)))
  midpoint <- (head(edges, -1) + tail(edges, -1)) / 2
  coordinate <- if (component == "sin") sin(midpoint * 2 * pi / 86400) else cos(midpoint * 2 * pi / 86400)
  keep <- switch(comparison, `>` = coordinate > values[1], `>=` = coordinate >= values[1],
    `<` = coordinate < values[1], `<=` = coordinate <= values[1],
    between = coordinate >= values[1] & coordinate <= values[2],
    `not between` = coordinate < values[1] | coordinate > values[2],
    `=` = abs(coordinate - values[1]) < 1e-12, `==` = abs(coordinate - values[1]) < 1e-12,
    `!=` = abs(coordinate - values[1]) >= 1e-12, `<>` = abs(coordinate - values[1]) >= 1e-12,
    rep(NA, length(midpoint)))
  if (anyNA(keep)) return(list(Text = "Encoded clock comparison; see native component details", Intervals = NULL))
  intervals <- cbind(Start = head(edges, -1)[keep], End = tail(edges, -1)[keep])
  if (!nrow(intervals)) {
    # Equality constraints can describe isolated clock positions rather than arcs.
    candidates <- sort(unique(roots))
    if (comparison %in% c("=", "==") && length(candidates)) return(list(Text = paste(ScidrTimeFormat(candidates, "clock"), collapse = " or "), Intervals = cbind(Start = candidates, End = candidates)))
    return(list(Text = "No clock interval retained", Intervals = intervals))
  }
  merged <- list(intervals[1, ])
  if (nrow(intervals) > 1L) for (i in 2:nrow(intervals)) {
    last <- length(merged)
    if (abs(merged[[last]][2] - intervals[i, 1]) < 1e-7) merged[[last]][2] <- intervals[i, 2] else merged[[last + 1L]] <- intervals[i, ]
  }
  intervals <- do.call(rbind, merged)
  list(Text = paste(apply(intervals, 1, function(x) paste(ScidrTimeFormat(x, "clock"), collapse = " to ")), collapse = " or "), Intervals = intervals)
}

ScidrTimeFinding <- function(finding, model_variable, feature_map, temporal, original, row) {
  # The backend's conditions and statistics come straight from C++, so they are
  # repaired here, at the one point where a conditional finding enters the
  # package, rather than at each of the places that later render them.
  finding <- ScidrRecoverValues(finding)
  Reverse <- function(feature) {
    if (is.numeric(feature)) feature <- unlist(feature_map, use.names = FALSE)[feature]
    source <- names(feature_map)[vapply(feature_map, function(x) feature %in% x, logical(1))]
    if (length(source) != 1L) stop("Unknown temporal model feature: ", feature)
    source
  }
  variable <- Reverse(model_variable)
  conditions <- lapply(finding$conditions, function(condition) {
    feature <- condition$column
    if (is.numeric(feature)) feature <- unlist(feature_map, use.names = FALSE)[feature]
    source <- Reverse(feature)
    declaration <- temporal[match(source, temporal$Variable), , drop = FALSE]
    condition$column <- source
    if (!is.na(declaration$TimeType) && declaration$TimeType == "clock") {
      region <- ScidrTimeRegion(if (grepl("_sin$", feature)) "sin" else "cos", condition$comparison, condition$value_comp)
      condition$native_component <- list(Feature = feature, Comparison = condition$comparison, Value = condition$value_comp)
      condition$comparison <- "in clock region"
      condition$value_comp <- region$Text
      condition$clock_intervals <- region$Intervals
      seconds <- ScidrTimeSeconds(original[[source]][row], declaration)
      condition$value_this <- ScidrTimeFormat(seconds, "clock")
    }
    condition
  })
  statistics <- finding$group_statistics
  declaration <- temporal[match(variable, temporal$Variable), , drop = FALSE]
  if (!is.na(declaration$TimeType) && declaration$TimeType == "clock") {
    regions <- character()
    component <- if (grepl("_sin$", model_variable)) "sin" else "cos"
    if (length(statistics$lower_thr)) regions <- c(regions, ScidrTimeRegion(component, "<", statistics$lower_thr)$Text)
    if (length(statistics$upper_thr)) regions <- c(regions, ScidrTimeRegion(component, ">", statistics$upper_thr)$Text)
    statistics <- list(n_obs = statistics$n_obs, TimeType = "clock", ClockRegion = paste(regions, collapse = " or "), ComponentStatistics = statistics)
  }
  list(Variable = variable, Conditions = conditions, Statistics = statistics)
}

ScidrTimeObserved <- function(data, variable, rows, temporal = NULL) {
  if (is.null(temporal)) return(ScidrReviewFormat(data[[variable]][rows]))
  declaration <- temporal[match(variable, temporal$Variable), , drop = FALSE]
  if (is.null(temporal) || !nrow(declaration) || is.na(declaration$TimeType)) return(ScidrReviewFormat(data[[variable]][rows]))
  original <- data[[variable]][rows]
  seconds <- tryCatch(ScidrTimeSeconds(original, declaration, strict = FALSE), error = function(e) NULL)
  if (is.null(seconds)) return(ScidrReviewFormat(original))
  output <- ScidrTimeFormat(seconds, declaration$TimeType)
  invalid <- !is.na(seconds) & (!is.finite(seconds) | (declaration$TimeType == "clock" & (seconds < 0 | seconds >= 86400)))
  invalid <- invalid | (!is.na(original) & is.na(seconds))
  output[invalid] <- as.character(original[invalid])
  output
}

ScidrTimeCollapseCells <- function(df_Cells, temporal) {
  if (!nrow(df_Cells) || !any(df_Cells$TemporalType == "clock", na.rm = TRUE)) return(df_Cells)
  pieces <- list(df_Cells[df_Cells$TemporalType != "clock" | is.na(df_Cells$TemporalType), ])
  clock_rows <- which(!is.na(df_Cells$TemporalType) & df_Cells$TemporalType == "clock")
  groups <- split(clock_rows, interaction(df_Cells$.Row[clock_rows], df_Cells$Variable[clock_rows], drop = TRUE))
  for (indices in groups) {
    evidence <- df_Cells[indices, ]
    chosen <- if (all(is.na(evidence$Score))) 1L else which.max(replace(evidence$Score, is.na(evidence$Score), -Inf))
    row <- evidence[chosen, ]
    row$Flag <- if (any(evidence$Flag %in% TRUE)) TRUE else if (any(evidence$Flag %in% FALSE)) FALSE else NA
    row$Expected <- NA_real_
    row$ExpectedDisplay <- NA_character_
    row$NativeComponents <- list(as.data.frame(evidence[setdiff(names(evidence), "NativeComponents")]))
    row$Explanation <- "Clock time flagged on its circular representation; component predictions are technical evidence, not an expected replacement clock time."
    pieces[[length(pieces) + 1L]] <- row
  }
  dplyr::bind_rows(pieces)
}
