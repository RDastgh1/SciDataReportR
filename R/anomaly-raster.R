# Anomaly map: one cell per variable/record, coloured by the kind of anomaly
# carried by that cell. The categories are deliberately coarser than the anomaly
# type: a reader scanning a wall of cells can distinguish blank, too high, too
# low, conditional and longitudinal at a glance, which ten types do not allow.
ScidrAnomalyRasterCategories <- function() data.frame(
  Key = c("missing", "high", "low", "conditional", "change", "other"),
  Category = c("Blank where comparable records have a value", "Higher than expected",
    "Lower than expected", "Unusual given other fields", "Unusual change over time",
    "Other flagged value"),
  Colour = c("#e08214", "#b2182b", "#2166ac", "#762a83", "#1b7837", "#8c7a00"),
  stringsAsFactors = FALSE)

# A cell takes the category of its most serious finding, so a variable carrying
# both a directional and a generic finding is drawn as the directional one.
ScidrAnomalyRasterCategoryKey <- function(type, observed, lower, upper) {
  if (identical(type, "Unexpected missingness")) return("missing")
  value <- suppressWarnings(as.numeric(observed))
  lower <- suppressWarnings(as.numeric(lower))
  upper <- suppressWarnings(as.numeric(upper))
  if (length(value) == 1L && is.finite(value)) {
    if (length(upper) == 1L && is.finite(upper) && value > upper) return("high")
    if (length(lower) == 1L && is.finite(lower) && value < lower) return("low")
  }
  if (identical(type, "Conditional outlier")) return("conditional")
  if (type %in% c("Unexpected change over time", "Uncommon transition")) return("change")
  "other"
}

# The payload is sparse: variables, records and the flagged cells only. The grey
# background matrix is built in the browser, so an all-grey frame of 100,000
# cells costs nothing in the saved file.
ScidrAnomalyRasterData <- function(object, display, analysis) {
  source <- ScidrReviewDisplaySource(object, analysis)
  df_Data <- source$ReviewData
  variables <- intersect(source$Metadata$Variables, names(df_Data))
  catalog <- source$Metadata$VariableMetadata
  labels <- if (is.null(catalog)) variables else {
    index <- match(variables, catalog$Variable)
    ifelse(is.na(index) | is.na(catalog$Label[index]), variables, catalog$Label[index])
  }
  rows <- df_Data$.Row
  Context <- function(variable) {
    if (!length(variable) || is.na(variable) || !variable %in% names(df_Data)) return(NULL)
    as.list(ScidrReviewFormat(df_Data[[variable]]))
  }
  df_Evidence <- display$Evidence[display$Evidence$Analysis == analysis, , drop = FALSE]
  placed <- !is.na(df_Evidence$Variable) & df_Evidence$Variable %in% variables &
    !is.na(df_Evidence$SourceRow) & df_Evidence$SourceRow %in% rows
  unplaced <- length(unique(df_Evidence$IssueID[!placed]))
  df_Evidence <- df_Evidence[placed, , drop = FALSE]
  if (nrow(df_Evidence) && "ReviewTier" %in% names(df_Evidence))
    df_Evidence <- df_Evidence[order(df_Evidence$ReviewTier, df_Evidence$EvidenceID), , drop = FALSE]
  df_Categories <- ScidrAnomalyRasterCategories()
  cells <- list()
  seen <- character()
  for (i in seq_len(nrow(df_Evidence))) {
    e <- df_Evidence[i, ]
    key <- paste(e$SourceRow, e$Variable, sep = "\r")
    if (key %in% seen) next
    seen <- c(seen, key)
    category <- ScidrAnomalyRasterCategoryKey(e$FindingType, e$Observed, e$ReferenceLower, e$ReferenceUpper)
    value <- if (identical(category, "missing")) "(blank)" else ScidrReviewRoundedValue(e$Observed)
    cells[[length(cells) + 1L]] <- list(
      row = match(e$SourceRow, rows) - 1L, variable = match(e$Variable, variables) - 1L,
      category = match(category, df_Categories$Key),
      type = e$FindingType, value = value, reason = e$WhyFlagged)
  }
  list(variables = as.list(variables), labels = as.list(unname(labels)), rows = as.list(rows),
    id = source$Metadata$IDVariable, idValues = Context(source$Metadata$IDVariable),
    time = source$Metadata$TimeVariable, timeValues = Context(source$Metadata$TimeVariable),
    categories = as.list(df_Categories$Category), colours = as.list(df_Categories$Colour),
    cells = cells, unplaced = unplaced)
}

ScidrAnomalyRasterJS <- function() paste0('function(element, map) {
  if (typeof Plotly === "undefined") return;
  const rows = map.rows.length, variables = map.variables.length;
  if (!rows || !variables) return;
  const z = [], text = [];
  for (let v = 0; v < variables; v++) {
    z.push(new Array(rows).fill(0));
    text.push(new Array(rows).fill(""));
  }
  const esc = function(x) {return String(x == null ? "" : x).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");};
  const record = function(i) {
    const parts = ["Source row " + map.rows[i]];
    if (map.id && map.idValues) parts.push(esc(map.id) + ": " + esc(map.idValues[i]));
    if (map.time && map.timeValues) parts.push(esc(map.time) + ": " + esc(map.timeValues[i]));
    return parts.join("<br>");
  };
  const wrap = function(x, width) {
    const words = String(x == null ? "" : x).split(" "), lines = [];
    let line = "";
    words.forEach(function(word) {
      if ((line + " " + word).trim().length > width) {lines.push(line.trim()); line = word;}
      else line = line + " " + word;
    });
    if (line.trim().length) lines.push(line.trim());
    return lines.join("<br>");
  };
  map.cells.forEach(function(cell) {
    z[cell.variable][cell.row] = cell.category;
    text[cell.variable][cell.row] = record(cell.row) + "<br><b>" + esc(map.labels[cell.variable]) +
      "</b> (" + esc(map.variables[cell.variable]) + ") = " + esc(cell.value) +
      "<br><b>" + esc(map.categories[cell.category - 1]) + "</b><br>" +
      wrap(esc(cell.type) + ": " + esc(cell.reason), 70);
  });
  for (let v = 0; v < variables; v++) for (let i = 0; i < rows; i++)
    if (!text[v][i]) text[v][i] = record(i) + "<br><b>" + esc(map.labels[v]) + "</b> (" +
      esc(map.variables[v]) + ")<br>No anomaly reported for this value.";
  const levels = map.categories.length;
  const palette = ["#e6ebef"].concat(map.colours);
  const colorscale = [];
  for (let k = 0; k <= levels; k++) {
    colorscale.push([k / (levels + 1), palette[k]]);
    colorscale.push([(k + 1) / (levels + 1), palette[k]]);
  }
  const tickvals = [], ticktext = [];
  for (let k = 0; k <= levels; k++) {
    tickvals.push(k);
    ticktext.push(k === 0 ? "No anomaly reported" : map.categories[k - 1]);
  }
  Plotly.newPlot(element, [{
    type: "heatmap", z: z, text: text, hoverinfo: "text", xgap: 0, ygap: 0,
    zmin: -0.5, zmax: levels + 0.5, colorscale: colorscale,
    colorbar: {tickmode: "array", tickvals: tickvals, ticktext: ticktext, thickness: 14, len: 1, outlinewidth: 0}
  }], {
    margin: {l: variables <= 60 ? 190 : 12, r: 12, t: 10, b: 34},
    height: Math.max(260, Math.min(28 * variables + 70, 900)),
    xaxis: {title: "Records, in saved source-row order", showticklabels: false, ticks: "", showgrid: false, zeroline: false},
    yaxis: {autorange: "reversed", showgrid: false, zeroline: false, ticks: "",
      showticklabels: variables <= 60, tickmode: "array",
      tickvals: map.variables.map(function(_, i) {return i;}),
      ticktext: map.labels.map(function(label) {return String(label).length > 28 ? String(label).slice(0, 27) + "\\u2026" : String(label);})},
    hoverlabel: {align: "left"}
  }, {displayModeBar: false, responsive: true});
}')

ScidrAnomalyRasterTag <- function(object, display, analysis, multiple) {
  map <- ScidrAnomalyRasterData(object, display, analysis)
  if (!length(map$variables) || !length(map$rows)) return(NULL)
  id <- paste0("scidr-raster-", substr(digest::digest(list(analysis, length(map$cells), length(map$rows))), 1, 10))
  payload <- jsonlite::toJSON(map, auto_unbox = TRUE, null = "null", na = "null", digits = NA)
  htmltools::tagList(
    if (multiple) htmltools::tags$h4(analysis),
    htmltools::tags$div(id = id, class = "scidr-raster"),
    htmltools::tags$p(class = "scidr-raster-note",
      paste0(length(map$cells), " of ", length(map$rows) * length(map$variables),
        " values carry a reported anomaly across ", length(map$variables), " variables and ",
        length(map$rows), " records.",
        if (map$unplaced > 0) paste0(" ", map$unplaced,
          " further finding", if (map$unplaced > 1L) "s" else "",
          " concern a whole record or a column outside the modelled set and cannot be placed on a cell.") else "",
        if (length(map$variables) > 60L) " Variable names are omitted above 60 variables; hover a cell for its name." else "")),
    htmltools::tags$script(htmltools::HTML(paste0(
      "(function(){const draw = ", ScidrAnomalyRasterJS(), ";\n",
      "const map = ", payload, ";\n",
      "let attempts = 0;\n",
      "const start = function(){\n",
      "  const element = document.getElementById('", id, "');\n",
      "  if (!element) return;\n",
      "  if (typeof Plotly === 'undefined') { if (attempts++ < 100) setTimeout(start, 100); return; }\n",
      "  draw(element, map);\n};\n",
      "if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();})();"))))
}
