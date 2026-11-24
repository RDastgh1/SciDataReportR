ScidrReviewExplanationMeta <- function(object) {
  df_Table <- ScidrReviewContributions(object)
  output <- list()
  for (i in which(object$Evidence$Method == "isolation")) {
    e <- object$Evidence[i, ]
    source <- object$Sources[[e$Analysis]]
    df_Values <- df_Table[df_Table$Analysis == e$Analysis & df_Table$.Row == e$SourceRow, ]
    df_Summary <- source$ExplanationSummary
    status <- if (is.null(df_Summary)) "Column explanations unavailable - rerun detection to generate them" else {
      index <- match(e$SourceRow, df_Summary$.Row)
      if (is.na(index)) "Unavailable" else if (df_Summary$Status[index] == "Unavailable" && "Reason" %in% names(df_Summary)) paste0("Unavailable: ", df_Summary$Reason[index]) else as.character(df_Summary$Status[index])
    }
    if (!nrow(df_Values)) {
      output[[e$IssueID]] <- list(status = status, method = if (is.null(source$ExplanationMetadata$Method)) "shap" else source$ExplanationMetadata$Method, values = list())
      next
    }
    df_Summary <- df_Summary[match(e$SourceRow, df_Summary$.Row), ]
    linked <- unique(unlist(strsplit(df_Values$AdditionalEvidence[nzchar(df_Values$AdditionalEvidence)], "; ", fixed = TRUE)))
    df_Related <- ScidrReviewCategoryFields(object, object$Evidence[object$Evidence$IssueID %in% linked & object$Evidence$Method != "isolation", ])
    output[[e$IssueID]] <- list(status = status, method = if (is.null(source$ExplanationMetadata$Method)) "shap" else source$ExplanationMetadata$Method, score = df_Summary$Score, baseline = df_Summary$Baseline,
      sum = df_Summary$SumContributions, residual = df_Summary$ReconstructionResidual,
      values = lapply(seq_len(nrow(df_Values)), function(j) lapply(as.list(df_Values[j, ]), unname)),
      related = lapply(seq_len(nrow(df_Related)), function(j) as.list(ScidrReviewFlat(df_Related[j, intersect(
        c("IssueID", "Method", "Variable", "Label", "Observed", "Explanation", "Conditions", "ReferenceN", "ReferenceLower", "ReferenceUpper"), names(df_Related)), drop = FALSE]))))
  }
  output
}

ScidrAnomalyContributionsJS <- function() 'function(panel, values, meta) {
  const row = {}; meta.columns.forEach(function(name,i) {row[name]=values[i];});
  const info = (meta.explanations || {})[row.IssueID]; if(!info) return;
  const block=document.createElement("section"), heading=document.createElement("h4");
  const sensitivity=info.method === "sensitivity";
  const effect=function(x) {return sensitivity ? x.ScoreDecrease : x.Contribution;};
  heading.textContent=sensitivity ? "Columns worth checking" : "Columns contributing to this score";block.appendChild(heading);
  const note=document.createElement("p");note.textContent=sensitivity ? "Each column was replaced separately with values from reference records in the same stratum. Positive score decrease means these replacements reduced the isolation score; negative means they increased it. Effects are not additive and do not identify errors or corrections. Reference values may contain errors, and replacements can create unrealistic combinations, especially for correlated columns. Replacement-score ranges are not confidence intervals or error probabilities." : "These approximate contributions explain model behavior, not incorrect values or error probabilities. Positive raises the score relative to the background; negative lowers it. Background records may contain errors. Correlated or derived columns can share attribution. Simulation SD/range describe Monte Carlo sensitivity, not confidence that a value is wrong.";block.appendChild(note);
  if(info.status !== "Explained" || !info.values.length) {
    const state=document.createElement("p");state.textContent=info.status === "Not requested" ? "Not requested: this record was outside the explanation selection or explanations were disabled. Rerun detection with explanation_args$rows to request it." : info.status;
    block.appendChild(state);panel.prepend(block);return;
  }
  const all=info.values.slice(), positive=all.filter(function(x) {return effect(x) > 0;}).sort(function(a,b) {return effect(b)-effect(a) || a.Variable.localeCompare(b.Variable);});
  const negative=all.filter(function(x) {return effect(x) < 0;}).sort(function(a,b) {return effect(a)-effect(b) || a.Variable.localeCompare(b.Variable);});
  const shown=positive.slice(0,5).concat(negative.slice(0,5));
  if(!positive.length) {const message=document.createElement("p");message.textContent="No positive contributors were found; no problem columns are inferred.";block.appendChild(message);}
  const table=function(items) {
    const el=document.createElement("table");el.className="scidr-static-table";
    const head=document.createElement("tr");(sensitivity ? ["Measure / variable","Codebook category","Recorded value","Score decrease","Mean replacement score","Replacement-score range / count"] : ["Measure / variable","Codebook category","Recorded value","Signed contribution","Simulation SD","Simulation range"]).forEach(function(text) {const cell=document.createElement("th");cell.textContent=text;head.appendChild(cell);});
    const thead=document.createElement("thead");thead.appendChild(head);el.appendChild(thead);
    const tbody=document.createElement("tbody");
    const fmt=function(x) {return x == null ? "Unavailable" : Number(x).toPrecision(5);};
    items.forEach(function(x) {const tr=document.createElement("tr");[x.Label+" ("+x.Variable+")",x.CodebookCategory || "Uncategorized",x.Observed == null ? "Missing" : x.Observed,fmt(effect(x)),fmt(sensitivity ? x.ReplacementScoreMean : x.SimulationSD),sensitivity ? fmt(x.ReplacementScoreMin)+" to "+fmt(x.ReplacementScoreMax)+" (n="+x.ReplacementN+")" : fmt(x.SimulationMin)+" to "+fmt(x.SimulationMax)].forEach(function(text) {const cell=document.createElement("td");cell.textContent=text == null ? "" : String(text);tr.appendChild(cell);});tbody.appendChild(tr);});
    el.appendChild(tbody);return el;
  };
  block.appendChild(table(shown));
  const details=document.createElement("details"), title=document.createElement("summary");title.textContent="All "+all.length+(sensitivity ? " column screening results" : " feature contributions");details.appendChild(title);let rendered=false;details.addEventListener("toggle",function() {if(details.open && !rendered) {details.appendChild(table(all));rendered=true;}});block.appendChild(details);
  if(!sensitivity) {
  const total=document.createElement("p");total.textContent="Baseline "+Number(info.baseline).toPrecision(6)+" + summed contributions "+Number(info.sum).toPrecision(6)+" = reconstructed score "+Number(info.baseline+info.sum).toPrecision(6)+"; original score "+Number(info.score).toPrecision(6)+". Numerical residual: "+Number(info.residual).toPrecision(3)+".";block.appendChild(total);
  }
  panel.prepend(block);
  if(typeof Plotly !== "undefined") {
    const chart=document.createElement("div");chart.style.height="400px";chart.style.width="100%";block.appendChild(chart);
    if(sensitivity) {
      const esc=function(text) {return String(text).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");};
      Plotly.newPlot(chart,[{type:"bar",x:shown.map(effect),y:shown.map(function(x) {return esc(x.Label+" ("+x.Variable+")");}),orientation:"h",marker:{color:shown.map(function(x) {return effect(x)>0 ? "#bc3f36" : "#608da8";})},hovertemplate:"%{y}<br>Score decrease: %{x:.5f}<extra></extra>"}],
        {height:400,autosize:true,margin:{t:20,l:180,r:20,b:65},xaxis:{title:"Score decrease after replacing one column",zeroline:true},yaxis:{automargin:true,autorange:"reversed"},showlegend:false},{responsive:true,displaylogo:false,modeBarButtonsToRemove:["toImage","sendDataToCloud"]});
    } else {
    const included=new Set(shown.map(function(x) {return x.Variable;}));
    const remaining=all.filter(function(x) {return !included.has(x.Variable);}).reduce(function(sum,x) {return sum+x.Contribution;},0);
    const esc=function(text) {return String(text).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");};
    const labels=["Background baseline"].concat(shown.map(function(x) {return esc(x.Label+" ("+x.Variable+")");}),["Remaining features","Reconstructed score"]);
    Plotly.newPlot(chart,[{type:"waterfall",orientation:"v",x:labels,y:[info.baseline].concat(shown.map(function(x) {return x.Contribution;}),[remaining,0]),measure:["absolute"].concat(shown.map(function() {return "relative";}),["relative","total"]),increasing:{marker:{color:"#bc3f36"}},decreasing:{marker:{color:"#608da8"}}}],
      {height:400,autosize:true,margin:{t:20,l:65,r:20,b:120},yaxis:{title:"Isolation model score"},xaxis:{automargin:true},showlegend:false},{responsive:true,displaylogo:false,modeBarButtonsToRemove:["toImage","sendDataToCloud"]});
  }
  }
  if(info.related && info.related.length) {
    const additional=document.createElement("details"), label=document.createElement("summary");label.textContent="Additional conditional/date/cell evidence (not confirmation)";additional.appendChild(label);
    info.related.forEach(function(e) {const text=document.createElement("p");text.textContent=e.Method+": "+(e.Label || e.Variable || "")+" - "+(e.Explanation || "")+" [issue "+e.IssueID+"]";additional.appendChild(text);});block.appendChild(additional);
  }
}'
