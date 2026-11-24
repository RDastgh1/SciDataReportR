ScidrReviewPlotSources <- function(object) {
  stats::setNames(lapply(names(object$Sources), function(analysis) {
    source <- ScidrReviewDisplaySource(object, analysis)
    condition_variables <- unlist(lapply(object$Evidence$Conditions, function(conditions)
      vapply(conditions, function(condition) condition$column, character(1))), use.names = FALSE)
    variables <- unique(c(source$Metadata$Variables, object$Evidence$Variable, object$Evidence$OtherVariable, condition_variables, source$Metadata$IDVariable,
      source$Metadata$TimeVariable, source$Metadata$StrataVariables))
    variables <- intersect(variables, names(source$ReviewData))
    df_Data <- source$ReviewData
    declarations <- source$Metadata$TimeConversions
    kinds <- stats::setNames(rep(NA_character_, length(variables)), variables)
    if (isTRUE(nrow(declarations) > 0)) kinds[intersect(variables, declarations$Variable)] <- declarations$TimeType[match(intersect(variables, declarations$Variable), declarations$Variable)]
    list(modelVariables = as.list(intersect(source$Metadata$Variables, names(df_Data))), rows = as.list(df_Data$.Row), id = source$Metadata$IDVariable,
      time = source$Metadata$TimeVariable, strata = as.list(source$Metadata$StrataVariables),
      labels = as.list(source$Metadata$Labels),
      clocks = as.list(stats::setNames(!is.na(kinds) & kinds == "clock", variables)),
      durations = as.list(stats::setNames(!is.na(kinds) & kinds == "duration", variables)),
      numeric = as.list(stats::setNames(vapply(df_Data[variables], is.numeric, logical(1)) & (is.na(kinds) | kinds != "clock") | (!is.na(kinds) & kinds == "duration"), variables)),
      dates = as.list(vapply(df_Data[variables], ScidrAnomalyIsDate, logical(1))),
      values = stats::setNames(lapply(variables, function(variable) {
        x <- df_Data[[variable]]
        if (!is.na(kinds[variable])) {
          declaration <- declarations[match(variable, declarations$Variable), ]
          values <- tryCatch(ScidrTimeSeconds(x, declaration, strict = FALSE), error = function(e) rep(NA_real_, length(x)))
          if (kinds[variable] == "clock") values[!is.na(values) & (!is.finite(values) | values < 0 | values >= 86400)] <- NA_real_
        } else values <- if (ScidrAnomalyIsDate(x)) {
          if (inherits(x, "Date")) as.numeric(x) else as.numeric(x) / 86400
        } else if (is.numeric(x)) x else as.character(x)
        as.list(values)
      }), variables))
  }), names(object$Sources))
}

ScidrReviewPlotJS <- function() paste0('function(panel, values, meta) {
  const clockPlot = ', ScidrAnomalyClockPlotJS(), ';
  if (!meta.plotSources || typeof Plotly === "undefined") return;
  const row = {}; meta.columns.forEach(function(name, i) {row[name] = values[i];});
  const source = meta.plotSources[row.Analysis];
  if (!source) return;
  if (!row.Variable || !source.values[row.Variable]) {
    const index = source.rows.indexOf(Number(row.SourceRow));
    if (index < 0) return;
    const retained = source.rows.map(function(_, i) {return i;}).filter(function(i) {
      return !source.strata || source.strata.every(function(name) {return source.values[name][i] === source.values[name][index];});
    });
    const keys = (source.modelVariables || Object.keys(source.values)).filter(function(name) {
      return name !== source.id && name !== source.time && (!source.strata || source.strata.indexOf(name) < 0) && source.values[name] && new Set(retained.map(function(i) {return source.values[name][i];}).filter(function(v) {return v != null;})).size > 1;
    });
    const heading = document.createElement("h4"); heading.textContent = "Interactive record comparison";
    const note = document.createElement("p"); note.textContent = "Select variables to compare this record with retained records in the same stratum. These views are descriptive, not isolation-forest feature attributions or evidence that these fields are wrong.";

    if (!keys.length) {note.textContent += " No suitable varying model variables are available for plotting."; panel.prepend(note); panel.prepend(heading); return;}
    const controls = document.createElement("div"), chart = document.createElement("div"); chart.style.height = "350px";
    const label = function(name) {return (source.labels[name] || name) + " (" + name + ")" + ((source.durations || {})[name] ? " (seconds)" : "");};
    const esc = function(x) {return String(x).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;");};
    const makeSelect = function(title, blank) {
      const element = document.createElement("select"); element.setAttribute("aria-label",title);
      const add = function(value,text) {const o=document.createElement("option");o.value=value;o.textContent=text;element.appendChild(o);};
      if(blank) add("","Distribution only"); keys.forEach(function(name) {add(name,label(name));});
      const wrapper=document.createElement("label");wrapper.textContent=title+" ";wrapper.appendChild(element);controls.appendChild(wrapper);return element;
    };
    const x = makeSelect("Comparison X",false), y = makeSelect("Comparison Y",true);
    const numeric = keys.filter(function(name) {return source.numeric[name] && !source.dates[name];});
    const explanation=(meta.explanations || {})[row.IssueID];
    const effect=function(v) {return explanation && explanation.method === "sensitivity" ? v.ScoreDecrease : v.Contribution;};
    const contributors=explanation && explanation.status === "Explained" ? explanation.values.filter(function(v) {return effect(v) > 0 && keys.indexOf(v.Variable) >= 0;}).sort(function(a,b) {return effect(b)-effect(a) || a.Variable.localeCompare(b.Variable);}).map(function(v) {return v.Variable;}) : [];
    if(contributors.length) {
      const numericContributors=contributors.filter(function(name) {return numeric.indexOf(name) >= 0;});
      x.value=numericContributors.length >= 2 ? numericContributors[0] : contributors[0];
      y.value=numericContributors.length >= 2 ? numericContributors[1] : "";
      note.textContent += " Initial variables are leading positive score-screening results; they are not identified errors.";
    } else if(explanation && explanation.status === "Explained") {
      note.textContent += " No usable positive contributors were found. Choose comparison variables manually; no problem columns are inferred.";
      x.selectedIndex=-1;y.value="";
    } else {
      x.value=numeric[0] || keys[0];y.value=numeric.length > 1 ? numeric[1] : "";
      note.textContent += " Column explanations are unavailable for this record; initial axes are descriptive defaults.";
    }
    controls.style.display="flex";controls.style.gap="12px";
    const block=document.createElement("section");[heading,note,controls,chart].forEach(function(el) {block.appendChild(el);});panel.prepend(block);
    const baseNote = note.textContent;
    const draw = function() {
      const a=x.value,b=y.value;
      if(!a) return;
      const clock = (source.clocks || {})[a] ? a : ((source.clocks || {})[b] ? b : null);
      if(clock) {note.textContent=baseNote+" Clock selections show a 24-hour distribution, not a joint two-variable plot. Thirty-minute display bins are not detection thresholds.";clockPlot(chart,source,clock,retained,index,esc(label(clock)),350);return;}
      note.textContent = baseNote + (b && !(source.numeric[a] && source.numeric[b] && !source.dates[a] && !source.dates[b]) ? " Two-variable scatterplots require numeric axes; showing the X variable distribution instead." : "");
      if (source.values[a][index] == null || (b && source.values[b][index] == null)) note.textContent += " The selected record has missing values on the chosen axes; a point or highlight may be unavailable.";
      let data, layout={height:350,autosize:true,margin:{t:20,r:20,b:70,l:65},xaxis:{title:esc(label(a)),automargin:true},yaxis:{title:"Observations",automargin:true},showlegend:true};
      if(b && source.numeric[a] && source.numeric[b] && !source.dates[a] && !source.dates[b]) {
        const eligible=retained.filter(function(i) {return source.values[a][i] != null && source.values[b][i] != null;});
        data=[{type:"scatter",mode:"markers",name:"Reference records",x:eligible.map(function(i) {return source.values[a][i];}),y:eligible.map(function(i) {return source.values[b][i];}),marker:{color:"#608da8"}},
          {type:"scatter",mode:"markers",name:"Selected record",x:[source.values[a][index]],y:[source.values[b][index]],marker:{color:"#bc3f36",size:12}}];
        layout.yaxis.title=esc(label(b));
      } else if(source.numeric[a] || source.dates[a]) {
        const value=function(v) {if (!source.dates[a]) return v; if (v == null) return null; const d=new Date(v*86400000); return Number.isFinite(d.getTime()) ? d.toISOString() : null;};
        data=[{type:"histogram",x:retained.map(function(i) {return value(source.values[a][i]);}),name:"Retained values",marker:{color:"#608da8"}}];
        const current=value(source.values[a][index]);if(current != null) layout.shapes=[{type:"line",x0:current,x1:current,yref:"paper",y0:0,y1:1,line:{color:"#bc3f36",width:3}}];
      } else {
        const counts=Object.create(null);retained.forEach(function(i) {const v=source.values[a][i];if(v != null) counts[v]=(counts[v] || 0)+1;});
        const categories=Object.keys(counts);
        data=[{type:"bar",x:categories.map(esc),y:categories.map(function(v) {return counts[v];}),marker:{color:categories.map(function(v) {return v===source.values[a][index] ? "#bc3f36" : "#608da8";})}}];
      }
      layout.showlegend = data.length > 1;
      Plotly.react(chart,data,layout,{responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:["toImage","sendDataToCloud"]});
    };
    x.addEventListener("change",draw);y.addEventListener("change",draw);draw();return;
  }
  const variable = row.Variable, index = source.rows.indexOf(Number(row.SourceRow));
  if (index < 0) return;
  const escape = function(x) {return String(x == null ? "" : x).replace(/[&<>"\x27]/g, function(c) {return {"&":"&amp;","<":"&lt;",">":"&gt;","\\\"":"&quot;","\x27":"&#39;"}[c];});};
  const date = function(x) {if (x == null) return null; const d = new Date(x * 86400000); return Number.isFinite(d.getTime()) ? d.toISOString() : null;};
  const defined = function(x) {return x != null && x !== "";};
  const label = escape((source.labels[variable] || variable) + " (" + variable + ")" + ((source.durations || {})[variable] ? " (seconds)" : ""));
  let data = [], axis = label, shapes = [], annotations = [], note = "All retained observations in this analysis group; this distribution is not restricted to the learned conditional subgroup.";
  let indices = source.rows.map(function(_, i) {return i;}).filter(function(i) {
    return !source.strata || source.strata.every(function(name) {return source.values[name][i] === source.values[name][index];});
  });
  if ((source.clocks || {})[variable]) {
    const heading=document.createElement("h4");heading.textContent="Clock-time distribution";
    const note=document.createElement("p");note.textContent="Circular 24-hour distribution within this stratum. Midnight-adjacent values remain close. Thirty-minute display bins are not detection thresholds or clinical ranges.";
    const chart=document.createElement("div");chart.style.height="350px";panel.prepend(chart);panel.prepend(note);panel.prepend(heading);
    clockPlot(chart,source,variable,indices,index,label,350);return;
  }
  const conditions = (meta.plotConditions || {})[row.IssueID] || [];
  const supported = ["=", "==", "!=", "<>", "<", ">", "<=", ">=", "in", "not in"];
  if (conditions.length && conditions.every(function(c) {return source.values[c.column] && supported.indexOf(c.comparison) >= 0;})) {
    // Conditions arrive as text while the plotted columns keep their own types, so
    // a numeric column compared with "=" or "in" never matched and the panel drew
    // an empty chart. Put both sides in the column\'s own type before comparing.
    const cast = function(column, value) {
      if (value == null) return value;
      if (Array.isArray(value)) return value.map(function(v) {return cast(column, v);});
      return source.numeric[column] && !source.dates[column] ? Number(value) : String(value);
    };
    indices = indices.filter(function(i) {return conditions.every(function(c) {
      const raw = source.values[c.column][i];
      if (raw == null) return false;
      const x = cast(c.column, raw), y = cast(c.column, c.value_comp);
      if (typeof x === "number" && isNaN(x)) return false;
      switch (c.comparison) {
        case "=": case "==": return x === y;
        case "!=": case "<>": return x !== y;
        case "<": return x < y;
        case ">": return x > y;
        case "<=": return x <= y;
        case ">=": return x >= y;
        case "in": return (Array.isArray(y) ? y : [y]).indexOf(x) >= 0;
        case "not in": return (Array.isArray(y) ? y : [y]).indexOf(x) < 0;
      }
    });});
    note = indices.length
      ? indices.length + " record" + (indices.length === 1 ? "" : "s") + " here meet the learned conditions, including the highlighted record. This count can differ from the reference count in the finding, because exclusions applied during detection are not reproduced in this panel."
      : "No record in this view meets the learned conditions, so there is nothing to plot against. This happens when the conditions refer to values that detection saw but this panel does not; use the reference count reported with the finding instead.";
  }
  const current = source.values[variable][index], other = source.values[row.OtherVariable];
  if (String(row.Methods || row.Method).split(', ').indexOf("missingness") >= 0) {
    const entries=indices.map(function(i) {return source.values[variable][i];});
    const absent=entries.filter(function(v) {return v == null;}).length;
    data=[{type:"bar",x:["Observed","Missing"],y:[entries.length-absent,absent],marker:{color:["#608da8","#bc3f36"]},name:"Retained records"}];
    axis="Completion status";
    note="Observed versus missing values within this stratum. The selected finding is a blank, not a numeric outlier; learned reference counts are shown separately.";
  } else if (source.time && source.id && source.values[source.time] && source.numeric[variable]) {
    const id = source.values[source.id][index];
    const history_indices = source.rows.map(function(_, i) {return i;}).filter(function(i) {
      return source.values[source.id][i] === id && (!source.strata || source.strata.every(function(name) {return source.values[name][i] === source.values[name][index];}));
    }).sort(function(a,b) {return source.values[source.time][a] - source.values[source.time][b];});
    const time = function(i) {const x = source.values[source.time][i]; return source.dates[source.time] ? date(x) : x;};
    const selected = history_indices.filter(function(i) {return source.rows[i] === Number(row.SourceRow) || source.rows[i] === Number(row.OtherRow);});
    data = [{type:"scatter",mode:"lines+markers",name:"Retained history",x:history_indices.map(time),y:history_indices.map(function(i) {return source.values[variable][i];}),connectgaps:false},
      {type:"scatter",mode:"markers",name:"Finding endpoints",x:selected.map(time),y:selected.map(function(i) {return source.values[variable][i];}),marker:{color:"#bc3f36",size:11}}];
    axis = escape((source.labels[source.time] || source.time) + " (" + source.time + ")");
    note = "Complete retained participant history within this analysis group; highlighted points identify the finding, not which endpoint is incorrect.";
  } else if (other && source.dates[variable] && source.dates[row.OtherVariable]) {
    const intervals = indices.map(function(i) {const x = source.values[variable][i]; return x == null || other[i] == null ? null : other[i] - x;}).filter(defined);
    const actual = current == null || other[index] == null ? null : other[index] - current;
    data = [{type:"histogram",x:intervals,name:"Retained intervals",marker:{color:"#608da8"}}];
    if (actual != null) shapes = [{type:"line",x0:actual,x1:actual,yref:"paper",y0:0,y1:1,line:{color:"#bc3f36",width:3}}];
    axis = escape((source.labels[row.OtherVariable] || row.OtherVariable) + " (" + row.OtherVariable + ") minus " + (source.labels[variable] || variable) + " (" + variable + ") (days)");
    note = "Signed date intervals across all retained complete pairs; red line is the selected interval. Negative intervals are not automatically invalid.";
  } else if (source.numeric[variable] || source.dates[variable]) {
    const xs = indices.map(function(i) {return source.values[variable][i];}).filter(defined).map(function(x) {return source.dates[variable] ? date(x) : x;});
    const actual = source.dates[variable] ? date(current) : current;
    data = [{type:"histogram",x:xs,name:"Retained values",marker:{color:"#608da8"}}];
    if (actual != null) shapes = [{type:"line",x0:actual,x1:actual,yref:"paper",y0:0,y1:1,line:{color:"#bc3f36",width:3}}];
    // The boundary the finding is stated against, so "above the subgroup boundary
    // of 285" is visible rather than a number the reader has to picture.
    const bounds = (meta.plotBounds || {})[row.IssueID];
    if (bounds && !source.dates[variable]) ["lower","upper"].forEach(function(side) {
      if (bounds[side] == null) return;
      shapes.push({type:"line",x0:bounds[side],x1:bounds[side],yref:"paper",y0:0,y1:1,
        line:{color:"#6b7a86",width:2,dash:"dash"}});
      annotations.push({x:bounds[side],yref:"paper",y:1,yanchor:"bottom",showarrow:false,
        font:{size:11,color:"#6b7a86"},text:"learned " + side + " boundary"});
    });
  } else {
    const counts = Object.create(null); indices.map(function(i) {return source.values[variable][i];}).filter(defined).forEach(function(x) {counts[x] = (counts[x] || 0) + 1;});
    const categories = Object.keys(counts);
    data = [{type:"bar",x:categories.map(escape),y:categories.map(function(x) {return counts[x];}),marker:{color:categories.map(function(x) {return x === current ? "#bc3f36" : "#608da8";})}}];
  }
  const heading = document.createElement("h4"); heading.textContent = "Interactive view";
  const caption = document.createElement("p"); caption.textContent = note + " Hover for values; drag to zoom; double-click to reset.";
  const chart = document.createElement("div"); chart.style.height = "330px"; chart.style.width = "100%";
  panel.prepend(chart); panel.prepend(caption); panel.prepend(heading);
  Plotly.newPlot(chart,data,{height:330,autosize:true,margin:{t:20,r:20,b:70,l:65},xaxis:{title:axis,automargin:true},yaxis:{title:data.length > 1 ? label : "Observations",automargin:true},shapes:shapes,annotations:annotations,showlegend:data.length > 1},
    {responsive:true,displaylogo:false,scrollZoom:true,modeBarButtonsToRemove:["toImage","sendDataToCloud"]});
}')

ScidrAnomalyClockPlotJS <- function() 'function(chart,source,variable,indices,index,label,height) {
  const counts=Array(48).fill(0);
  indices.forEach(function(i) {const value=source.values[variable][i];if(value != null && value>=0 && value<86400) counts[Math.floor(value/1800)]++;});
  const text=counts.map(function(_,i) {const h=Math.floor(i/2),m=i%2 ? "30" : "00";return String(h).padStart(2,"0")+":"+m;});
  const data=[{type:"barpolar",theta:counts.map(function(_,i) {return (i+0.5)*7.5;}),r:counts,width:7.5,text:text,name:"Reference observations",marker:{color:"#608da8"},hovertemplate:"%{text}<br>30-minute bin: %{r} observations<extra></extra>"}];
  const current=source.values[variable][index];
  if(current != null && current>=0 && current<86400) data.push({type:"scatterpolar",theta:[current/240],r:[Math.max(...counts,1)+1],mode:"markers",name:"Selected record",marker:{color:"#bc3f36",size:12},hovertemplate:"Selected clock time<extra></extra>"});
  Plotly.newPlot(chart,data,{height:height,autosize:true,margin:{t:25,r:55,b:35,l:55},polar:{angularaxis:{rotation:90,direction:"clockwise",tickvals:[0,90,180,270],ticktext:["00:00","06:00","12:00","18:00"]},radialaxis:{title:"Count",rangemode:"tozero"}},annotations:[{text:label,x:0.5,y:-0.10,xref:"paper",yref:"paper",showarrow:false}],showlegend:true},
    {responsive:true,displaylogo:false,modeBarButtonsToRemove:["toImage","sendDataToCloud"]});
}'
