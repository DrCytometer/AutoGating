convert_gates_to_table <- function(gate_list) {
  gate_table <- data.frame(
    gate.number = integer(),
    gate.name = character(),
    popul.name = character(),
    gate.marker.x = character(),
    gate.marker.y = character(),
    parent.gate = character(),
    parent.popul = character(),
    popul.label = character(),
    popul.label.pos = character(),
    gate.type = character(),
    gate.algorithm = character(),
    gate.param = character(),
    density.quantile.x = numeric(),
    density.quantile.y = numeric(),
    stringsAsFactors = FALSE
  )

  gate_number <- 0
  processed_gates <- list() # Track unique gates

  for (gate_name in names(gate_list)) {
    gate_info <- gate_list[[gate_name]]
    gate_type <- if (inherits(gate_info, "polygonGate")) "polygonal" else "rectangular"
    dimensions <- parameters(gate_info)

    if (length(dimensions) < 1) next # Skip empty gates

    gate_marker_x <- ifelse(length(dimensions) > 1, dimensions[2], NA)
    gate_marker_y <- dimensions[1]

    parent_gate <- sub("/[^/]+$", "", gate_name)  # Get parent gate path
    parent_popul <- sub(".*/", "", parent_gate)   # Extract population name from parent gate

    gate_name_clean <- gsub("^\\$`|`$", "", gate_name)
    gate_name_clean <- sub(".*/", "", gate_name_clean)
    gate_name_clean <- gsub(" ", ".", gate_name_clean)

    popul_label <- gate_info@filterId

    # Define position based on label pattern
    popul_label_pos <- if (grepl("\\*", popul_label)) {
      "*"
    } else if (length(strsplit(popul_label, ";")[[1]]) == 1) {
      "4"
    } else {
      "2;4"
    }

    gate_algorithm <- if (gate_type == "free" && grepl("All", gate_name)) {
      "all"
    } else if (gate_type == "free" && grepl("singlets", gate_name)) {
      "singlets"
    } else if (gate_type == "1dsep") {
      "bimodal"
    } else {
      "free"
    }

    # Create a unique key for gates that share the same parameter and parent
    gate_key <- paste(parent_gate, gate_marker_y, sep = "_")

    if (!gate_key %in% names(processed_gates)) {
      # Create new gate entry
      gate_number <- gate_number + 1
      processed_gates[[gate_key]] <- gate_number

      gate_table <- rbind(gate_table, data.frame(
        gate.number = gate_number,
        gate.name = gate_name_clean,
        popul.name = gate_name_clean,
        gate.marker.x = gate_marker_x,
        gate.marker.y = gate_marker_y,
        parent.gate = gsub("^\\$`|`$", "", parent_gate),
        parent.popul = gsub("^\\$`|`$", "", parent_popul),
        popul.label = popul_label,
        popul.label.pos = popul_label_pos,
        gate.type = gate_type,
        gate.algorithm = gate_algorithm,
        gate.param = paste0("quantile.x=0;quantile.y=0"),
        density.quantile.x = 0,
        density.quantile.y = 0,
        stringsAsFactors = FALSE
      ))

    } else {
      # Append population name to existing gate entry
      existing_gate_number <- processed_gates[[gate_key]]
      gate_table[gate_table$gate.number == existing_gate_number, "popul.name"] <- paste(
        gate_table[gate_table$gate.number == existing_gate_number, "popul.name"],
        gate_name_clean,
        sep = ";"
      )
    }
  }

  return(gate_table)
}

# Example usage
gate_table <- convert_gates_to_table(gate.list)
print(gate_table)
