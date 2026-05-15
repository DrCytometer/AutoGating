# ==============================================================================
# autogating_gui/app.R
#
# Interactive Shiny GUI for building AutoGating gate definitions.
#
# Tabs
#   1. Setup          — load FCS files, select cytometer, load existing gates
#   2. Transformations — per-channel biexp tuning with ridgeline plots
#   3. Gate Editor    — interactive gate definition (population, scatter/histogram)
#   4. Gate Tree      — view / reorder / delete defined gates
#   5. Export         — save gate definitions as R source or RDS
#
# Dependencies:
#   shiny, bslib, shinyjs, flowCore, ggplot2, dplyr, DT
# ==============================================================================

library(shiny)
library(bslib)
library(shinyjs)
library(flowCore)
library(ggplot2)
library(dplyr)
library(DT)

# Allow FCS files up to 2 GB
options(shiny.maxRequestSize = 2 * 1024^3)

# Suppress DT/jsonlite named-vector deprecation warning
options(DT.TOJSON_ARGS = list(keep_vec_names = FALSE))

# ---- source AutoGating package helpers (adjust path as needed) ---------------
# source("../get_autogating_param.r")
# source("../get_autogating_param_aurora.r")
# source("../get_autogating_param_id7000.r")
# source("../get_autogating_param_symphony.r")
# source("../analyze_aux.r")
# source("../calculate_threshold.R")
# source("../calculate_threshold_kde_min.R")
# source("../calculate_threshold_mixture.R")
# source("../calculate_threshold_otsu.R")
# source("../resolve_gate_param.R")
# source("../calculate_gate.r")
# source("../apply_gate.R")

# ==============================================================================
# Constants
# ==============================================================================

CYTOMETERS <- c(
  "BD LSRFortessa"  = "bd_fortessa",
  "Cytek Aurora"    = "aurora",
  "Sony ID7000"     = "id7000",
  "BD FACSymphony"  = "symphony"
)

GATE_TYPES <- c(
  "Free / all-events"      = "free",
  "Singlets"               = "singlets",
  "1-D separator (1dsep)"  = "1dsep",
  "2-D separator (2dsep)"  = "2dsep",
  "Replicate from gate"    = "replicate"
)

ALGORITHMS_1D <- c(
  "Tail probability" = "tail",
  "Bimodal KDE"      = "bimodal",
  "KDE minimum"      = "kde.min",
  "Otsu"             = "otsu",
  "Mixture model"    = "mixture"
)

ALGORITHMS_2D <- c(
  "Bimodal / Bimodal" = "bimodal:bimodal",
  "Tail / Tail"       = "tail:tail",
  "Tail / Bimodal"    = "tail:bimodal",
  "Bimodal / Tail"    = "bimodal:tail"
)

MAX_XFORM_EVENTS   <- 30000L   # events shown in transformation histogram
MAX_SCATTER_EVENTS <- 50000L   # events shown in gate editor scatter
IGNORE_TOKEN      <- "*"

# ==============================================================================
# UI
# ==============================================================================

ui <- page_navbar(
  title     = "AutoGating Gate Editor",
  theme     = bs_theme(bootswatch = "flatly", primary = "#2C6FAC"),
  header    = useShinyjs(),

  # --------------------------------------------------------------------------
  # TAB 1 — Setup
  # --------------------------------------------------------------------------
  nav_panel(
    "Setup",
    icon = icon("gear"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320,
        h5("1. Load FCS Files"),
        fileInput("fcs_files", "Select FCS file(s)",
                  multiple = TRUE, accept = ".fcs"),
        h5("2. Load Panel CSV"),
        fileInput("panel_csv", "Panel CSV (dye / antigen / type)",
                  accept = ".csv"),
        h5("3. Cytometer"),
        selectInput("cytometer", "Cytometer",
                    choices = CYTOMETERS, selected = "bd_fortessa"),
        hr(),
        h5("4. Load Existing Gates (optional)"),
        fileInput("pop_gates_file",  "Population gates (.r / .rds)"),
        fileInput("act_gates_file",  "Activation gates (.r / .rds)"),
        hr(),
        h5("5. Compensation CSV (optional)"),
        fileInput("comp_csv", "Compensation CSV", accept = ".csv"),
        hr(),
        actionButton("btn_load", "Load & Initialise",
                     class = "btn-primary btn-block",
                     icon  = icon("play"))
      ),
      card(
        card_header("Setup Status"),
        verbatimTextOutput("setup_status"),
        card_body(
          h6("Panel"),
          DTOutput("panel_table")
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # TAB 2 — Transformations
  # --------------------------------------------------------------------------
  nav_panel(
    "Transformations",
    icon = icon("sliders"),
    layout_sidebar(
      sidebar = sidebar(
        width = 300,
        selectInput("ridge_population", "Display population",
                    choices  = c("All events (ungated)" = "__all__"),
                    selected = "__all__"),
        hr(),
        h6("Channel"),
        uiOutput("transform_channel_ui"),
        fluidRow(
          column(6, actionButton("btn_prev_ch", "← Prev",
                                 class = "btn-outline-secondary btn-sm w-100")),
          column(6, actionButton("btn_next_ch", "Next →",
                                 class = "btn-outline-secondary btn-sm w-100"))
        ),
        hr(),
        h6("Width Basis"),
        sliderInput("xform_width", NULL,
                    min = -5000, max = 0, value = -100, step = 10),
        numericInput("xform_width_num", NULL, value = -100),
        h6("Pos"),
        numericInput("xform_pos", "Pos",  value = 4.418, step = 0.1),
        h6("Neg"),
        numericInput("xform_neg", "Neg",  value = 0,     step = 0.1),
        hr(),
        actionButton("btn_apply_xform",  "Apply to channel",
                     class = "btn-warning btn-sm"),
        actionButton("btn_reset_xform",  "Reset to defaults",
                     class = "btn-outline-secondary btn-sm"),
        hr(),
        downloadButton("btn_save_xform", "Save transformation CSV",
                       class = "btn-success btn-sm")
      ),
      card(
        card_header(uiOutput("xform_plot_header")),
        plotOutput("xform_plot", height = "480px")
      )
    )
  ),

  # --------------------------------------------------------------------------
  # TAB 3 — Gate Editor
  # --------------------------------------------------------------------------
  nav_panel(
    "Gate Editor",
    icon = icon("crop"),
    layout_columns(
      col_widths = c(3, 5, 4),

      # ---- left panel: gate definition fields --------------------------------
      card(
        card_header("Gate Definition"),
        card_body(
          textInput("gate_name", "Gate name", placeholder = "e.g. leukocytes"),
          selectInput("gate_type", "Gate type", choices = GATE_TYPES, selected = "free"),
          selectInput("parent_gate", "Parent gate",
                      choices = c("(none — root)" = "__root__")),
          selectInput("parent_popul", "Parent population",
                      choices = c()),
          hr(),
          selectInput("gate_marker_x", "X-axis marker", choices = c()),
          selectInput("gate_marker_y", "Y-axis marker", choices = c()),
          hr(),
          uiOutput("algorithm_ui"),
          hr(),
          h6("Sub-populations"),
          p(class = "text-muted small",
            "For 1dsep: neg then pos; 2dsep: DN, Xpos, Ypos, DP.",
            "Use * to ignore a slot."),
          uiOutput("popul_fields_ui"),
          hr(),
          h6("Reference / negative control (optional)"),
          selectInput("ref_gate",  "Reference gate",  choices = c("(none)" = "")),
          selectInput("ref_popul", "Reference population", choices = c()),
          hr(),
          h6("Search bounds (fraction of range)"),
          sliderInput("region_min", "Lower search bound",
                      min = 0, max = 0.5, value = 0.01, step = 0.01),
          sliderInput("region_max", "Upper search bound",
                      min = 0.5, max = 1,  value = 0.99, step = 0.01),
          hr(),
          h6("Stats parent"),
          checkboxInput("stats_parent_immediate", "Use immediate parent", value = TRUE),
          conditionalPanel(
            condition = "!input.stats_parent_immediate",
            selectInput("stats_parent_gate",  "Stats gate",   choices = c()),
            selectInput("stats_parent_popul", "Stats popul",  choices = c())
          ),
          hr(),
          actionButton("btn_calculate_gate", "Calculate gate",
                       class = "btn-primary",
                       icon  = icon("calculator")),
          actionButton("btn_add_gate", "Add gate to definition",
                       class = "btn-success",
                       icon  = icon("plus")),
          br(), br(),
          verbatimTextOutput("gate_calc_status")
        )
      ),

      # ---- centre: scatter / histogram plot ----------------------------------
      card(
        card_header(
          "Parent population scatter / histogram",
          actionButton("btn_refresh_plot", "", icon = icon("sync"),
                       class = "btn-outline-secondary btn-sm float-end")
        ),
        card_body(
          plotOutput("gate_plot",
                     click    = "gate_plot_click",
                     dblclick = "gate_plot_dblclick",
                     height   = "440px"),
          p(class = "text-muted small mt-1",
            "Click to set the initial gate position hint. ",
            "Click again to reposition. Double-click to clear."),
          fluidRow(
            column(6, verbatimTextOutput("click_info")),
            column(6,
                   sliderInput("max_scatter_pts", "Max display points",
                               min = 5000, max = 100000, value = 50000,
                               step = 5000))
          )
        )
      ),

      # ---- right: live gate parameter review ---------------------------------
      card(
        card_header("Current gate parameters (preview)"),
        card_body(
          verbatimTextOutput("gate_param_preview")
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # TAB 4 — Gate Tree
  # --------------------------------------------------------------------------
  nav_panel(
    "Gate Tree",
    icon = icon("sitemap"),
    layout_columns(
      col_widths = c(7, 5),
      card(
        card_header("Defined Gates"),
        card_body(
          DTOutput("gate_tree_table"),
          br(),
          actionButton("btn_delete_gate", "Delete selected gate",
                       class = "btn-danger btn-sm",
                       icon  = icon("trash"))
        )
      ),
      card(
        card_header("Gate detail (read-only R list)"),
        card_body(
          verbatimTextOutput("gate_detail_preview")
        )
      )
    )
  ),

  # --------------------------------------------------------------------------
  # TAB 5 — Export
  # --------------------------------------------------------------------------
  nav_panel(
    "Export",
    icon = icon("download"),
    card(
      card_header("Export Gate Definitions"),
      card_body(
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("Population Gates"),
            selectInput("export_pop_format", "Format",
                        choices = c("R source file (.r)" = "r",
                                    "RDS file (.rds)"    = "rds")),
            textInput("export_pop_varname", "Variable name",
                      value = "fcs.population.gates"),
            downloadButton("btn_export_pop", "Download population gates",
                           class = "btn-primary")
          ),
          card(
            card_header("Activation Gates"),
            selectInput("export_act_format", "Format",
                        choices = c("R source file (.r)" = "r",
                                    "RDS file (.rds)"    = "rds")),
            textInput("export_act_varname", "Variable name",
                      value = "fcs.activation.gates"),
            downloadButton("btn_export_act", "Download activation gates",
                           class = "btn-primary")
          )
        ),
        hr(),
        h6("Population gate definitions (R source preview)"),
        verbatimTextOutput("export_pop_preview"),
        h6("Transformation CSV preview"),
        DTOutput("export_xform_table"),
        downloadButton("btn_export_xform_csv",
                       "Download transformation CSV",
                       class = "btn-success")
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================

server <- function(input, output, session) {

  # --------------------------------------------------------------------------
  # Reactive state
  # --------------------------------------------------------------------------

  rv <- reactiveValues(
    # raw data
    fcs_list          = list(),   # named list of flowFrame
    fcs_expr          = NULL,     # combined matrix, all channels
    panel             = NULL,     # data.frame: dye / antigen / type
    ap                = NULL,     # autogating param list
    initialized       = FALSE,

    # transformations: named list, one row per fluorescence channel
    xform_params      = NULL,     # data.frame: dye, pos, neg, width, length, max.range
    xform_cache       = list(),   # named by channel: transformed expression vector

    # gate definitions
    population_gates  = list(),
    activation_gates  = list(),

    # gate calculation results for current gate-editor gate
    fcs_gates_temp    = list(),   # running index of idx per gate/popul
    current_boundary  = NULL,
    calc_gate_type    = NULL,     # gate type used for the most recent boundary calculation
    click_hint        = NULL,     # data.frame(x,y) from user click

    gate_number_counter = 1L,
    comp_applied        = FALSE
  )

  # --------------------------------------------------------------------------
  # Helpers
  # --------------------------------------------------------------------------

  get_scatter_params <- function() {
    req(rv$panel)
    rv$panel$antigen[ rv$panel$type == "scatter" ]
  }

  get_fluor_params <- function() {
    req(rv$panel)
    rv$panel$antigen[ rv$panel$type %in% c("lineage", "activation") ]
  }

  get_all_params <- function() {
    req(rv$panel)
    rv$panel$antigen
  }

  # Return transformed expression for a channel (cached).
  get_transformed <- function(channel) {
    if (is.null(rv$fcs_expr) || is.null(rv$xform_params) || is.null(rv$panel))
      stop("Data not yet loaded.")
    if (!is.null(rv$xform_cache[[channel]])) return(rv$xform_cache[[channel]])
    if (channel %in% get_scatter_params()) {
      vals <- rv$fcs_expr[, channel]
    } else {
      row <- rv$xform_params[ rv$xform_params$antigen == channel, ]
      if (nrow(row) == 0) return(rv$fcs_expr[, channel])
      xf  <- flowjo_biexp(
        channelRange = row$length[1],
        maxValue     = row$max.range[1],
        pos          = row$pos[1],
        neg          = row$neg[1],
        widthBasis   = row$width[1],
        inverse      = FALSE
      )
      vals <- xf(rv$fcs_expr[, channel])
    }
    rv$xform_cache[[channel]] <- vals
    vals
  }

  # Invalidate all xform_cache entries (called when params change)
  invalidate_xform_cache <- function() rv$xform_cache <- list()

  # Retrieve event indices for a given gate/population combo
  get_gate_idx <- function(gate, popul) {
    if (is.null(gate) || gate == "__root__" || gate == "") {
      return(seq_len(nrow(rv$fcs_expr)))
    }
    rv$fcs_gates_temp[[ gate ]][[ popul ]]$idx
  }

  # Build a data.frame for the scatter/histogram plot
  build_plot_data <- function(gate, popul, marker_x, marker_y) {
    if (is.null(rv$fcs_expr)) return(NULL)
    idx <- get_gate_idx(gate, popul)
    if (is.null(idx)) return(NULL)
    if (length(idx) > input$max_scatter_pts)
      idx <- sample(idx, input$max_scatter_pts)
    data.frame(
      x = get_transformed(marker_x)[idx],
      y = get_transformed(marker_y)[idx]
    )
  }

  # Existing gate names
  defined_gate_names <- function() names(rv$population_gates)

  # Populations within a gate
  populations_in_gate <- function(gname) {
    g <- rv$population_gates[[ gname ]]
    if (is.null(g)) return(character(0))
    names(g$populations)[ names(g$populations) != IGNORE_TOKEN ]
  }

  # --------------------------------------------------------------------------
  # TAB 1: Setup
  # --------------------------------------------------------------------------

  observeEvent(input$btn_load, {
    req(input$fcs_files, input$panel_csv)

    tryCatch({
      # ---- Reset all derived state before (re-)loading ----------------------
      # This ensures a clean slate if the user loads new FCS/compensation files.
      rv$fcs_list         <- list()
      rv$fcs_expr         <- NULL
      rv$panel            <- NULL
      rv$ap               <- NULL
      rv$xform_params     <- NULL
      rv$xform_cache      <- list()
      rv$population_gates <- list()
      rv$activation_gates <- list()
      rv$fcs_gates_temp   <- list()
      rv$current_boundary <- NULL
      rv$calc_gate_type   <- NULL
      rv$click_hint       <- NULL
      rv$gate_number_counter <- 1L
      rv$comp_applied     <- FALSE
      rv$initialized      <- FALSE

      # Load panel
      panel <- read.csv(input$panel_csv$datapath, stringsAsFactors = FALSE)
      stopifnot(all(c("dye", "antigen", "type") %in% colnames(panel)))
      rv$panel <- panel

      # Load FCS
      fcs_list <- setNames(
        lapply(seq_len(nrow(input$fcs_files)), function(i)
          suppressWarnings(read.FCS(input$fcs_files$datapath[i], transformation = FALSE))),
        input$fcs_files$name
      )
      rv$fcs_list <- fcs_list

      # Build combined expression matrix (intersection of channels)
      all_channels <- panel$dye[ panel$dye %in% colnames(exprs(fcs_list[[1]])) ]
      expr_list <- lapply(fcs_list, function(ff) exprs(ff)[, all_channels, drop = FALSE])
      combined  <- do.call(rbind, expr_list)
      colnames(combined) <- panel$antigen[ match(all_channels, panel$dye) ]

      # ---- Apply compensation if a CSV was provided --------------------------
      if (!is.null(input$comp_csv) && nrow(input$comp_csv) > 0) {
        tryCatch({
          fluor_idx      <- which(panel$type %in% c("lineage", "activation"))
          fluor_dyes     <- panel$dye    [ fluor_idx ]
          fluor_antigens <- panel$antigen[ fluor_idx ]

          # Read exactly as read.compensation.csv() does: (this could be changed in both)
          #   - no row.names argument (keep all columns)
          #   - drop the first (label) column with [, -1]
          comp_raw <- read.csv(input$comp_csv$datapath, stringsAsFactors = FALSE)
          comp_mat <- data.matrix(comp_raw[, -1, drop = FALSE])

          # Match panel dye names → CSV column names via make.names()
          dye_names <- make.names(fluor_dyes)
          col_names  <- make.names(colnames(comp_mat))
          reorder    <- match(dye_names, col_names)

          if (any(is.na(reorder))) {
            missing <- fluor_dyes[ is.na(reorder) ]
            showNotification(
              paste0("Compensation: could not match dye(s): ",
                     paste(missing, collapse = ", "), ". NOT applied."),
              type = "warning"
            )
          } else {
            # Reorder rows and columns to match panel fluorescence order,
            # then name with antigen names (matching combined column names)
            comp_sub <- comp_mat[ reorder, reorder ]
            rownames(comp_sub) <- fluor_antigens
            colnames(comp_sub) <- fluor_antigens

            # Apply: compensated = raw %*% solve(comp) ## probably needs to be matched to names, so should use flowCore
            # (equivalent to flowCore::compensate with a named spillover matrix)
            fluor_cols <- fluor_antigens[ fluor_antigens %in% colnames(combined) ]
            if (length(fluor_cols) > 0) {
              comp_inv <- solve(comp_sub[ fluor_cols, fluor_cols ])
              combined[, fluor_cols] <- combined[, fluor_cols] %*% comp_inv
              rv$comp_applied <- TRUE
              showNotification(
                sprintf("Compensation applied (%d channels).", length(fluor_cols)),
                type = "message"
              )
            }
          }
        }, error = function(e) {
          showNotification(
            paste("Compensation error:", conditionMessage(e), "— NOT applied."),
            type = "warning"
          )
          message("Compensation error: ", conditionMessage(e))
        })
      }

      # Store expression matrix (compensated if applicable)
      # xform_cache starts empty so all transforms are computed fresh from this data.
      rv$fcs_expr <- combined

      # Initialise ap
      rv$ap <- build_default_ap(input$cytometer)

      # Initialise transformation params
      fluor_idx <- which(panel$type %in% c("lineage", "activation"))
      defaults  <- rv$ap$transformation.defaults
      rv$xform_params <- data.frame(
        dye       = panel$dye    [ fluor_idx ],
        antigen   = panel$antigen[ fluor_idx ],
        pos       = defaults$pos,
        neg       = defaults$neg,
        width     = defaults$width,
        length    = defaults$length,
        max.range = defaults$max.range,
        stringsAsFactors = FALSE
      )

      # Initialise root gate
      rv$fcs_gates_temp[["all"]] <- list(
        all = list(
          parent.gate  = NULL,
          parent.popul = NULL,
          boundary     = NULL,
          idx          = seq_len(nrow(combined))
        )
      )
      rv$population_gates[["all"]] <- list(
        gate.number   = 0L,
        gate.name     = "all",
        gate.type     = "free",
        gate.marker.x = colnames(combined)[1],
        gate.marker.y = colnames(combined)[2],
        parent.gate   = NULL,
        parent.popul  = NULL,
        populations   = list(all = list(label = "All", label.pos = 0L)),
        algorithm     = "all",
        gate.param    = list(quantile.x = 0, quantile.y = 0),
        stats.parent  = list(gate = "*", popul = "*")
      )
      rv$gate_number_counter <- 1L
      rv$initialized <- TRUE

      # Update UI selectors
      update_all_selectors(session, rv, combined)

    }, error = function(e) {
      showNotification(paste("Load error:", conditionMessage(e)), type = "error")
      message("Load error: ", conditionMessage(e))
    })
  })

  output$setup_status <- renderText({
    if (!rv$initialized) return("Not initialised. Load files and click 'Load & Initialise'.")
    n_files <- length(rv$fcs_list)
    n_chan  <- ncol(rv$fcs_expr)
    n_event <- nrow(rv$fcs_expr)
    comp_str <- if (rv$comp_applied) "Yes" else "No"
    sprintf(
      "Initialised.\n  Files: %d\n  Channels: %d\n  Total events: %d\n  Cytometer: %s\n  Compensation applied: %s",
      n_files, n_chan, n_event, input$cytometer, comp_str
    )
  })

  output$panel_table <- renderDT({
    req(rv$panel)
    datatable(rv$panel, selection = "none",
              options  = list(pageLength = 10, dom = "tp"))
  })

  # --------------------------------------------------------------------------
  # TAB 2: Transformations
  # --------------------------------------------------------------------------

  # Sync slider <-> numeric input (debounced to avoid loop)
  observeEvent(input$xform_width, {
    if (!isTRUE(all.equal(input$xform_width, input$xform_width_num)))
      updateNumericInput(session, "xform_width_num", value = input$xform_width)
  }, ignoreInit = TRUE)
  observeEvent(input$xform_width_num, {
    if (!isTRUE(all.equal(input$xform_width_num, input$xform_width)))
      updateSliderInput(session, "xform_width", value = input$xform_width_num)
  }, ignoreInit = TRUE)

  output$transform_channel_ui <- renderUI({
    req(rv$xform_params)
    selectInput("xform_channel", NULL,
                choices  = rv$xform_params$antigen,
                selected = isolate(input$xform_channel) %||% rv$xform_params$antigen[1])
  })

  # Prev / Next channel buttons
  observeEvent(input$btn_prev_ch, {
    req(rv$xform_params, input$xform_channel)
    chs <- rv$xform_params$antigen
    i   <- match(input$xform_channel, chs)
    if (!is.na(i) && i > 1L)
      updateSelectInput(session, "xform_channel", selected = chs[i - 1L])
  })
  observeEvent(input$btn_next_ch, {
    req(rv$xform_params, input$xform_channel)
    chs <- rv$xform_params$antigen
    i   <- match(input$xform_channel, chs)
    if (!is.na(i) && i < length(chs))
      updateSelectInput(session, "xform_channel", selected = chs[i + 1L])
  })

  # When channel selection changes, populate controls from stored params
  observeEvent(input$xform_channel, {
    req(rv$xform_params)
    row <- rv$xform_params[ rv$xform_params$antigen == input$xform_channel, ]
    if (nrow(row) == 0) return()
    updateSliderInput(session,  "xform_width",     value = row$width[1])
    updateNumericInput(session, "xform_width_num", value = row$width[1])
    updateNumericInput(session, "xform_pos",       value = row$pos[1])
    updateNumericInput(session, "xform_neg",       value = row$neg[1])
  })

  observeEvent(input$btn_apply_xform, {
    req(rv$xform_params, input$xform_channel)
    ch <- input$xform_channel
    rv$xform_params[ rv$xform_params$antigen == ch, "width" ] <- input$xform_width
    rv$xform_params[ rv$xform_params$antigen == ch, "pos"   ] <- input$xform_pos
    rv$xform_params[ rv$xform_params$antigen == ch, "neg"   ] <- input$xform_neg
    rv$xform_cache[[ ch ]] <- NULL   # invalidate cache for this channel only
    showNotification(paste("Transformation updated for", ch), type = "message")
  })

  observeEvent(input$btn_reset_xform, {
    req(rv$ap, rv$xform_params, input$xform_channel)
    ch       <- input$xform_channel
    defaults <- rv$ap$transformation.defaults
    rv$xform_params[ rv$xform_params$antigen == ch, "width" ] <- defaults$width
    rv$xform_params[ rv$xform_params$antigen == ch, "pos"   ] <- defaults$pos
    rv$xform_params[ rv$xform_params$antigen == ch, "neg"   ] <- defaults$neg
    rv$xform_cache[[ ch ]] <- NULL
    updateSliderInput(session,  "xform_width",     value = defaults$width)
    updateNumericInput(session, "xform_width_num", value = defaults$width)
    updateNumericInput(session, "xform_pos",       value = defaults$pos)
    updateNumericInput(session, "xform_neg",       value = defaults$neg)
  })

  # Dynamic plot header
  output$xform_plot_header <- renderUI({
    ch  <- input$xform_channel %||% "—"
    pop <- input$ridge_population %||% "__all__"
    pop_label <- if (pop == "__all__") "all events" else gsub("::", " / ", pop)
    sprintf("Transformation preview  ·  %s  ·  %s  (≤ %s events)",
            ch, pop_label, format(MAX_XFORM_EVENTS, big.mark = ","))
  })

  # Single-channel transformation histogram
  # Axis follows the same convention as build.gate.plot():
  #   - data stored and plotted in transformed space
  #   - breaks are defined in raw space, passed through the transform for
  #     their on-axis positions, labelled with original (raw) numbers
  output$xform_plot <- renderPlot({
    req(rv$fcs_expr, rv$xform_params, rv$ap, input$xform_channel)

    ch  <- input$xform_channel
    row <- rv$xform_params[ rv$xform_params$antigen == ch, ]
    if (nrow(row) == 0) return(NULL)

    # Build a live transform from current slider/field values (not from cache),
    # so the plot updates immediately as the user moves the slider.
    xf <- tryCatch(
      flowjo_biexp(
        channelRange = row$length[1],
        maxValue     = row$max.range[1],
        pos          = input$xform_pos,
        neg          = input$xform_neg,
        widthBasis   = input$xform_width,
        inverse      = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(xf)) {
      return(ggplot() + theme_void() +
               labs(title = "Invalid transformation parameters"))
    }

    # Get events for the selected population
    pop_sel <- input$ridge_population %||% "__all__"
    if (pop_sel == "__all__") {
      idx <- seq_len(nrow(rv$fcs_expr))
    } else {
      parts <- strsplit(pop_sel, "::")[[1]]
      idx   <- get_gate_idx(parts[1], parts[2])
    }
    if (is.null(idx) || length(idx) == 0) return(NULL)

    # Downsample
    if (length(idx) > MAX_XFORM_EVENTS)
      idx <- sample(idx, MAX_XFORM_EVENTS)

    # Raw values → transformed values for plotting
    raw_vals  <- rv$fcs_expr[ idx, ch ]
    tran_vals <- xf(raw_vals)

    # ---- Axis breaks / labels / limits (matching build.gate.plot approach) ----
    # Breaks defined in raw space using the cytometer's break.fluor defaults.
    fl  <- rv$ap$figure.limits
    ap  <- rv$ap

    # Raw break positions (same as ap$gate.defaults$break.fluor if present,
    # otherwise fall back to break.x * break.x.factor)
    raw_breaks <- if (!is.null(ap$gate.defaults$break.fluor)) {
      ap$gate.defaults$break.fluor
    } else {
      ap$gate.defaults$break.x * ap$gate.defaults$break.x.factor
    }

    # Transform the break positions and the axis limits
    tran_breaks <- xf(raw_breaks)
    tran_lim    <- xf(c(fl$expr.min, fl$expr.max))

    # Human-readable labels: use K / M suffix like build.gate.plot's suffix()
    fmt_break <- function(v) {
      if      (abs(v) >= 1e6)  paste0(v / 1e6, "M")
      else if (abs(v) >= 1e3)  paste0(v / 1e3, "K")
      else                     as.character(v)
    }
    axis_labels <- sapply(raw_breaks, fmt_break)

    # ---- Plot -----------------------------------------------------------------
    df <- data.frame(x = tran_vals)

    ggplot(df, aes(x)) +
      geom_histogram(bins = 256, fill = "#2C6FAC", colour = NA, alpha = 0.85,
                     na.rm = TRUE) +
      scale_x_continuous(
        name   = ch,
        breaks = tran_breaks,
        labels = axis_labels,
        limits = tran_lim
      ) +
      scale_y_continuous(name = "Count", expand = expansion(mult = c(0, 0.05))) +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text          = element_text(size = 12),
        axis.title         = element_text(size = 13)
      )
  })

  output$btn_save_xform <- downloadHandler(
    filename = function() "transformation_params.csv",
    content  = function(f) {
      req(rv$xform_params)
      write.csv(rv$xform_params[, c("dye","pos","neg","width","length","max.range")],
                f, row.names = FALSE)
    }
  )

  # --------------------------------------------------------------------------
  # TAB 3: Gate Editor
  # --------------------------------------------------------------------------

  # Gate type → algorithm UI
  output$algorithm_ui <- renderUI({
    gt <- req(input$gate_type)
    if (gt == "1dsep") {
      tagList(
        selectInput("algorithm_1d", "Algorithm", choices = ALGORITHMS_1D),
        selectInput("sep_axis", "Separation axis",
                    choices = c("X-axis" = "1", "Y-axis" = "2"))
      )
    } else if (gt == "2dsep") {
      selectInput("algorithm_2d", "Algorithm (x:y)", choices = ALGORITHMS_2D)
    } else if (gt == "replicate") {
      selectInput("origin_gate", "Replicate from gate",
                  choices = defined_gate_names())
    } else if (gt == "singlets") {
      p(class = "text-muted small",
        "Singlet gate: fits a robust linear model to FSC-A vs FSC-H (or the ",
        "selected X vs Y markers) and returns a polygon boundary. ",
        "X-axis marker is used as 'area', Y-axis as 'height'.")
    } else {
      p(class = "text-muted small", "No algorithm selection for free/all-events gate.")
    }
  })

  # Gate type → population name fields
  output$popul_fields_ui <- renderUI({
    gt <- req(input$gate_type)
    if (gt == "1dsep") {
      tagList(
        textInput("popul_neg", "Negative population name", value = "neg"),
        textInput("popul_pos", "Positive population name", value = "pos")
      )
    } else if (gt == "2dsep") {
      tagList(
        textInput("popul_dn",   "Double-negative (DN)",  value = "dn"),
        textInput("popul_xpos", "X-positive (Xpos)",     value = "xpos"),
        textInput("popul_ypos", "Y-positive (Ypos)",     value = "ypos"),
        textInput("popul_dp",   "Double-positive (DP)",  value = "dp")
      )
    } else if (gt == "singlets") {
      textInput("popul_single", "Population name", value = "singlets")
    } else {
      textInput("popul_single", "Population name", value = "all")
    }
  })

  # Update parent gate → parent population cascade
  observeEvent(defined_gate_names(), {
    gnames <- defined_gate_names()
    updateSelectInput(session, "parent_gate",
                      choices  = c("(none — root)" = "__root__", gnames),
                      selected = if (length(gnames) > 0) tail(gnames, 1) else "__root__")
    updateSelectInput(session, "ref_gate",
                      choices  = c("(none)" = "", gnames))
    updateSelectInput(session, "stats_parent_gate", choices = gnames)
  }, ignoreNULL = FALSE)

  observeEvent(input$parent_gate, {
    gn <- input$parent_gate
    if (gn == "__root__" || is.null(gn)) {
      pops <- "(root)"
      updateSelectInput(session, "parent_popul", choices = pops)
    } else {
      pops <- populations_in_gate(gn)
      updateSelectInput(session, "parent_popul",
                        choices  = pops,
                        selected = if (length(pops)) pops[1] else NULL)
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$ref_gate, {
    gn <- input$ref_gate
    if (is.null(gn) || gn == "") {
      updateSelectInput(session, "ref_popul", choices = c())
    } else {
      pops <- populations_in_gate(gn)
      updateSelectInput(session, "ref_popul", choices = pops)
    }
  })

  # Update marker selectors when panel loads — default to FSC-A / FSC-H
  observeEvent(rv$panel, {
    req(rv$panel)
    all_m   <- get_all_params()
    scatter <- get_scatter_params()
    default_x <- {
      fsc_a <- scatter[ grepl("FSC-A", scatter, ignore.case = TRUE) ]
      if (length(fsc_a)) fsc_a[1] else if (length(scatter)) scatter[1] else all_m[1]
    }
    default_y <- {
      fsc_h <- scatter[ grepl("FSC-H", scatter, ignore.case = TRUE) ]
      if (length(fsc_h)) fsc_h[1] else if (length(scatter) >= 2) scatter[2] else all_m[2]
    }
    updateSelectInput(session, "gate_marker_x", choices = all_m, selected = default_x)
    updateSelectInput(session, "gate_marker_y", choices = all_m, selected = default_y)
  })

  # ---- Gate plot: only recalculates on explicit button press -----------------
  # Driven by a counter incremented by "Calculate gate" and "Refresh plot".
  # Selector changes (marker, parent, gate type) do NOT re-render the plot.

  gate_plot_trigger <- reactiveVal(0L)
  observeEvent(input$btn_calculate_gate, { gate_plot_trigger(gate_plot_trigger() + 1L) })
  observeEvent(input$btn_refresh_plot,   { gate_plot_trigger(gate_plot_trigger() + 1L) })

  # Snapshot inputs at the moment of the trigger
  gate_plot_snap <- eventReactive(gate_plot_trigger(), {
    list(
      marker_x     = input$gate_marker_x,
      marker_y     = input$gate_marker_y,
      parent_gate  = input$parent_gate,
      parent_popul = input$parent_popul,
      gate_type    = input$gate_type,
      max_pts      = input$max_scatter_pts
    )
  }, ignoreInit = TRUE)

  output$gate_plot <- renderPlot({
    snap <- req(gate_plot_snap())
    req(rv$initialized, rv$ap)

    mx <- snap$marker_x
    my <- snap$marker_y
    pg <- snap$parent_gate %||% "__root__"
    pp <- snap$parent_popul
    if (pg == "__root__") { pg <- "all"; pp <- "all" }

    # ---- events --------------------------------------------------------------
    idx <- get_gate_idx(pg, pp)
    if (is.null(idx) || length(idx) == 0) {
      return(ggplot() + theme_void() +
               labs(title = "No events in selected parent population"))
    }
    if (length(idx) > snap$max_pts)
      idx <- sample(idx, snap$max_pts)

    df <- data.frame(
      x = get_transformed(mx)[idx],
      y = get_transformed(my)[idx]
    )

    # ---- axis scaling: same convention as build.gate.plot() ------------------
    # Breaks are defined in raw (pre-transform) space, then passed through the
    # transform to get their on-axis positions, but labelled with the raw values.

    fmt_break <- function(v) {
      if      (abs(v) >= 1e6)  paste0(v / 1e6, "M")
      else if (abs(v) >= 1e3)  paste0(v / 1e3, "K")
      else                     as.character(v)
    }

    axis_scale <- function(marker) {
      ap  <- rv$ap
      fl  <- ap$figure.limits
      gd  <- ap$gate.defaults
      is_sc <- marker %in% get_scatter_params()

      if (is_sc) {
        # Scatter channels are linear — no transform applied
        brk_raw <- gd$break.scatter * gd$break.scatter.factor
        lim_raw <- c(fl$scatter.min, fl$scatter.max)
        list(breaks = brk_raw,
             labels = sapply(brk_raw, fmt_break),
             limits = lim_raw)
      } else {
        # Fluorescence: data is already in transformed space via get_transformed()
        # so we must transform the break positions and limits too.
        brk_raw <- if (!is.null(gd$break.fluor)) gd$break.fluor
                   else gd$break.x * gd$break.x.factor
        lim_raw <- c(fl$expr.min, fl$expr.max)
        row <- rv$xform_params[ rv$xform_params$antigen == marker, ]
        if (nrow(row) == 0) {
          # No transform info: fall back to untransformed linear labels
          return(list(breaks = brk_raw,
                      labels = sapply(brk_raw, fmt_break),
                      limits = lim_raw))
        }
        xf <- tryCatch(
          flowjo_biexp(
            channelRange = row$length[1], maxValue = row$max.range[1],
            pos = row$pos[1], neg = row$neg[1], widthBasis = row$width[1],
            inverse = FALSE
          ),
          error = function(e) NULL
        )
        if (is.null(xf)) {
          return(list(breaks = brk_raw,
                      labels = sapply(brk_raw, fmt_break),
                      limits = lim_raw))
        }
        list(breaks = xf(brk_raw),
             labels = sapply(brk_raw, fmt_break),
             limits = xf(lim_raw))
      }
    }

    sc_x <- axis_scale(mx)
    sc_y <- axis_scale(my)

    # ---- plot ----------------------------------------------------------------
    # use_histogram only for genuine 1D cases: same marker on both axes,
    # or a true 1dsep gate type. free/singlets/2dsep always use hex-bin.
    use_histogram <- (mx == my) || snap$gate_type == "1dsep"

    if (use_histogram) {
      p <- ggplot(df, aes(x)) +
        geom_histogram(bins = 256, fill = "#2C6FAC", colour = NA, alpha = 0.85,
                       na.rm = TRUE) +
        scale_x_continuous(name   = mx,
                           breaks = sc_x$breaks,
                           labels = sc_x$labels,
                           limits = sc_x$limits) +
        scale_y_continuous(name   = "Count",
                           expand = expansion(mult = c(0, 0.05)))
    } else {
      p <- ggplot(df, aes(x, y)) +
        geom_hex(bins = 120, na.rm = TRUE) +
        scale_fill_viridis_c(option = "magma", trans = "log1p") +
        scale_x_continuous(name   = mx,
                           breaks = sc_x$breaks,
                           labels = sc_x$labels,
                           limits = sc_x$limits) +
        scale_y_continuous(name   = my,
                           breaks = sc_y$breaks,
                           labels = sc_y$labels,
                           limits = sc_y$limits)
    }

    p <- p +
      labs(title = sprintf("Parent: %s / %s  (%s events)",
                           pg, pp, format(length(idx), big.mark = ","))) +
      theme_bw(base_size = 13) +
      theme(legend.position  = "none",
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())

    # Overlay click hint (coordinates in transformed/display space)
    if (!is.null(rv$click_hint)) {
      p <- p + geom_vline(xintercept = rv$click_hint$x,
                          colour = "#e07b39", linewidth = 1, linetype = "dashed")
      if (snap$gate_type == "2dsep")
        p <- p + geom_hline(yintercept = rv$click_hint$y,
                             colour = "#e07b39", linewidth = 1, linetype = "dashed")
    }

    # Overlay calculated boundary — isolate() so boundary changes don't
    # trigger an extra render outside of explicit button presses.
    # Only overlay if the stored calc_gate_type matches the current plot type.
    bd        <- isolate(rv$current_boundary)
    bd_gt     <- isolate(rv$calc_gate_type)   # "free", "1dsep", or "2dsep"
    show_bd   <- !is.null(bd) && !is.null(bd_gt) && bd_gt != "free"

    if (show_bd && !use_histogram) {
      # 2D boundary: polygon/path drawn on hex-bin scatter
      bd_df <- as.data.frame(bd)
      colnames(bd_df) <- c("bx", "by")
      p <- p + geom_path(data = bd_df, aes(x = bx, y = by),
                         colour = "#27ae60", linewidth = 1.2,
                         inherit.aes = FALSE)
    } else if (show_bd && use_histogram && bd_gt == "1dsep") {
      # 1D boundary: draw the threshold as vertical lines on histogram
      # Boundary is a 4×2 rectangle; x-threshold is at bd[1,1] and bd[3,1]
      x_thr <- sort(unique(round(bd[, 1], 6)))
      if (length(x_thr) >= 2) {
        p <- p + geom_vline(xintercept = x_thr[2],   # the gate threshold
                            colour = "#27ae60", linewidth = 1.2)
      }
    }

    p
  })

  # Capture click → store hint (in display/transformed space — same as plot coords)
  observeEvent(input$gate_plot_click, {
    rv$click_hint <- data.frame(x = input$gate_plot_click$x,
                                y = input$gate_plot_click$y)
  })
  observeEvent(input$gate_plot_dblclick, {
    rv$click_hint <- NULL
  })

  output$click_info <- renderText({
    if (is.null(rv$click_hint)) return("No gate hint set.")
    sprintf("Gate hint: x = %.3f, y = %.3f", rv$click_hint$x, rv$click_hint$y)
  })

  # Build gate param list from UI
  build_gate_param <- function() {
    p <- list(
      region.min = input$region_min,
      region.max = input$region_max
    )
    if (input$gate_type == "1dsep") {
      p$sep.axis <- as.integer(input$sep_axis %||% "1")
    }
    if (!is.null(rv$click_hint)) {
      p$user.hint.x <- rv$click_hint$x
      p$user.hint.y <- rv$click_hint$y
    }
    p
  }

  build_populations_list <- function() {
    gt <- input$gate_type
    if (gt == "1dsep") {
      neg_name <- input$popul_neg %||% "neg"
      pos_name <- input$popul_pos %||% "pos"
      populs <- list()
      populs[[ neg_name ]] <- list(label = neg_name, label.pos = 2L)
      populs[[ pos_name ]] <- list(label = pos_name, label.pos = 4L)
      populs
    } else if (gt == "2dsep") {
      nms <- c(input$popul_dn   %||% "dn",
               input$popul_xpos %||% "xpos",
               input$popul_ypos %||% "ypos",
               input$popul_dp   %||% "dp")
      pos  <- c(1L, 3L, 7L, 9L)
      setNames(
        lapply(seq_along(nms), function(i)
          list(label = nms[i], label.pos = pos[i])),
        nms
      )
    } else if (gt == "singlets") {
      nm <- input$popul_single %||% "singlets"
      setNames(list(list(label = nm, label.pos = 4L)), nm)
    } else {
      nm <- input$popul_single %||% "all"
      setNames(list(list(label = nm, label.pos = 0L)), nm)
    }
  }

  build_algorithm <- function() {
    gt <- input$gate_type
    if      (gt == "1dsep")    input$algorithm_1d   %||% "tail"
    else if (gt == "2dsep") {
      parts <- strsplit(input$algorithm_2d %||% "bimodal:bimodal", ":")[[1]]
      parts
    }
    else if (gt == "replicate") "replicate"
    else if (gt == "singlets")  "singlets"
    else                        "all"
  }

  # Calculate gate button
  observeEvent(input$btn_calculate_gate, {
    if (!isTRUE(rv$initialized)) {
      output$gate_calc_status <- renderText("Not initialised — load data first.")
      return()
    }
    output$gate_calc_status <- renderText("Calculating…")

    tryCatch({
      # ---- resolve parent population indices ---------------------------------
      pg <- input$parent_gate  %||% "__root__"
      pp <- input$parent_popul %||% "all"
      if (pg == "__root__") { pg <- "all"; pp <- "all" }

      idx <- rv$fcs_gates_temp[[ pg ]][[ pp ]]$idx
      if (is.null(idx)) idx <- seq_len(nrow(rv$fcs_expr))
      if (length(idx) == 0L)
        stop(sprintf("Parent '%s / %s' contains no events.", pg, pp))

      # ---- markers -----------------------------------------------------------
      mx <- input$gate_marker_x
      my <- input$gate_marker_y
      if (is.null(mx) || !nzchar(mx)) stop("X marker not selected.")
      if (is.null(my) || !nzchar(my)) stop("Y marker not selected.")
      if (!mx %in% colnames(rv$fcs_expr))
        stop(sprintf("Marker '%s' not found in expression matrix.", mx))
      if (!my %in% colnames(rv$fcs_expr))
        stop(sprintf("Marker '%s' not found in expression matrix.", my))

      # ---- build 2-column gate data matrix -----------------------------------
      gate_data <- matrix(
        c(get_transformed(mx)[idx], get_transformed(my)[idx]),
        ncol = 2L,
        dimnames = list(NULL, c(mx, my))
      )

      # ---- gate definition inputs --------------------------------------------
      gt          <- input$gate_type %||% "free"
      popul_names <- names(build_populations_list())
      algo        <- build_algorithm()
      gp          <- build_gate_param()

      if (length(popul_names) == 0L)
        stop("No population names defined.")

      # "singlets" is a GUI alias — pipeline uses gate.type="free" + algo="singlets"
      calc_gate_type <- if (gt == "singlets") "free" else gt

      # ---- build ap ----------------------------------------------------------
      ap_calc <- rv$ap
      if (is.null(ap_calc))
        stop("ap not initialised — load data first.")
      for (k in names(gp)) ap_calc$gate.defaults[[k]] <- gp[[k]]
      if (gt == "singlets")
        ap_calc$singlet.scatter.parameter <- c(mx, my)

      p_resolved <- modifyList(ap_calc$gate.defaults, gp)

      # ---- call calculate.gate.* ---------------------------------------------
      calc_fn_name <- paste0("calculate.gate.", calc_gate_type)
      if (!exists(calc_fn_name, mode = "function"))
        stop(sprintf(
          "%s() not found — source AutoGating R files before running the app.",
          calc_fn_name))

      boundaries <- do.call(calc_fn_name, list(
        gate.data   = gate_data,
        popul.names = popul_names,
        algorithm   = algo,
        p           = p_resolved,
        gate.number = rv$gate_number_counter,
        gate.name   = input$gate_name %||% "unnamed",
        fcs.gates   = rv$fcs_gates_temp,
        figure.dir  = tempdir(),
        ap          = ap_calc
      ))

      rv$current_boundary <- boundaries[[ popul_names[1] ]]
      rv$calc_gate_type   <- calc_gate_type   # "free", "1dsep", or "2dsep"
      output$gate_calc_status <- renderText(
        sprintf("OK — %d population(s), %s events in parent.",
                length(boundaries), format(length(idx), big.mark = ","))
      )

    }, error = function(e) {
      msg <- conditionMessage(e)
      output$gate_calc_status <- renderText(paste0("Error: ", msg))
      message("Gate calculation error: ", msg)   # also print to R console
    })
  })

  # Gate param preview
  output$gate_param_preview <- renderText({
    gt <- input$gate_type %||% "?"
    gd_preview <- list(
      gate.name    = input$gate_name   %||% "(unnamed)",
      gate.type    = if (gt == "singlets") "free" else gt,
      gate.marker.x = input$gate_marker_x %||% "?",
      gate.marker.y = input$gate_marker_y %||% "?",
      parent.gate  = input$parent_gate  %||% NULL,
      parent.popul = input$parent_popul %||% NULL,
      algorithm    = build_algorithm(),
      gate.param   = build_gate_param(),
      populations  = build_populations_list()
    )
    deparse_compact(gd_preview)
  })

  # Add gate to population_gates list
  observeEvent(input$btn_add_gate, {
    req(rv$initialized, input$gate_name != "")

    gname <- trimws(input$gate_name)
    if (gname %in% names(rv$population_gates) && gname != "all") {
      showNotification(
        paste("Gate", gname, "already exists. Delete first or rename."),
        type = "warning"
      )
      return()
    }

    pg <- input$parent_gate
    pp <- input$parent_popul
    if (is.null(pg) || pg == "__root__") { pg <- NULL; pp <- NULL }

    stats_parent <- if (input$stats_parent_immediate) {
      list(gate = "*", popul = "*")
    } else {
      list(gate  = input$stats_parent_gate  %||% "*",
           popul = input$stats_parent_popul %||% "*")
    }

    gt <- input$gate_type

    new_gate <- list(
      gate.number   = rv$gate_number_counter,
      gate.name     = gname,
      # "singlets" is a GUI convenience — the pipeline uses gate.type = "free"
      gate.type     = if (gt == "singlets") "free" else gt,
      gate.marker.x = input$gate_marker_x,
      gate.marker.y = input$gate_marker_y,
      parent.gate   = pg,
      parent.popul  = pp,
      populations   = build_populations_list(),
      algorithm     = build_algorithm(),
      gate.param    = build_gate_param(),
      stats.parent  = stats_parent
    )

    rv$population_gates[[ gname ]] <- new_gate
    rv$gate_number_counter         <- rv$gate_number_counter + 1L
    rv$current_boundary            <- NULL
    rv$calc_gate_type              <- NULL
    rv$click_hint                  <- NULL

    # Propagate idx into fcs_gates_temp if boundary was calculated
    # (In full implementation, call apply.gate() here)

    # Update the ridge-plot population selector
    gate_popul_choices <- lapply(names(rv$population_gates), function(gn) {
      pops <- populations_in_gate(gn)
      setNames(paste0(gn, "::", pops), paste0(gn, " / ", pops))
    })
    all_choices <- c("All events (ungated)" = "__all__",
                     unlist(gate_popul_choices))
    updateSelectInput(session, "ridge_population",
                      choices = all_choices, selected = "__all__")

    showNotification(paste("Gate", gname, "added."), type = "message")
  })

  # --------------------------------------------------------------------------
  # TAB 4: Gate Tree
  # --------------------------------------------------------------------------

  output$gate_tree_table <- renderDT({
    req(rv$initialized)
    gates <- rv$population_gates
    if (length(gates) == 0)
      return(datatable(data.frame(message = "No gates defined.")))

    df <- do.call(rbind, lapply(gates, function(g) {
      data.frame(
        number   = g$gate.number,
        name     = g$gate.name,
        type     = g$gate.type,
        parent   = paste(g$parent.gate   %||% "—",
                         g$parent.popul  %||% "", sep = " / "),
        x        = g$gate.marker.x,
        y        = g$gate.marker.y,
        pops     = paste(names(g$populations), collapse = ", "),
        algo     = paste(g$algorithm, collapse = ":"),
        stringsAsFactors = FALSE
      )
    }))
    df <- df[ order(df$number), ]
    datatable(df, selection = "single",
              options  = list(pageLength = 20, dom = "tp"))
  }, server = FALSE)

  output$gate_detail_preview <- renderText({
    sel <- input$gate_tree_table_rows_selected
    if (is.null(sel) || length(rv$population_gates) == 0) return("Select a row.")
    gnames <- names(rv$population_gates)
    gnames_ordered <- gnames[ order(sapply(rv$population_gates, `[[`, "gate.number")) ]
    if (sel > length(gnames_ordered)) return("")
    g <- rv$population_gates[[ gnames_ordered[ sel ] ]]
    deparse_compact(g)
  })

  observeEvent(input$btn_delete_gate, {
    sel <- input$gate_tree_table_rows_selected
    req(sel)
    gnames <- names(rv$population_gates)
    gnames_ordered <- gnames[ order(sapply(rv$population_gates, `[[`, "gate.number")) ]
    if (sel > length(gnames_ordered)) return()
    del_name <- gnames_ordered[ sel ]
    if (del_name == "all") {
      showNotification("Cannot delete root gate.", type = "warning")
      return()
    }
    rv$population_gates[[ del_name ]] <- NULL
    rv$fcs_gates_temp[[ del_name ]]   <- NULL
    showNotification(paste("Deleted gate:", del_name), type = "message")
  })

  # --------------------------------------------------------------------------
  # TAB 5: Export
  # --------------------------------------------------------------------------

  output$export_pop_preview <- renderText({
    if (length(rv$population_gates) == 0) return("# No gates defined.")
    varname <- input$export_pop_varname %||% "fcs.population.gates"
    gates_list <- unname(rv$population_gates)
    gates_ordered <- gates_list[ order(sapply(gates_list, `[[`, "gate.number")) ]
    code <- sprintf(
      "%s <- %s",
      varname,
      deparse_list(gates_ordered)
    )
    code
  })

  output$btn_export_pop <- downloadHandler(
    filename = function() {
      ext <- if (input$export_pop_format == "rds") "rds" else "r"
      paste0("population_gates.", ext)
    },
    content  = function(f) {
      gates_list    <- unname(rv$population_gates)
      gates_ordered <- gates_list[ order(sapply(gates_list, `[[`, "gate.number")) ]
      if (input$export_pop_format == "rds") {
        saveRDS(gates_ordered, f)
      } else {
        varname <- input$export_pop_varname %||% "fcs.population.gates"
        code    <- sprintf("%s <- %s\n", varname, deparse_list(gates_ordered))
        writeLines(code, f)
      }
    }
  )

  output$btn_export_act <- downloadHandler(
    filename = function() {
      ext <- if (input$export_act_format == "rds") "rds" else "r"
      paste0("activation_gates.", ext)
    },
    content  = function(f) {
      if (input$export_act_format == "rds") {
        saveRDS(rv$activation_gates, f)
      } else {
        varname <- input$export_act_varname %||% "fcs.activation.gates"
        code    <- sprintf("%s <- %s\n", varname, deparse_list(rv$activation_gates))
        writeLines(code, f)
      }
    }
  )

  output$export_xform_table <- renderDT({
    req(rv$xform_params)
    datatable(rv$xform_params, selection = "none",
              options  = list(pageLength = 15, dom = "tp"))
  })

  output$btn_export_xform_csv <- downloadHandler(
    filename = function() "transformation_params.csv",
    content  = function(f) {
      req(rv$xform_params)
      write.csv(
        rv$xform_params[, c("dye","pos","neg","width","length","max.range")],
        f, row.names = FALSE
      )
    }
  )
}

# ==============================================================================
# Helper utilities
# ==============================================================================

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0 && a != "") a else b

# Build a minimal ap structure for a given cytometer string
build_default_ap <- function(cytometer) {
  defaults_by_cytometer <- list(
    bd_fortessa = list(
      transformation.defaults = list(length=256, max.range=262144,
                                     pos=4.418539922, neg=0, width=-100),
      figure.limits = list(scatter.min=0, scatter.max=265000,
                           expr.min=-2000, expr.max=265000)
    ),
    aurora = list(
      transformation.defaults = list(length=256, max.range=4194304,
                                     pos=5.62, neg=0, width=-1000),
      figure.limits = list(scatter.min=0, scatter.max=4200000,
                           expr.min=-5000, expr.max=4200000)
    ),
    id7000 = list(
      transformation.defaults = list(length=256, max.range=1000000,
                                     pos=5, neg=0, width=-250),
      figure.limits = list(scatter.min=0, scatter.max=1050000,
                           expr.min=-2000, expr.max=1050000)
    ),
    symphony = list(
      transformation.defaults = list(length=256, max.range=262144,
                                     pos=4.418539922, neg=0, width=-100),
      figure.limits = list(scatter.min=0, scatter.max=265000,
                           expr.min=-2000, expr.max=265000)
    )
  )
  base <- list(
    cytometer    = cytometer,
    ignore       = "*",
    parent.string = "parent",
    stats.separator = " :in: ",
    gate.number.width = 3L,
    gate.defaults = list(
      quantile.x = 0, quantile.y = 0.0001,
      limit.xmin = -0.01, limit.xmax = 1.01,
      limit.ymin = -0.01, limit.ymax = 1.01,
      region.xmin = 0.01, region.xmax = 0.99,
      region.ymin = 0.01, region.ymax = 0.99,
      region.min  = 0.01, region.max  = 0.99,
      density.quantile.x = 0, density.quantile.y = 0.0001,
      density.trimneg.x = TRUE, density.trimneg.y = FALSE,
      break.x = c(-1,0,1,10,100), break.x.factor = 1000,
      break.y = c(-1,0,1,10,100), break.y.factor = 1000,
      label.size = 4,
      algorithm.1dsep = "tail",
      algorithm.2dsep = c("bimodal","bimodal"),
      default.tail.probability = 0.99,
      sep.axis = 1L
    )
  )
  overrides <- defaults_by_cytometer[[ cytometer ]]
  if (!is.null(overrides)) {
    base$transformation.defaults <- overrides$transformation.defaults
    base$figure.limits           <- overrides$figure.limits
  }
  base
}

# Update all selector inputs after data is loaded
update_all_selectors <- function(session, rv, combined) {
  all_m   <- colnames(combined)
  scatter <- rv$panel$antigen[ rv$panel$type == "scatter" ]

  # Default x = FSC-A (or first scatter), y = FSC-H (or second scatter)
  # This is the standard starting point for singlet gating.
  default_x <- {
    fsc_a <- scatter[ grepl("FSC-A", scatter, ignore.case = TRUE) ]
    if (length(fsc_a)) fsc_a[1] else if (length(scatter)) scatter[1] else all_m[1]
  }
  default_y <- {
    fsc_h <- scatter[ grepl("FSC-H", scatter, ignore.case = TRUE) ]
    if (length(fsc_h)) fsc_h[1] else if (length(scatter) >= 2) scatter[2] else all_m[2]
  }

  updateSelectInput(session, "gate_marker_x", choices = all_m, selected = default_x)
  updateSelectInput(session, "gate_marker_y", choices = all_m, selected = default_y)
  updateSelectInput(session, "parent_gate",
                    choices  = c("(none — root)" = "__root__", "all"),
                    selected = "__root__")
  updateSelectInput(session, "ridge_population",
                    choices  = c("All events (ungated)" = "__all__"),
                    selected = "__all__")
}

# Compact deparse for display
deparse_compact <- function(x) {
  paste(deparse(x, control = c("keepNA","keepInteger")), collapse = "\n")
}

# Deparse a list as a formatted R list literal
deparse_list <- function(x) {
  paste(deparse(x, control = c("keepNA","keepInteger"), width.cutoff = 72),
        collapse = "\n")
}

# ==============================================================================
# Run
# ==============================================================================

shinyApp(ui = ui, server = server)
