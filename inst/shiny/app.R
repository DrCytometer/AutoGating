# app.R
#
# The AutoGating app: an interactive editor for gating models (see
# ?autogating.app).  autogating.app() evaluates this file in an environment
# whose parent is the package namespace, with app.start holding the
# experiment, parameters and model passed in, and runs the shinyApp it
# returns.
#
# State: the experiment (fcs.setup), the parameters (ap), the training set
# (pooled training events), the model and the gating result of the model on
# the training events.  Every change to the model (a gate edited, added,
# deleted or recalculated, a transform changed) recalculates the gates
# without a stored boundary on the training events; the previous models are
# kept for Undo.

library( shiny )
library( bslib )
library( DT )
library( ggplot2 )
library( shinyFiles )

start <- if ( exists( "app.start", envir = environment(), inherits = FALSE ) ) app.start else list()
start.ap <- start$ap %||% start$fcs.setup$ap
start.value <- function( key, default = "" ) {
    v <- start.ap[[ key ]]
    if ( is.null( v ) ) default else v
}

region.id <- function( region ) gsub( "+", "p", gsub( "-", "m", region, fixed = TRUE ), fixed = TRUE )
region.labels <- c( below = "Below the threshold", above = "Above the threshold",
                    `x-y-` = "X- Y-", `x+y-` = "X+ Y-", `x-y+` = "X- Y+", `x+y+` = "X+ Y+",
                    inside = "Inside" )
path.types <- list(
    panel_file = "csv", experiment_file = "csv", compensation_file = "csv", transformation_file = "csv",
    agmodel_file = "agmodel", flowjo_file = c( "wsp", "wspt" ),
    population_gates_file = c( "r", "R", "rds" ), activation_gates_file = c( "r", "R", "rds" ),
    tr_file = "csv", model_file = "agmodel" )

path.input <- function( kind, id, label, value = "" ) {
    button <- paste0( id, "_browse" )
    title <- paste( "Choose:", label )
    div( class = "d-flex gap-1 align-items-end",
         div( style = "flex: 1", textInput( id, label, value, width = "100%" ) ),
         div( class = "mb-3",
              switch( kind,
                      dir  = shinyDirButton( button, "Browse", title ),
                      file = shinyFilesButton( button, "Browse", title, multiple = FALSE ),
                      save = shinySaveButton( button, "Browse", title,
                                              filetype = stats::setNames( list( path.types[[ id ]] ), path.types[[ id ]][ 1L ] ) ) ) ) )
}

default.algorithm <- c( `1dsep` = "mixture", `2dsep` = "mixture", singlets = "singlets",
                        free = "cells", replicate = "replicate" )


# ==============================================================================
# User interface
# ==============================================================================

ui <- page_navbar(
    title = "AutoGating",
    theme = bs_theme( version = 5 ),
    fillable = FALSE,

    nav_panel(
        "Data",
        layout_columns(
            col_widths = c( 6, 6 ),
            card(
                card_header( "Experiment" ),
                selectInput( "cytometer", "Cytometer",
                             c( Aurora = "aurora", ID7000 = "id7000", Symphony = "symphony",
                                Fortessa = "fortessa" ),
                             selected = tolower( start.value( "cytometer", "aurora" ) ) ),
                path.input( "dir", "data_dir", "Folder of FCS files", start.value( "data.dir" ) ),
                path.input( "file", "panel_file", "Panel CSV (dye, antigen, type)", start.value( "panel.filename" ) ),
                path.input( "file", "experiment_file", "Experiment CSV (filename, sample, gate.calculation)",
                           start.value( "experiment.filename" ) ),
                path.input( "file", "compensation_file", "Compensation CSV (optional)",
                           start.value( "compensation.filename" ) ),
                path.input( "file", "transformation_file", "Transformation CSV (optional)",
                           start.value( "transformation.filename" ) ),
                numericInput( "events_per_sample", "Training events per sample",
                              start.value( "events.per.sample", 50000L ), min = 1000, step = 1000 ),
                div( class = "d-flex flex-wrap gap-1",
                     actionButton( "load_data", "Load experiment", class = "btn-primary" ),
                     actionButton( "resample", "Resample training events" ) )
            ),
            card(
                card_header( "Gating model" ),
                radioButtons( "model_source", NULL,
                              c( "Empty" = "empty", ".agmodel file" = "agmodel",
                                 "FlowJo workspace" = "flowjo", "R gate lists" = "rlist" ) ),
                conditionalPanel(
                    "input.model_source == 'agmodel'",
                    path.input( "file", "agmodel_file", ".agmodel file", start.value( "model.filename" ) ) ),
                conditionalPanel(
                    "input.model_source == 'flowjo'",
                    path.input( "file", "flowjo_file", "FlowJo workspace (.wsp or .wspt)" ),
                    checkboxInput( "flowjo_convert", "Convert imported gates to calculated gates", TRUE ) ),
                conditionalPanel(
                    "input.model_source == 'rlist'",
                    path.input( "file", "population_gates_file", "Population gates (.r or .rds)",
                               start.value( "population.gates.filename" ) ),
                    path.input( "file", "activation_gates_file", "Activation gates (optional)",
                               start.value( "activation.gates.filename" ) ) ),
                actionButton( "load_model", "Load model", class = "btn-primary" ),
                hr(),
                verbatimTextOutput( "data_summary" )
            )
        ),
        card( card_header( "Messages" ), verbatimTextOutput( "log" ) )
    ),

    nav_panel(
        "Transforms",
        layout_sidebar(
            sidebar = sidebar(
                width = 320,
                selectInput( "tr_channel", "Channel", NULL ),
                selectInput( "tr_population", "Events shown", NULL ),
                radioButtons( "tr_id", "Transform", c( Logicle = "1", Log = "2", Linear = "0" ), inline = TRUE ),
                numericInput( "tr_t", "Top of the scale, T (raw units)", NA ),
                conditionalPanel(
                    "input.tr_id == '1'",
                    numericInput( "tr_w", "Linear width, W (decades)", NA, step = 0.05 ),
                    numericInput( "tr_m", "Width of the display, M (decades)", NA, step = 0.1 ) ),
                conditionalPanel(
                    "input.tr_id == '2'",
                    numericInput( "tr_logm", "Width of the display, M (decades)", NA, step = 0.1 ) ),
                conditionalPanel(
                    "input.tr_id == '0'",
                    numericInput( "tr_a", "Offset, A (raw units)", 0 ) ),
                actionButton( "tr_apply", "Apply to channel", class = "btn-primary" ),
                actionButton( "tr_reset", "Reset" ),
                hr(),
                path.input( "save", "tr_file", "Save the transforms as", "transformation.csv" ),
                actionButton( "tr_save", "Save CSV" )
            ),
            plotOutput( "tr_plot", height = "420px" ),
            textOutput( "tr_note" ),
            DTOutput( "tr_table" )
        )
    ),

    nav_panel(
        "Gates",
        layout_columns(
            col_widths = c( 4, 5, 3 ),
            card(
                card_header( "Gate hierarchy" ),
                div( class = "d-flex flex-wrap gap-1",
                     actionButton( "gate_new", "New gate", class = "btn-sm" ),
                     actionButton( "gate_delete", "Delete", class = "btn-sm" ),
                     actionButton( "gate_recalc", "Recalculate", class = "btn-sm" ),
                     actionButton( "gate_recalc_all", "Recalculate all", class = "btn-sm" ),
                     actionButton( "undo", "Undo", class = "btn-sm" ) ),
                DTOutput( "gate_table" )
            ),
            card(
                card_header( "Training events" ),
                plotOutput( "gate_plot", height = "560px", click = "gate_click" ),
                textOutput( "click_note" )
            ),
            card(
                card_header( "Edit gate" ),
                uiOutput( "editor" ),
                div( class = "d-flex flex-wrap gap-1",
                     actionButton( "gate_apply", "Apply", class = "btn-primary" ),
                     actionButton( "gate_recommend", "Recommend algorithm" ) )
            )
        )
    ),

    nav_panel(
        "Save and run",
        layout_columns(
            col_widths = c( 6, 6 ),
            card(
                card_header( "Save" ),
                path.input( "save", "model_file", "Gating model (.agmodel)", "model.agmodel" ),
                actionButton( "save_model", "Save model" ),
                hr(),
                path.input( "dir", "training_figure_dir", "Folder of training figures",
                           file.path( start.value( "figure.dir", "figures" ), "training" ) ),
                actionButton( "save_figures", "Write training figures" )
            ),
            card(
                card_header( "Run on every sample" ),
                path.input( "dir", "statistics_dir", "Statistics folder", start.value( "statistics.dir", "statistics" ) ),
                checkboxInput( "sample_figures", "Figures of every gate on every sample", FALSE ),
                path.input( "dir", "sample_figure_dir", "Folder of sample figures", start.value( "figure.dir", "figures" ) ),
                actionButton( "run_all", "Run", class = "btn-primary" ),
                DTOutput( "run_status" )
            )
        ),
        card(
            card_header( "Finish" ),
            p( "Close the app and return the model and the experiment (with any transform changes) to R." ),
            actionButton( "finish", "Finish", class = "btn-success" )
        )
    )
)


# ==============================================================================
# Server
# ==============================================================================

server <- function( input, output, session ) {

    rv <- reactiveValues( fcs.setup = NULL, ap = NULL, ts = NULL, model = NULL, result = NULL,
                          history = list(), selected = NULL, log = character( 0 ), run = NULL )
    form <- reactiveVal( NULL )

    note <- function( text ) rv$log <- c( rv$log, text )

    run.step <- function( expr ) {
        withCallingHandlers(
            tryCatch( { expr; TRUE }, error = function( e ) {
                if ( inherits( e, "shiny.silent.error" ) ) return( FALSE )
                showNotification( conditionMessage( e ), type = "error", duration = 10 )
                note( paste( "Error:", conditionMessage( e ) ) )
                FALSE
            } ),
            message = function( m ) {
                note( trimws( conditionMessage( m ) ) )
                invokeRestart( "muffleMessage" )
            } )
    }

    path.of <- function( x ) {
        x <- trimws( x %||% "" )
        if ( nzchar( x ) ) path.expand( x ) else NULL
    }

    existing <- function( x, what ) {
        f <- path.of( x )
        if ( is.null( f ) ) stop( sprintf( "%s: enter a path.", what ) )
        if ( !file.exists( f ) ) stop( sprintf( "%s not found: %s", what, f ) )
        f
    }

    roots <- c( Home = normalizePath( "~" ), `Working folder` = normalizePath( getwd() ), getVolumes()() )
    path.kinds <- c( data_dir = "dir", panel_file = "file", experiment_file = "file", compensation_file = "file",
                     transformation_file = "file", agmodel_file = "file", flowjo_file = "file",
                     population_gates_file = "file", activation_gates_file = "file", tr_file = "save",
                     model_file = "save", training_figure_dir = "dir", statistics_dir = "dir",
                     sample_figure_dir = "dir" )
    for ( id in names( path.kinds ) ) local( {
        field <- id
        kind <- path.kinds[[ id ]]
        button <- paste0( field, "_browse" )
        switch( kind,
                dir  = shinyDirChoose( input, button, roots = roots ),
                file = shinyFileChoose( input, button, roots = roots, filetypes = path.types[[ field ]] ),
                save = shinyFileSave( input, button, roots = roots, filetypes = path.types[[ field ]] ) )
        observeEvent( input[[ button ]], {
            chosen <- input[[ button ]]
            path <- switch( kind,
                            dir  = parseDirPath( roots, chosen ),
                            file = parseFilePaths( roots, chosen )$datapath,
                            save = parseSavePath( roots, chosen )$datapath )
            if ( length( path ) == 1L && nzchar( path ) )
                updateTextInput( session, field, value = as.character( path ) )
        } )
    } )

    defs <- function() rv$model$gate.definitions %||% list()

    recompute <- function( model, history = TRUE ) {
        result <- withProgress( message = "Calculating gates", value = 0.5,
            train.boundaries( model$gate.definitions, rv$ts$data, done = model$trained.boundaries,
                              transforms = model$transforms ) )
        model <- record.training( model, result, rv$ts )
        if ( history && !is.null( rv$model ) )
            rv$history <- c( utils::tail( rv$history, 19L ), list( rv$model ) )
        rv$model <- model
        rv$result <- result
        failed <- vapply( result$status, function( s ) isTRUE( s$flag %in% c( "failed", "skipped" ) ), logical( 1 ) )
        for ( g in names( result$status )[ failed ] )
            note( sprintf( "Gate '%s': %s", g, result$status[[ g ]]$message ) )
    }

    set.model <- function( model, history = TRUE ) {
        prepared <- prepare.model( model, rv$fcs.setup )
        for ( msg in prepared$messages ) note( msg )
        recompute( prepared$model, history )
    }

    install.data <- function( fcs.setup, ap, model ) {
        rv$fcs.setup <- fcs.setup
        rv$ap <- ap
        ts <- training.set( fcs.setup, ap )
        for ( msg in ts$messages ) note( msg )
        rv$ts <- ts
        note( sprintf( "%d training events from %d sample(s).", length( ts$rows ), length( ts$samples ) ) )
        set.model( model %||% new.agmodel( list() ), history = FALSE )
    }

    if ( !is.null( start$fcs.setup ) )
        run.step( {
            model <- start$model
            if ( is.null( model ) && ( !is.null( start.ap$model.filename ) ||
                                       !is.null( start.ap$population.gates.filename ) ) )
                model <- model.from.files( start.ap, start$fcs.setup$fcs.panel )
            isolate( install.data( start$fcs.setup, start.ap, model ) )
        } )

    channel.choices <- reactive( {
        req( rv$fcs.setup )
        ch <- colnames( rv$ts$data )
        labels <- rv$fcs.setup$fcs.marker.label.figure
        stats::setNames( ch, ifelse( ch %in% names( labels ), labels[ ch ], ch ) )
    } )

    axes <- reactive( {
        req( rv$fcs.setup, rv$ts )
        figure.axes( rv$fcs.setup, rv$ap, rv$ts$data )
    } )

    population.rows <- function( key ) {
        if ( !has.value( key ) ) return( seq_len( nrow( rv$ts$data ) ) )
        k <- parse.population.key( key )
        rv$result$events[[ k[ 1L ] ]][[ k[ 2L ] ]] %||% integer( 0 )
    }


    # --------------------------------------------------------------------------
    # Data
    # --------------------------------------------------------------------------

    observeEvent( input$load_data, run.step( {
        optional <- function( x, what ) if ( is.null( path.of( x ) ) ) NULL else existing( x, what )
        ap <- get.autogating.param(
            input$cytometer, data.dir = existing( input$data_dir, "Folder of FCS files" ),
            panel.filename = existing( input$panel_file, "Panel CSV" ),
            experiment.filename = existing( input$experiment_file, "Experiment CSV" ),
            compensation.filename = optional( input$compensation_file, "Compensation CSV" ),
            transformation.filename = optional( input$transformation_file, "Transformation CSV" ),
            events.per.sample = input$events_per_sample,
            statistics.dir = input$statistics_dir )
        fcs.setup <- withProgress( message = "Reading FCS files", value = 0.5, setup.autogating( ap ) )
        rv$history <- list()
        install.data( fcs.setup, ap, rv$model )
    } ) )

    observeEvent( input$resample, run.step( {
        if ( is.null( rv$fcs.setup ) ) stop( "Load an experiment first." )
        ap <- rv$ap
        ap$events.per.sample <- as.integer( input$events_per_sample )
        rv$ap <- ap
        ts <- training.set( rv$fcs.setup, ap )
        rv$ts <- ts
        note( sprintf( "%d training events from %d sample(s).", length( ts$rows ), length( ts$samples ) ) )
        set.model( rv$model )
    } ) )

    observeEvent( input$load_model, run.step( {
        if ( is.null( rv$fcs.setup ) ) stop( "Load an experiment first." )
        panel <- rv$fcs.setup$fcs.panel
        model <- switch( input$model_source,
            empty = new.agmodel( list() ),
            agmodel = read.agmodel( existing( input$agmodel_file, ".agmodel file" ), panel ),
            flowjo = {
                imported <- import.flowjo.workspace( existing( input$flowjo_file, "FlowJo workspace" ),
                                                     rv$fcs.setup, rv$ap )
                gate.defs <- imported$gate.definitions
                if ( isTRUE( input$flowjo_convert ) ) {
                    transforms <- rv$fcs.setup$ap$transform.params
                    converted <- withProgress( message = "Converting imported gates", value = 0.5,
                        convert.imported.gates( gate.defs, rv$ts$data, transforms[ lengths( transforms ) > 0L ], panel ) )
                    gate.defs <- converted$gate.definitions
                    note( sprintf( "%d imported gate(s) converted to calculated gates.", nrow( converted$proposals ) ) )
                }
                new.agmodel( gate.defs )
            },
            rlist = {
                activation <- path.of( input$activation_gates_file )
                up <- convert.gate.definitions(
                    load.gate.definitions( existing( input$population_gates_file, "Population gates" ) ),
                    load.gate.definitions( if ( is.null( activation ) ) NULL else existing( activation, "Activation gates" ) ) )
                for ( msg in up$messages ) note( msg )
                new.agmodel( up$gate.definitions )
            } )
        rv$history <- list()
        rv$selected <- NULL
        form( NULL )
        set.model( model, history = FALSE )
        note( sprintf( "Model loaded: %d gate(s).", length( model$gate.definitions ) ) )
    } ) )

    output$data_summary <- renderText( {
        if ( is.null( rv$fcs.setup ) ) return( "No experiment loaded." )
        fs <- rv$fcs.setup
        trained <- length( rv$model$trained.boundaries )
        sprintf( paste( "Samples: %d (%d for training)", "Events: %s", "Training events: %s",
                        "Channels: %d", "Gates: %d (%d with a boundary)", sep = "\n" ),
                 fs$fcs.sample.n, length( rv$ts$samples ), format( fs$fcs.event.n, big.mark = "," ),
                 format( length( rv$ts$rows ), big.mark = "," ), ncol( fs$fcs.expr.data ),
                 length( defs() ), trained )
    } )

    output$log <- renderText( paste( utils::tail( rv$log, 200L ), collapse = "\n" ) )


    # --------------------------------------------------------------------------
    # Transforms
    # --------------------------------------------------------------------------

    observeEvent( rv$fcs.setup, {
        params <- rv$fcs.setup$ap$transform.params
        ch <- channel.choices()
        ch <- ch[ ch %in% names( params )[ lengths( params ) > 0L ] ]
        keep <- if ( isTRUE( input$tr_channel %in% ch ) ) input$tr_channel else ch[ 1L ]
        updateSelectInput( session, "tr_channel", choices = ch, selected = keep )
    } )

    observeEvent( rv$model, {
        keep <- input$tr_population %||% ""
        choices <- c( "Training events" = "", population.keys( defs() ) )
        updateSelectInput( session, "tr_population", choices = choices,
                           selected = if ( keep %in% choices ) keep else "" )
    } )

    fill.transform <- function() {
        p <- rv$fcs.setup$ap$transform.params[[ input$tr_channel ]]
        req( length( p ) > 0L )
        updateRadioButtons( session, "tr_id", selected = as.character( p$id ) )
        updateNumericInput( session, "tr_t", value = p$scale.t )
        updateNumericInput( session, "tr_w", value = p$logicle.w %||% 0.5 )
        updateNumericInput( session, "tr_m", value = p$logicle.m %||% 4.5 )
        updateNumericInput( session, "tr_logm", value = p$log.m %||% 4.5 )
        updateNumericInput( session, "tr_a", value = p$linear.a %||% 0 )
    }
    observeEvent( input$tr_channel, { req( input$tr_channel ); fill.transform() } )
    observeEvent( input$tr_reset, { req( input$tr_channel ); fill.transform() } )

    candidate.transform <- reactive( {
        req( input$tr_id, is.finite( input$tr_t ), input$tr_t > 0 )
        id <- as.integer( input$tr_id )
        if ( id == 1L ) {
            req( is.finite( input$tr_w ), is.finite( input$tr_m ), input$tr_w >= 0, input$tr_m > 2 * input$tr_w )
            list( id = 1L, scale.t = input$tr_t, logicle.w = input$tr_w, logicle.m = input$tr_m, logicle.a = 0 )
        } else if ( id == 2L ) {
            req( is.finite( input$tr_logm ), input$tr_logm > 0 )
            list( id = 2L, scale.t = input$tr_t, log.m = input$tr_logm )
        } else {
            list( id = 0L, scale.t = input$tr_t, linear.a = if ( is.finite( input$tr_a ) ) input$tr_a else 0 )
        }
    } )

    transform.preview <- reactive( {
        req( rv$ts, input$tr_channel %in% colnames( rv$ts$data ) )
        ch <- input$tr_channel
        p <- candidate.transform()
        tr <- make.transform( p )
        old.inverse <- rv$fcs.setup$ap$transform.inv[[ ch ]]
        raw <- old.inverse( rv$ts$data[ population.rows( input$tr_population ), ch ] )
        tmp <- rv$fcs.setup
        tmp$ap$transform[[ ch ]] <- tr$forward
        tmp$ap$transform.inv[[ ch ]] <- tr$inverse
        tmp$ap$transform.params[[ ch ]] <- p
        ax <- figure.axes( tmp, rv$ap, rv$ts$data )[[ ch ]]
        list( x = tr$forward( raw ), raw = raw, zero = tr$forward( 0 ), axis = ax )
    } )

    output$tr_plot <- renderPlot( {
        pv <- transform.preview()
        lim <- pv$axis$limits
        x <- pv$x
        x[ !is.finite( x ) ] <- lim[ 1 ]
        x <- pmin( pmax( x, lim[ 1 ] ), lim[ 2 ] )
        ggplot( data.frame( x = x ), aes( x = x ) ) +
            geom_histogram( bins = 256L, boundary = lim[ 1 ], fill = "grey45", colour = NA ) +
            geom_vline( xintercept = pv$zero, linetype = "dashed", colour = "#C0392B" ) +
            scale_x_continuous( pv$axis$title, breaks = pv$axis$breaks, labels = pv$axis$labels,
                                limits = lim, expand = expansion( mult = 0.01 ) ) +
            labs( y = "Events" ) +
            theme_bw( base_size = 13 ) + theme( panel.grid = element_blank() )
    }, res = 96 )

    output$tr_note <- renderText( {
        pv <- transform.preview()
        lim <- pv$axis$limits
        sprintf( "%s events; %.1f%% below zero (dashed line); %.1f%% at or below the bottom of the axis.",
                 format( length( pv$raw ), big.mark = "," ), 100 * mean( pv$raw < 0 ),
                 100 * mean( !is.finite( pv$x ) | pv$x <= lim[ 1 ] ) )
    } )

    output$tr_table <- renderDT( {
        req( rv$fcs.setup )
        params <- rv$fcs.setup$ap$transform.params
        params <- params[ lengths( params ) > 0L ]
        cell <- function( p, key ) { v <- p[[ key ]]; if ( is.null( v ) ) NA_real_ else as.numeric( v ) }
        tbl <- data.frame( channel = names( params ),
                           transform = c( "0" = "linear", "1" = "logicle", "2" = "log" )[ as.character( vapply( params, `[[`, numeric( 1 ), "id" ) ) ],
                           T = vapply( params, cell, numeric( 1 ), "scale.t" ),
                           W = vapply( params, cell, numeric( 1 ), "logicle.w" ),
                           M = vapply( params, function( p ) cell( p, if ( identical( as.integer( p$id ), 2L ) ) "log.m" else "logicle.m" ), numeric( 1 ) ),
                           A = vapply( params, cell, numeric( 1 ), "linear.a" ),
                           stringsAsFactors = FALSE )
        datatable( tbl, rownames = FALSE, selection = "none", class = "compact stripe",
                   options = list( dom = "t", paging = FALSE ) ) |>
            formatRound( c( "W", "M" ), 3 ) |>
            formatRound( c( "T", "A" ), 0 )
    } )

    observeEvent( input$tr_apply, run.step( {
        ch <- input$tr_channel
        p <- candidate.transform()
        fs <- retransform.channel( rv$fcs.setup, ch, p )
        ts <- rv$ts
        ts$data <- fs$fcs.expr.data[ ts$rows, , drop = FALSE ]
        rv$fcs.setup <- fs
        rv$ts <- ts
        set.model( rv$model )
        note( sprintf( "Transform of %s changed; boundaries on it were converted.", ch ) )
    } ) )

    observeEvent( input$tr_save, run.step( {
        file <- path.of( input$tr_file )
        write.transformation.csv( rv$fcs.setup, file )
        note( sprintf( "Transforms saved to %s.", file ) )
    } ) )


    # --------------------------------------------------------------------------
    # Gate hierarchy
    # --------------------------------------------------------------------------

    gate.rows <- reactive( {
        req( rv$model )
        gate.table( rv$model, rv$result )
    } )

    output$gate_table <- renderDT( {
        tbl <- gate.rows()
        sel <- which( tbl$gate == ( isolate( rv$selected ) %||% "" ) )
        depth <- integer( 0 )
        parent.gate <- vapply( strsplit( tbl$parent, "/", fixed = TRUE ), function( k ) if ( length( k ) ) k[ 1L ] else "", "" )
        for ( i in seq_len( nrow( tbl ) ) )
            depth[ tbl$gate[ i ] ] <- if ( parent.gate[ i ] %in% names( depth ) ) depth[[ parent.gate[ i ] ]] + 1L else 0L
        shown <- data.frame( gate = sprintf( "<span style=\"padding-left: %.1fem; white-space: nowrap\">%s</span>",
                                             0.8 * depth[ tbl$gate ], htmltools::htmlEscape( tbl$gate ) ),
                             populations = tbl$populations, algorithm = tbl$algorithm, status = tbl$status,
                             stringsAsFactors = FALSE )
        datatable( shown, rownames = FALSE, class = "compact stripe", escape = FALSE,
                   selection = list( mode = "single", selected = sel ),
                   options = list( dom = "ft", paging = FALSE, scrollY = "600px", scrollCollapse = TRUE,
                                   ordering = FALSE ) ) |>
            formatStyle( names( shown ), fontSize = "0.85rem" )
    } )

    edit.gate <- function( name ) {
        g <- gates.by.name( defs() )[[ name ]]
        req( g )
        mode <- rv$model$application$mode.per.gate[[ name ]] %||% "fixed"
        form( c( gate.form.values( g, mode ), old.name = name ) )
    }

    observeEvent( input$gate_table_rows_selected, {
        tbl <- gate.rows()
        name <- tbl$gate[ input$gate_table_rows_selected ]
        if ( length( name ) != 1L ) return()
        if ( !identical( name, rv$selected ) || is.null( form() ) || !identical( form()$old.name, name ) ) {
            rv$selected <- name
            edit.gate( name )
        }
    } )

    observeEvent( input$gate_new, {
        req( rv$model )
        v <- gate.form.values( NULL )
        g <- gates.by.name( defs() )[[ rv$selected %||% "" ]]
        if ( !is.null( g ) ) {
            pops <- names( effective.gate.def( g, gates.by.name( defs() ) )$populations )
            if ( length( pops ) ) v$parent <- population.key( g$gate.name, pops[ length( pops ) ] )
        }
        form( c( v, list( old.name = NULL ) ) )
    } )

    observeEvent( input$gate_delete, {
        req( rv$selected )
        d <- drop.gate( rv$model, rv$selected )
        showModal( modalDialog(
            title = "Delete gate",
            p( sprintf( "Delete %s?", paste( d$removed, collapse = ", " ) ) ),
            if ( length( d$broken ) )
                p( sprintf( "These gates use a deleted gate's threshold or population and will not run: %s.",
                            paste( d$broken, collapse = ", " ) ) ),
            footer = tagList( modalButton( "Cancel" ), actionButton( "confirm_delete", "Delete", class = "btn-danger" ) ) ) )
    } )

    observeEvent( input$confirm_delete, run.step( {
        removeModal()
        d <- drop.gate( rv$model, rv$selected )
        recompute( d$model )
        rv$selected <- NULL
        form( NULL )
        note( sprintf( "Deleted %s.", paste( d$removed, collapse = ", " ) ) )
    } ) )

    observeEvent( input$gate_recalc, run.step( {
        req( rv$selected )
        model <- rv$model
        model$trained.boundaries[[ rv$selected ]] <- NULL
        recompute( model )
    } ) )

    observeEvent( input$gate_recalc_all, run.step( {
        req( rv$model )
        model <- rv$model
        model$trained.boundaries <- list()
        recompute( model )
    } ) )

    observeEvent( input$undo, run.step( {
        if ( length( rv$history ) == 0L ) stop( "Nothing to undo." )
        previous <- rv$history[[ length( rv$history ) ]]
        rv$history <- rv$history[ -length( rv$history ) ]
        set.model( previous, history = FALSE )
        if ( !is.null( rv$selected ) && !( rv$selected %in% names( gates.by.name( defs() ) ) ) )
            rv$selected <- NULL
        if ( is.null( rv$selected ) ) form( NULL ) else edit.gate( rv$selected )
    } ) )


    # --------------------------------------------------------------------------
    # Gate figure
    # --------------------------------------------------------------------------

    output$gate_plot <- renderPlot( {
        req( rv$model, rv$selected )
        by.name <- gates.by.name( defs() )
        g <- by.name[[ rv$selected ]]
        req( g )
        fig <- gate.figure( effective.gate.def( g, by.name ), rv$result, rv$ts$data, axes(), rv$ap )
        validate( need( !is.null( fig ), sprintf( "Gate '%s' has no boundary on the training events: %s",
                                                  rv$selected, rv$result$status[[ rv$selected ]]$message %||% "" ) ) )
        fig
    }, res = 96 )

    output$click_note <- renderText( {
        f <- form()
        if ( is.null( f ) || is.null( f$old.name ) || !( f$type %in% ag.threshold.types ) ) return( "" )
        "Click the figure to set a manual threshold, then Apply; clear the field to calculate it again."
    } )

    observeEvent( input$gate_click, {
        f <- form()
        click <- input$gate_click
        if ( is.null( f ) || is.null( f$old.name ) || !identical( f$old.name, rv$selected ) ) return()
        by.name <- gates.by.name( defs() )
        eff <- effective.gate.def( by.name[[ f$old.name ]], by.name )
        if ( identical( eff$gate.type, "1dsep" ) ) {
            v <- if ( separation.axis( eff ) == "y" ) click$y else click$x
            if ( is.null( v ) ) return()
            updateNumericInput( session, "ed_fixed_threshold_x", value = signif( v, 5 ) )
        } else if ( identical( eff$gate.type, "2dsep" ) ) {
            updateNumericInput( session, "ed_fixed_threshold_x", value = signif( click$x, 5 ) )
            updateNumericInput( session, "ed_fixed_threshold_y", value = signif( click$y, 5 ) )
        }
    } )


    # --------------------------------------------------------------------------
    # Gate editor
    # --------------------------------------------------------------------------

    collect.form <- function() {
        f <- form()
        val <- function( id, default ) { x <- input[[ id ]]; if ( is.null( x ) ) default else x }
        type <- val( "ed_type", f$type )
        regions <- form.regions( type )
        v <- list(
            name = val( "ed_name", f$name ), type = type, parent = val( "ed_parent", f$parent ),
            x = val( "ed_x", f$x ), y = val( "ed_y", f$y ), sep.axis = val( "ed_sep", f$sep.axis ),
            algorithm = val( "ed_alg", f$algorithm ), algorithm.x = val( "ed_alg_x", f$algorithm.x ),
            algorithm.y = val( "ed_alg_y", f$algorithm.y ),
            populations = stats::setNames( vapply( regions, function( r )
                as.character( val( paste0( "ed_pop_", region.id( r ) ), f$populations[ r ] %||% "" ) ),
                character( 1 ) ), regions ),
            width.sd = val( "ed_width", f$width.sd ), debris.algorithm = val( "ed_debris", f$debris.algorithm ),
            edge.density = val( "ed_edge", f$edge.density ), origin = val( "ed_origin", f$origin ),
            mode = val( "ed_mode", f$mode ) )
        for ( axis in c( "x", "y" ) ) {
            for ( key in c( "axis.lo", "axis.hi", "tail.fraction", "threshold.factor", "fixed.threshold" ) ) {
                id <- paste0( "ed_", sub( ".", "_", key, fixed = TRUE ), "_", axis )
                v[[ paste0( key, ".", axis ) ]] <- val( id, f[[ paste0( key, ".", axis ) ]] )
            }
            v[[ paste0( "replicate.", axis ) ]] <- val( paste0( "ed_rep_", axis ), f[[ paste0( "replicate.", axis ) ]] )
            v[[ paste0( "replicate.algorithm.", axis ) ]] <-
                val( paste0( "ed_repalg_", axis ), f[[ paste0( "replicate.algorithm.", axis ) ]] )
        }
        c( v, list( old.name = f$old.name ) )
    }

    observeEvent( input$ed_type, {
        f <- form()
        req( f )
        if ( identical( input$ed_type, f$type ) ) return()
        v <- collect.form()
        type <- input$ed_type
        fresh <- gate.form.values( list( gate.type = type, populations = default.populations.for.type( type ) ) )
        v$type <- type
        v$populations <- fresh$populations
        v$algorithm <- v$algorithm.x <- v$algorithm.y <- default.algorithm[[ type ]]
        form( v )
    }, ignoreInit = TRUE )

    threshold.fields <- function( f, axis, title ) {
        ids <- function( key ) paste0( "ed_", key, "_", axis )
        value <- function( key ) f[[ paste0( key, ".", axis ) ]]
        gates <- vapply( Filter( function( g ) g$gate.type %in% ag.threshold.types && !identical( g$gate.name, f$old.name ),
                                 defs() ), function( g ) g$gate.name, character( 1 ) )
        sources <- c( "Choose" = "",
                      stats::setNames( paste0( "gate:", gates ), paste( "Gate", gates ) ),
                      stats::setNames( paste0( "population:", population.keys( defs() ) ),
                                       paste( "Population", population.keys( defs() ) ) ) )
        tags$details(
            tags$summary( title ),
            numericInput( ids( "fixed_threshold" ), "Manual threshold (display units)", value( "fixed.threshold" ), step = 0.001 ),
            numericInput( ids( "axis_lo" ), "Search from (display units)", value( "axis.lo" ), min = 0, max = 1, step = 0.01 ),
            numericInput( ids( "axis_hi" ), "Search to (display units)", value( "axis.hi" ), min = 0, max = 1, step = 0.01 ),
            numericInput( ids( "tail_fraction" ), "Tail fraction (tail algorithm)", value( "tail.fraction" ), min = 0, step = 1e-4 ),
            numericInput( ids( "threshold_factor" ), "Threshold factor (raw units)", value( "threshold.factor" ), min = 0, step = 0.05 ),
            selectInput( ids( "rep" ), "Replicate from (replicate algorithm)", sources, value( "replicate" ) ),
            selectInput( ids( "repalg" ), "Algorithm on the population", setdiff( ag.threshold.algorithms, "replicate" ),
                         value( "replicate.algorithm" ) ) )
    }

    output$editor <- renderUI( {
        f <- form()
        if ( is.null( f ) ) return( p( "Select a gate, or add a new one." ) )
        all.defs <- defs()
        excluded <- if ( is.null( f$old.name ) ) character( 0 ) else drop.gate( rv$model, f$old.name )$removed
        others <- Filter( function( g ) !( g$gate.name %in% excluded ), all.defs )
        parents <- c( "All events" = "", population.keys( others ) )
        chans <- channel.choices()
        type <- f$type
        algorithms <- function( current ) {
            a <- ag.algorithms.by.type[[ type ]]
            if ( !identical( current, "imported" ) ) a <- setdiff( a, "imported" )
            a
        }
        pop.inputs <- lapply( form.regions( type ), function( r )
            textInput( paste0( "ed_pop_", region.id( r ) ), region.labels[[ r ]], f$populations[[ r ]] %||% "" ) )
        tagList(
            tags$b( if ( is.null( f$old.name ) ) "New gate" else f$old.name ),
            textInput( "ed_name", "Name", f$name ),
            selectInput( "ed_type", "Type", ag.gate.types, type ),
            selectInput( "ed_parent", "Parent population", parents, f$parent ),
            if ( type == "replicate" ) selectInput( "ed_origin", "Copy of gate",
                setdiff( vapply( others, function( g ) g$gate.name, character( 1 ) ), f$old.name ), f$origin )
            else tagList(
                selectInput( "ed_x", if ( type == "1dsep" ) "Threshold channel" else "X channel", chans, f$x ),
                selectInput( "ed_y", if ( type == "1dsep" ) "Shown against (optional)" else "Y channel",
                             if ( type %in% c( "1dsep", "free" ) ) c( "None" = "", chans ) else chans, f$y ),
                if ( type == "1dsep" ) radioButtons( "ed_sep", "Threshold channel on the", c( "x axis" = "x", "y axis" = "y" ),
                                                     f$sep.axis, inline = TRUE ),
                if ( type == "2dsep" ) tagList(
                    selectInput( "ed_alg_x", "Algorithm, X", algorithms( f$algorithm.x ), f$algorithm.x ),
                    selectInput( "ed_alg_y", "Algorithm, Y", algorithms( f$algorithm.y ), f$algorithm.y ) )
                else selectInput( "ed_alg", "Algorithm", algorithms( f$algorithm ), f$algorithm ),
                tags$label( "Populations" ),
                if ( type == "2dsep" ) p( class = "small text-muted",
                    "Give quadrants the same name to join them; leave a name empty to leave the quadrant out." ),
                pop.inputs,
                if ( type == "1dsep" ) threshold.fields( f, "x", "Threshold settings" ),
                if ( type == "2dsep" ) tagList( threshold.fields( f, "x", "X threshold settings" ),
                                                threshold.fields( f, "y", "Y threshold settings" ) ),
                if ( type == "singlets" ) numericInput( "ed_width", "Band half-width (robust SDs)", f$width.sd, min = 0.5, step = 0.5 ),
                if ( type == "free" ) tagList(
                    selectInput( "ed_debris", "Debris cut algorithm (cells)",
                                 c( "Default" = "", setdiff( ag.threshold.algorithms, "replicate" ) ), f$debris.algorithm ),
                    numericInput( "ed_edge", "Edge density (cells)", f$edge.density, min = 0, step = 0.001 ) ) ),
            radioButtons( "ed_mode", "On each sample", c( "Fixed boundary" = "fixed", "Recalculate" = "recalculate" ),
                          f$mode, inline = TRUE )
        )
    } )

    observeEvent( input$gate_apply, run.step( {
        req( form() )
        v <- collect.form()
        old <- if ( is.null( v$old.name ) ) NULL else gates.by.name( defs() )[[ v$old.name ]]
        g <- gate.from.form( v, old )
        model <- put.gate( rv$model, g, v$old.name, mode = v$mode )
        recompute( model )
        rv$selected <- g$gate.name
        edit.gate( g$gate.name )
        note( sprintf( "Gate '%s' %s.", g$gate.name, if ( is.null( old ) ) "added" else "updated" ) )
    } ) )

    observeEvent( input$gate_recommend, run.step( {
        req( form() )
        v <- collect.form()
        rows <- population.rows( v$parent )
        data <- rv$ts$data[ rows, , drop = FALSE ]
        suggest <- function( ch ) recommend.algorithm( list( gate.type = "1dsep", gate.marker.x = ch ), data )
        if ( v$type == "1dsep" ) {
            a <- suggest( v$x )
            updateSelectInput( session, "ed_alg", selected = a )
            note( sprintf( "Recommended for %s: %s.", v$x, a ) )
        } else if ( v$type == "2dsep" ) {
            ax <- suggest( v$x )
            ay <- suggest( v$y )
            updateSelectInput( session, "ed_alg_x", selected = ax )
            updateSelectInput( session, "ed_alg_y", selected = ay )
            note( sprintf( "Recommended: %s on %s, %s on %s.", ax, v$x, ay, v$y ) )
        } else {
            note( "Algorithms are recommended for threshold gates only." )
        }
    } ) )


    # --------------------------------------------------------------------------
    # Save and run
    # --------------------------------------------------------------------------

    observeEvent( input$save_model, run.step( {
        file <- path.of( input$model_file )
        write.agmodel( rv$model, file, rv$fcs.setup$fcs.panel )
        note( sprintf( "Model saved to %s.", file ) )
    } ) )

    observeEvent( input$save_figures, run.step( {
        dir <- path.of( input$training_figure_dir )
        files <- withProgress( message = "Writing figures", value = 0.5,
            write.gate.figures( defs(), rv$result, rv$ts$data, rv$fcs.setup, rv$ap, dir, title = "Training" ) )
        note( sprintf( "%d figure(s) written to %s.", length( files ), dir ) )
    } ) )

    observeEvent( input$run_all, run.step( {
        ap <- rv$ap
        ap$statistics.dir <- path.of( input$statistics_dir )
        ap$figure.samples <- isTRUE( input$sample_figures )
        ap$figure.dir <- path.of( input$sample_figure_dir ) %||% ap$figure.dir
        rv$run <- withProgress( message = "Gating every sample", value = 0.5,
                                run.autogating( rv$fcs.setup, rv$model, ap ) )
    } ) )

    output$run_status <- renderDT( {
        req( rv$run )
        st <- rv$run$status
        flagged <- st[ !is.na( st$flag ), , drop = FALSE ]
        tbl <- if ( nrow( flagged ) ) flagged else
            data.frame( message = sprintf( "%d sample(s) gated; no gate flagged.", length( unique( st$sample ) ) ) )
        datatable( tbl, rownames = FALSE, selection = "none", options = list( dom = "tp", pageLength = 15 ) )
    } )

    observeEvent( input$finish, stopApp( list( model = rv$model, fcs.setup = rv$fcs.setup ) ) )
}

shinyApp( ui, server )
