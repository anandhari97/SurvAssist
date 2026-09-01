library(shiny)
library(survival)
library(survminer)
library(dplyr)
library(readxl)

# -----------------------------
# FUNCTIONS
# -----------------------------

#------------------------------------
#______Summary___________________
#------------------------------------
summary_continuous <- function(df, var){
  x <- df[[var]]
  data.frame(Var=var,
             Mean = round(mean(x, na.rm = TRUE), 2),
             Median = round(median(x, na.rm = TRUE), 2),
             Min = round(min(x, na.rm = TRUE), 2),
             Max = round(max(x, na.rm = TRUE), 2),
             SD = round(sd(x, na.rm = TRUE), 2)
  )
}

summary_categorical <- function(df, var){
  tab <- table(df[[var]])
  prop <- prop.table(tab) * 100
  
  data.frame(Var=var,
             Category = names(tab),
             Frequency = as.numeric(tab),
             Percent = round(as.numeric(prop), 2)
  )
}
#------------------------------------
#______Follow-up%___________________
#------------------------------------
follow_up_calc <- function(df, time, status, event_value, cutoff){
  valid <- !is.na(df[[time]]) & !is.na(df[[status]])
  df <- df[valid, ]
  
  rev_event <- ifelse(df[[status]] == event_value, 0, 1)
  
  fit <- survfit(
    Surv(df[[time]], rev_event) ~ 1
  )
  
  s <- summary(fit, times = cutoff, extend = TRUE)
  fup<- round(s$surv * 100, 2)
  data.frame(
    Cutoff = cutoff,
    Follow_up_percent = fup
  )
}
#------------------------------------
#______Survival Summary______________
#------------------------------------

Survival_summary <- function(df, status_var, event_value){
  
  n <- nrow(df)
  Event <- sum(df[[status_var]] == event_value, na.rm = TRUE)
  Cens <- sum(df[[status_var]] != event_value, na.rm = TRUE)
  Cens_pct <- (Cens * 100) / n
  
  data.frame(
    n = n,
    Events = Event,
    Censored = Cens,
    Censored_percent = round(Cens_pct, 2)
  )
}

Survival_summary_group <- function(df, status_var, event_value, group_var){
  
  df %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      n = n(),
      Events = sum(.data[[status_var]] == event_value, na.rm = TRUE),
      Censored = sum(.data[[status_var]] != event_value, na.rm = TRUE),
      Censored_percent = round((Censored * 100) / n, 2),
      .groups = "drop"
    )
}

#------------------------------------
#______KM Plot______________
#------------------------------------
get_km_plot <- function(obj, input){
  
  default_labels <- if(!is.null(obj$fit$strata)){
    gsub(".*=", "", names(obj$fit$strata))
  } else NULL
  
  legend_labels <- default_labels
  
  if(input$legend_labels != "" && !is.null(default_labels)){
    user_labels <- trimws(unlist(strsplit(input$legend_labels, ",")))
    if(length(user_labels) == length(default_labels)){
      legend_labels <- user_labels
    }
  }
  
  x_max <- ifelse(is.na(input$xmax), input$cutoff, input$xmax)
  
  legend_title <- if(input$legend_title == ""){
    if(input$use_group) input$group_var else "Overall Survival"
  } else input$legend_title
  
  x_label <- if(input$xlab == "") "Time" else input$xlab
  y_label <- if(input$ylab == "") "Survival Probability" else input$ylab
  
  p <- ggsurvplot(
    obj$fit,
    data = obj$df,
    conf.int = TRUE,
    pval = if(input$use_group) TRUE else FALSE,
    xlim = c(0, x_max),
    break.time.by = input$xbreak,
    xlab = x_label,
    ylab = y_label,
    ggtheme = theme_classic(),
    legend.title = legend_title,
    legend.labs = legend_labels
  )
  
  return(p)
}

#------------------------------------
#______Cox-summary______________
#------------------------------------
cox_summary <- function(fit, df){
  s <- summary(fit)
  
  res <- data.frame(
    Variable = rownames(s$coefficients),
    HR = round(s$coefficients[, "exp(coef)"], 3),
    Lower_95_CI = round(s$conf.int[, "lower .95"], 3),
    Upper_95_CI = round(s$conf.int[, "upper .95"], 3),
    P_value = signif(s$coefficients[, "Pr(>|z|)"], 3),
    stringsAsFactors = FALSE
  )
  
  final <- list()
  
  vars <- all.vars(fit$terms)[-1]  # remove survival part
  
  for (var in vars) {
    
    if (is.factor(df[[var]])) {
      
      levels_var <- levels(df[[var]])
      ref <- levels_var[1]
      
      # reference row
      ref_row <- data.frame(
        Variable = paste0(var, ": ", ref),
        HR = 1,
        Lower_95_CI = "-",
        Upper_95_CI = "-",
        P_value = "-"
      )
      
      # model rows for this variable
      var_rows <- res[grepl(paste0("^", var), res$Variable), ]
      
      # clean labels
      var_rows$Variable <- gsub(var, paste0(var, ": "), var_rows$Variable)
      
      final[[var]] <- rbind(ref_row, var_rows)
      
    } else {
      
      final[[var]] <- res[res$Variable == var, ]
    }
  }
  
  final_df <- do.call(rbind, final)
  rownames(final_df) <- NULL
  
  final_df
}

#------------------------------------
#______Cox-PH______________
#------------------------------------
get_cox_results <- function(df, input){
  
  df$Cox_time <- as.numeric(as.character(df[[input$time]]))
  status_vec <- as.character(df[[input$status]])
  df$Cox_event <- as.numeric(status_vec) == input$event_value
  
  df <- df[!is.na(df$Cox_time) & !is.na(df$Cox_event), ]
  
  df <- df %>%
    select(all_of(c("Cox_time","Cox_event", input$cox_vars))) %>%
    na.omit()
  
  if(!is.null(input$cox_cat_vars)){
    for(v in input$cox_cat_vars){
      if(v %in% names(df)){
        df[[v]] <- as.factor(df[[v]])
      }
    }
  }
  
  # ---- UNIVARIABLE ----
  if(input$cox_type == "Univariable"){
    
    results <- lapply(input$cox_vars, function(v){
      fmla <- as.formula(paste0("Surv(Cox_time, Cox_event) ~ ", v))
      fit <- coxph(fmla, data = df)
      cox_summary(fit, df)
    })
    
    res <- do.call(rbind, results)
    
  } else {
    
    # ---- MULTIVARIABLE ----
    fmla <- as.formula(
      paste0("Surv(Cox_time, Cox_event) ~ ",
             paste(input$cox_vars, collapse = " + "))
    )
    
    fit <- coxph(fmla, data = df)
    res <- cox_summary(fit, df)
  }
  
  res
}

# -----------------------------
# UI
# -----------------------------

ui <- fluidPage(
  titlePanel("SurvAssist"),
  hr(),
  h4("About"),
  helpText("An Interactive Survival Analysis Platform for Clinicians."),
  sidebarLayout(
    sidebarPanel(
      
      fileInput("file", "Upload Data (.csv or .xlsx)"),
      uiOutput("sheet_ui"),
      uiOutput("var_select"),
      
      checkboxInput("use_group", "Enable Group Comparison", value = FALSE),
      
      conditionalPanel(
        condition = "input.use_group == true",
        selectInput("group_var", "Group Variable", choices = NULL)
      ),
      
      numericInput("event_value", "Event Value", value = 1),
      numericInput("cutoff", "Time Point", value = NULL),
      hr(),
      h4("KM Plot Customization"),
      
      textInput("legend_title", "Legend Title", ""),
      textInput("legend_labels", "Legend Labels (comma separated)", ""),
      textInput("xlab", "X-axis Label", "Time"),
      textInput("ylab", "Y-axis Label", "Survival Probability"),
      
      numericInput("xmax", "X-axis Max", value = NA),
      numericInput("xbreak", "X-axis Break Interval", value = 6),
      checkboxGroupInput(
        "analysis_type",
        "Select Analysis",
        choices = c("Summary Statistics", "Follow-up%", "Kaplan-Meier","Cox Regression")
      ),
      
      hr(),
      
      h4("Summary Options"),
      selectInput("summary_var", "Variable", choices = NULL),
      radioButtons("var_type", "Variable Type",
                   choices = c("Continuous", "Categorical")),
      hr(),
      selectInput("cox_vars", "Cox Variables", 
                  choices = NULL, multiple = TRUE),
      selectInput("cox_cat_vars", "Categorical Variables (for Cox)", 
                  choices = NULL, multiple = TRUE),
      radioButtons("cox_type", "Model Type",
                   choices = c("Univariable", "Multivariable")),
      
      actionButton("run", "Run Analysis")
    ),
    
    mainPanel(
      conditionalPanel(
        condition = "input.analysis_type.includes('Summary Statistics')",
        tableOutput("summary"),
        downloadButton("download_summary","Download Summary"),
        downloadButton("download_all_summary", "Download All Summaries")
      ),
      
      conditionalPanel(
        condition = "input.analysis_type.includes('Follow-up%')",
        tableOutput("followup")
      ),
      
      conditionalPanel(
        condition = "input.analysis_type.includes('Kaplan-Meier')",
        
        h4("Survival Summary"),
        tableOutput("km_summary"),
        
        h4("Survival Probability at Time"),
        tableOutput("km_prob"),
        
        h4("Median Survival"),
        tableOutput("km_median"),
        downloadButton("download_km_summary", "Download KM Summary"),
        downloadButton("download_all_km", "Download All KM Results"),
        conditionalPanel(
          condition = "input.use_group == true",
          h4("Log-rank Test"),
          textOutput("km_pval")
        ),
        
        h4("Kaplan-Meier Curve"),
        plotOutput("km_plot"),
        
        br(),
        downloadButton("download_km", "Download KM Plot")
      ),
      conditionalPanel(
        condition = "input.analysis_type.includes('Cox Regression')",
        
        h4("Cox Regression"),
        
        tableOutput("cox_results"),
        downloadButton("download_cox", "Download Cox Results")
      )
    )
  ),
  hr(),
  
  div(
    style = "text-align:center; font-size:12px; color:gray;",
    "Developed by Anand Hari, Divya Dennis & Dr. Jagathnath Krishna K M | Division of Cancer Epidemiology & Biostatistics, Regional Cancer Centre, Trivandrum,India | Version 1.0 | © 2026"
  )
)

# -----------------------------
# SERVER
# -----------------------------

server <- function(input, output, session){
  
  # SHEETS
  sheets <- reactive({
    req(input$file)
    if (grepl("\\.xlsx$", input$file$name, ignore.case = TRUE)) {
      excel_sheets(input$file$datapath)
    } else NULL
  })
  
  # SHEET UI
  output$sheet_ui <- renderUI({
    req(sheets())
    selectInput("sheet", "Select Sheet", choices = sheets())
  })
  
  data <- reactive({
    req(input$file)
    
    if (grepl("\\.csv$", input$file$name, ignore.case = TRUE)) {
      df <- read.csv(input$file$datapath)
    } else {
      req(input$sheet)
      df <- read_excel(input$file$datapath, sheet = input$sheet)
    }
    
    df <- as.data.frame(df)
    names(df) <- make.names(names(df))
    
    df
  })
  
  output$var_select <- renderUI({
    df <- data()
    
    tagList(
      selectInput("time", "Time Variable", names(df)),
      selectInput("status", "Status Variable", names(df))
    )
  })
  
  observe({
    df <- data()
    updateSelectInput(session, "summary_var", choices = names(df))
    updateSelectInput(session, "group_var", choices = names(df))
    updateSelectInput(session, "cox_vars", choices = names(df))
    updateSelectInput(session, "cox_cat_vars", choices = names(df))
  })
  
  # -----------------------------
  # SUMMARY
  # -----------------------------
  output$summary <- renderTable({
    req(input$run)
    req("Summary Statistics" %in% input$analysis_type)
    
    df <- data()
    x <- df[[input$summary_var]]
    
    if(input$var_type == "Continuous"){
      if(!is.numeric(x)){
        return(data.frame(Error = "Variable is not numeric"))
      }
      summary_continuous(df, input$summary_var)
    } else {
      summary_categorical(df, input$summary_var)
    }
  })
  #-----------------------
  #Download summary
  #-----------------------
  output$download_summary <- downloadHandler(
    filename = function(){
      paste0("summary_", Sys.Date(), ".docx")
    },
    
    content = function(file){
      req(data(), input$summary_var)
      
      library(officer)
      library(flextable)
      df <- data()
      var_name <- as.character(input$summary_var)
      if(input$var_type == "Continuous"){
        res <- summary_continuous(df, var_name)
      } else {
        res <- summary_categorical(df, var_name)
      }
      doc <- read_docx()
      doc <- doc %>%
        body_add_par(paste("Summary:", var_name), style = "heading 1") %>%
        body_add_flextable(flextable(res) %>% autofit())
      print(doc, target = file)
    }
  )
  summary_store <- reactiveValues(data = list())
  observeEvent(input$run, {
    
    req("Summary Statistics" %in% input$analysis_type)
    req(input$summary_var)
    df <- data()
    var_name <-  as.character(input$summary_var)
    
    if(input$var_type == "Continuous"){
      if(!is.numeric(df[[var_name]])) return()
      res <- summary_continuous(df, var_name)
    } else {
      res <- summary_categorical(df, var_name)
    }
    
    key <- paste0(var_name)
    
    summary_store$data[[key]] <- res
  })
  output$download_all_summary <- downloadHandler(
    filename = function(){
      paste0("All_Summaries_", Sys.Date(), ".docx")
    },
    
    content = function(file){
      req(length(summary_store$data) > 0)
      library(officer)
      library(flextable)
      
      doc <- read_docx()
      for(name in names(summary_store$data)){
        
        res <- summary_store$data[[name]]
        
        #clean_name <- sub("_[0-9]{6}$", "", name)
        # clean_name <- gsub("_", " ", clean_name)
        
        doc <- doc %>%
          body_add_par(paste("Summary:", name), style = "heading 1") %>%
          body_add_flextable(flextable(res) %>% autofit()) %>%
          body_add_par("", style = "Normal")
      }
      
      print(doc, target = file)
      
    }
  )
  
  
  # -----------------------------
  # FOLLOW-UP
  # -----------------------------
  output$followup <- renderTable({
    req(input$run)
    req("Follow-up%" %in% input$analysis_type)
    
    follow_up_calc(
      data(),
      input$time,
      input$status,
      input$event_value,
      input$cutoff
    )
  })
  
  # -----------------------------
  # KM MODEL
  # -----------------------------
  surv_model <- reactive({
    req(input$run)
    req("Kaplan-Meier" %in% input$analysis_type)
    
    df <- data()
    
    df$KM_time <- as.numeric(df[[input$time]])
    df$KM_event <- as.numeric(df[[input$status]]) == input$event_value
    
    df <- df[!is.na(df$KM_time) & !is.na(df$KM_event), ]
    
    if(input$use_group){
      df[[input$group_var]] <- droplevels(as.factor(df[[input$group_var]]))
      fit <-surv_fit( as.formula(paste0("Surv(KM_time, KM_event) ~ ", input$group_var)), data = df )
      
    } else {
      fit <- survfit(Surv(KM_time, KM_event) ~ 1, data = df)
    }
    
    list(fit = fit, df = df)
  })
  
  # -----------------------------
  # KM SUMMARY TABLE
  # -----------------------------
  output$km_summary <- renderTable({
    req(surv_model())
    
    df <- surv_model()$df
    
    if(!input$use_group){
      Survival_summary(df, input$status, input$event_value)
    } else {
      Survival_summary_group(df, input$status, input$event_value, input$group_var)
    }
  })
  
  # -----------------------------
  # KM PROBABILITY
  # -----------------------------
  output$km_prob <- renderTable({
    req(surv_model())
    
    obj <- surv_model()
    s <- summary(obj$fit, times = input$cutoff, extend = TRUE)
    
    if(is.null(s$strata)){
      return(data.frame(
        Group = "Overall",
        Time = input$cutoff,
        Survival_Probability = round(s$surv, 4),
        Std_Error = round(s$std.err, 4)
      ))
    }
    
    groups <- gsub(".*=", "", s$strata)
    
    data.frame(
      Group = groups,
      Time = input$cutoff,
      Survival_Probability = round(s$surv, 4),
      Std_Error = round(s$std.err, 4)
    )
  })
  
  # -----------------------------
  # KM MEDIAN
  # -----------------------------
  output$km_median <- renderTable({
    req(surv_model())
    
    tab <- summary(surv_model()$fit)$table
    
    if(is.null(dim(tab))){
      return(data.frame(
        Group = "Overall",
        Median = ifelse(is.na(tab["median"]), "Not reached", round(tab["median"], 2)),
        Lower_CI = round(tab["0.95LCL"], 2),
        Upper_CI = round(tab["0.95UCL"], 2)
      ))
    }
    
    df <- as.data.frame(tab)
    df$Group <- gsub(".*=", "", rownames(df))
    df$Median <- ifelse(is.na(df$median), "Not reached", round(df$median, 2))
    df$Lower_CI <- round(df$`0.95LCL`, 2)
    df$Upper_CI <- round(df$`0.95UCL`, 2)
    df[, c("Group", "Median","Lower_CI", "Upper_CI")]
  })
  
  #---------------------------------------
  #Download KM summary
  #---------------------------------------
  output$download_km_summary <- downloadHandler(
    filename = function(){
      paste0("KM_Result_", Sys.Date(), ".docx")
    },
    
    content = function(file){
      
      req(surv_model())
      
      library(officer)
      library(flextable)
      obj <- surv_model()
      df <- obj$df
      
      # ---- SUMMARY ----
      if(!input$use_group){
        km_sum <- Survival_summary(df, input$status, input$event_value)
      } else {
        km_sum <- Survival_summary_group(df, input$status, input$event_value, input$group_var)
      }
      
      # ---- PROBABILITY ----
      s <- summary(obj$fit, times = input$cutoff, extend = TRUE)
      
      if(is.null(s$strata)){
        km_prob <- data.frame(
          Group = "Overall",
          Time = input$cutoff,
          Survival_Probability = round(s$surv,3),
          Std_Error = round(s$std.err,3)
        )
      } else {
        km_prob <- data.frame(
          Group = gsub(".*=", "", s$strata),
          Time = input$cutoff,
          Survival_Probability = round(s$surv,3),
          Std_Error = round(s$std.err,3)
        )
      }
      
      # ---- MEDIAN ----
      tab <- summary(obj$fit)$table
      
      if(is.null(dim(tab))){
        km_med <- data.frame(
          Group = "Overall",
          Median = round(tab["median"],2),
          Lower_CI = round(tab["0.95LCL"],2),
          Upper_CI = round(tab["0.95UCL"],2)
        )
      } else {
        df_med <- as.data.frame(tab)
        df_med$Group <- gsub(".*=", "", rownames(df_med))
        
        km_med <- df_med[, c("Group","median","0.95LCL","0.95UCL")]
        names(km_med) <- c("Group","Median","Lower_CI","Upper_CI")
        km_med[, c("Median", "Lower_CI", "Upper_CI")] <-
          lapply(
            km_med[, c("Median", "Lower_CI", "Upper_CI")],
            round,
            digits = 2)
      }
      
      # ---- CREATE DOC ----
      doc <- read_docx()
      
      doc <- doc %>%
        body_add_par("Kaplan-Meier Results", style = "heading 1")
      
      doc <- doc %>%
        body_add_par("Survival Summary", style = "heading 2") %>%
        body_add_flextable(flextable(km_sum)%>% autofit())
      
      doc <- doc %>%
        body_add_par("Survival Probability", style = "heading 2") %>%
        body_add_flextable(flextable(km_prob)%>% autofit())
      
      doc <- doc %>%
        body_add_par("Median Survival", style = "heading 2") %>%
        body_add_flextable(flextable(km_med)%>% autofit())
      
      
      p <- get_km_plot(obj, input)
      
      doc <- doc %>%
        body_add_par("Kaplan-Meier Curve", style = "heading 2") %>%
        body_add_gg(value = p$plot, width = 6, height = 5)
      
      # ---- SAVE ----
      print(doc, target = file)
    }
  )
  
  km_store <- reactiveValues(data = list())
  
  observeEvent(input$run, {
    
    req("Kaplan-Meier" %in% input$analysis_type)
    req(surv_model())
    obj <- surv_model()
    df <- obj$df
    
    # ---- KM SUMMARY ----
    if(!input$use_group){
      km_sum <- Survival_summary(df, input$status, input$event_value)
    } else {
      km_sum <- Survival_summary_group(df, input$status, input$event_value, input$group_var)
    }
    
    # ---- KM PROB ----
    s <- summary(obj$fit, times = input$cutoff, extend = TRUE)
    
    if(is.null(s$strata)){
      km_prob <- data.frame(
        Group = "Overall",
        Time = input$cutoff,
        Survival_Probability = round(s$surv,3),
        Std_Error = round(s$std.err,3)
      )
    } else {
      km_prob <- data.frame(
        Group = gsub(".*=", "", s$strata),
        Time = input$cutoff,
        Survival_Probability = round(s$surv,3),
        Std_Error = round(s$std.err,3)
      )
    }
    
    # ---- KM MEDIAN ----
    tab <- summary(obj$fit)$table
    
    if(is.null(dim(tab))){
      km_med <- data.frame(
        Group = "Overall",
        Median = round(tab["median"],2),
        Lower_CI = round(tab["0.95LCL"],2),
        Upper_CI = round(tab["0.95UCL"],2)
      )
    } else {
      df_med <- as.data.frame(tab)
      df_med$Group <- gsub(".*=", "", rownames(df_med))
      
      km_med <- df_med[, c("Group","median","0.95LCL","0.95UCL")]
      names(km_med) <- c("Group","Median","Lower_CI","Upper_CI")
      km_med[, c("Median", "Lower_CI", "Upper_CI")] <-
        lapply(
          km_med[, c("Median", "Lower_CI", "Upper_CI")],
          round,
          digits = 2)
    }
    # --- KM Plot ---
    p <- get_km_plot(obj, input)
    
    # ---- UNIQUE KEY ----
    key <- if(input$use_group){
      paste0(input$time, "_by_", input$group_var)
    } else {
      paste0(input$time, "_overall")
    }
    
    # ---- STORE ----
    km_store$data[[key]] <- list(
      Summary = km_sum,
      Probability = km_prob,
      Median = km_med,
      Plot = p$plot 
    )
  })
  
  output$download_all_km <- downloadHandler(
    filename = function(){
      paste0("KM_All_Results_", Sys.Date(), ".docx")
    },
    
    content = function(file){
      
      req(length(km_store$data) > 0)
      
      library(officer)
      library(flextable)
      
      doc <- read_docx()
      
      for(name in names(km_store$data)){
        
        block <- km_store$data[[name]]
        # ---- MAIN HEADING ----
        doc <- doc %>%
          body_add_par(paste("Analysis:", name), style = "heading 1")
        # ---- SUMMARY ----
        doc <- doc %>%
          body_add_par("Survival Summary", style = "heading 2") %>%
          body_add_flextable(flextable(block$Summary) %>% autofit())
        
        # ---- PROBABILITY ----
        doc <- doc %>%
          body_add_par("Survival Probability", style = "heading 2") %>%
          body_add_flextable(flextable(block$Probability) %>% autofit())
        
        # ---- MEDIAN ----
        doc <- doc %>%
          body_add_par("Median Survival", style = "heading 2") %>%
          body_add_flextable(flextable(block$Median) %>% autofit())
        #----KMplot----
        doc <- doc %>%
          body_add_par("Kaplan-Meier Curve", style = "heading 2") %>%
          body_add_gg(value = block$Plot, width = 6, height = 5)
        
        # spacing
        doc <- doc %>% body_add_par("", style = "Normal")
      }
      
      print(doc, target = file)
    }
  )
  
  # -----------------------------
  # LOG-RANK
  # -----------------------------
  output$km_pval <- renderText({
    req(input$use_group)
    df <- surv_model()$df
    
    fmla <- as.formula(
      paste0("Surv(KM_time, KM_event) ~ ", input$group_var)
    )
    
    sd <- survdiff(fmla, data = df)
    pval <- sd$pvalue
    
    paste("Log-rank p-value:", signif(pval, 3))
  })
  
  # -----------------------------
  # KM PLOT 
  # -----------------------------
  output$km_plot <- renderPlot({
    req(surv_model())
    
    p <- get_km_plot(surv_model(), input)
    print(p)
  })
  
  ##------------------
  #Cox-Regression
  #-------------------
  output$cox_results <- renderTable({
    req(input$run)
    req("Cox Regression" %in% input$analysis_type)
    req(input$cox_vars)
    
    get_cox_results(data(), input)
  })
  # -----------------------------
  # DOWNLOAD Cox output
  # -----------------------------
  output$download_cox <- downloadHandler(
    filename = function(){
      paste0("Cox_", input$cox_type, "_", Sys.Date(), ".docx")
    },
    
    content = function(file){
      
      req(input$cox_vars)
      
      library(officer)
      library(flextable)
      
      res <- get_cox_results(data(), input)
      
      doc <- read_docx()
      
      doc <- doc %>%
        body_add_par("Cox Regression Results", style = "heading 1") %>%
        body_add_par(input$cox_type, style = "heading 2") %>%
        body_add_flextable(flextable(res) %>% autofit())
      
      print(doc, target = file)
    }
  )
  
  # -----------------------------
  # DOWNLOAD KMplot
  # -----------------------------
  output$download_km <- downloadHandler(
    filename = function(){
      paste0("KM_plot_", Sys.Date(), ".png")
    },
    
    content = function(file){
      
      obj <- surv_model()
      default_labels <- if(!is.null(obj$fit$strata)){
        gsub(".*=", "", names(obj$fit$strata))
      } else NULL
      
      legend_labels <- default_labels
      if(input$legend_labels != "" && !is.null(default_labels)){
        user_labels <- trimws(unlist(strsplit(input$legend_labels, ",")))
        if(length(user_labels) == length(default_labels)){
          legend_labels <- user_labels
        }
      }
      x_max <- ifelse(is.na(input$xmax), input$cutoff, input$xmax)
      
      legend_title <- if(input$legend_title == ""){
        if(input$use_group) input$group_var else "Overall Survival"
      } else input$legend_title
      
      x_label <- if(input$xlab == "") "Time" else input$xlab
      y_label <- if(input$ylab == "") "Survival Probability" else input$ylab
      
      p <- ggsurvplot(
        obj$fit,
        data = obj$df,
        conf.int = TRUE,
        pval = if(input$use_group) TRUE else FALSE,
        xlim = c(0, x_max),
        break.time.by = input$xbreak,
        xlab = x_label,
        ylab = y_label,
        ggtheme = theme_classic(),
        legend.title = legend_title,
        legend.labs = legend_labels
      )
      
      ggsave(file, plot = p$plot, width = 7, height = 5,dpi = 300,bg="white")
    }
  )
}

shinyApp(ui = ui, server = server)