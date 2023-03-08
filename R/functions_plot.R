# Header #############################################################
# 
# Author: Lisa Nicvert
# Email:  lisa.nicvert@univ-lyon1.fr
# 
# Date: 2022-12-15
#
# Script Description: plot various data


# Libraries ---------------------------------------------------------------
library(ggplot2)
library(lubridate)
library(ggtext)
library(tidygraph)
library(ggraph)
library(patchwork) # to add ggplot objects together
library(pammtools) # for stepribbon

# Model -------------------------------------------------------------------
plot_interactions <- function(ue_df, 
                              scale = "days",
                              title = NA,
                              relative = FALSE,
                              silhouettes = NA,
                              timestep = NA,
                              ystep = NA,
                              textsize = 10,
                              linesize = .5, 
                              separate_self = FALSE
){
  # Plot interaction functions
  # ### Inputs
  # ue_df: dataframe with results of UnitEvents inference. Must have columns:
  #   time
  #   excitefunc
  #   from
  #   to
  # scale: days or hours, following whether we want the time axis
  #   graduated in days or hours
  # title: plot title
  # relative: plot intensity absolute or relative (divided by spont) value?
  # silhouettes: optional labels with animal silhouettes to replace 
  #   default labels.
  #   if it exists, must be a named vector with each name 
  #   corresponding to a species name in ue_df. The elements are 
  #   markdown codes containing a <img src='pato_to_image'/> element.
  # timestep: optional timestep for x-axis.
  # ystep: optional step for y-axis (function values).
  # textsize: text minimal size (for x and y axes)
  # linesize: linewidth
  # separate_self: whether to separate auto-interactions and plot them above
  # ### Output
  # ggplot object, a plot with the pairwise interaction functions 
  #   between species.
  
  ue_df_plot <- ue_df
  
  if(scale == "hours"){
    ue_df_plot$time <- ue_df_plot$time*24
  }
  
  if(relative){
    ue_df_plot$excitefunc <- ue_df_plot$excitefunc/ue_df_plot$spont
  }
  
  # If timestep not set
  if(is.na(timestep)){
    timestep <- ue_df_plot$time[2] - ue_df_plot$time[1]
  }
  
  if (separate_self) {
    self <- ue_df_plot %>% filter(from == to) 
    ue_df_plot$excitefunc[ue_df_plot$from == ue_df_plot$to] <- 0
    
    # Add dummy silhouette for self (will be set to element_blank later)
    silhouettes_self <- c(silhouettes,
                          "*Self-interactions*")
    names(silhouettes_self)[length(silhouettes_self)] <- "Self-interactions"
  }
  
  g <- ggplot(ue_df_plot) + 
    geom_step(aes(x=time, y=excitefunc), linewidth = linesize) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = linesize) + 
    scale_x_continuous(paste0("Time (", scale, ")"),
                       sec.axis = dup_axis(name = "To..."),
                       breaks = seq(0, max(ue_df_plot$time), by = timestep),
                       limits = c(0, max(ue_df_plot$time))) +
    ylab("From...") +
    theme_linedraw() +
    theme(axis.text.x.top = element_blank(),
          axis.ticks.x.top = element_blank(),
          axis.text.x = element_text(size = textsize),
          axis.text.y = element_text(size = textsize),
          axis.title.x.bottom = element_text(size = textsize*1.16, 
                                             face = "italic"),
          axis.title = element_text(size = textsize*1.5),
          title = element_text(size = textsize*1.5),
          panel.spacing = unit(0.5, "lines"))
  
  if (separate_self) {
    mar <- 0.2
    g <- g +
      geom_rect(data = subset(ue_df_plot, from == to),
                fill = "white", xmin = -Inf, xmax = Inf,
                ymin = -Inf, ymax = Inf) +
      ggtitle("Cross-species interactions")
    
    g2 <- ggplot(self) + 
      geom_step(aes(x=time, y=excitefunc), linewidth = linesize) +
      geom_hline(yintercept = 0, linetype = "dashed", linewidth = linesize) + 
      theme_linedraw() +
      theme(axis.text.x.top = element_blank(),
            axis.ticks.x.top = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = textsize),
            axis.text.y = element_text(size = textsize),
            axis.title = element_text(size = textsize*1.5),
            title = element_text(size = textsize*1.5),
            panel.spacing = unit(0.5, "lines")) +
      ggtitle("Auto-interactions")
  }
  
  
  if(!is.na(ystep)){
    mn <- floor(min(ue_df_plot$excitefunc)/ystep)*ystep
    mx <- ceiling(max(ue_df_plot$excitefunc)/ystep)*ystep
    
    g <- g + scale_y_continuous(breaks = seq(mn, 
                                             mx, 
                                             by = ystep))
  }
  
  # if there are silhouettes
  if(!all(is.na(silhouettes))){
    g <- g + 
      facet_grid(rows = vars(from),
                 cols = vars(to),
                 labeller = as_labeller(silhouettes),
                 switch = "y") +
      theme(strip.text.x = element_markdown(color = "black", 
                                            size = textsize),
            strip.text.y.left = element_markdown(color = "black", 
                                                 size = textsize,
                                                 angle = 0),
            strip.background = element_rect(fill="white"))
    
    if (separate_self) {
      g2 <- g2 + 
        facet_grid(cols = vars(to),
                   labeller = as_labeller(silhouettes_self)) +
        theme(strip.text.x = element_markdown(color = "black", 
                                              size = textsize),
              strip.text.y.left = element_markdown(color = "black",
                                                   size = textsize,
                                                   angle = 0),
              strip.placement = "outside",
              strip.background = element_rect(fill="white"))
    }
    
  }else{
    g <- g +
      facet_grid(rows = vars(from),
                 cols = vars(to),
                 switch = "y") +
      theme(strip.text = element_text(size = textsize,
                                      color = "black"),
            strip.background = element_rect(fill="white"))
    
    if (separate_self) {
      g2 <- g2 +
        facet_grid(cols = vars(to)) +
        theme(strip.text = element_text(size = textsize,
                                        color = "black"),
              strip.background = element_rect(fill="white"),
              strip.placement = "outside")
    }
  }
  
  # Add title in case single plot
  if (separate_self) {
    if(!is.na(title)){
      g <- g + ggtitle(title) +
        theme(strip.placement = "outside")
    }
  }
  
  nspecies <- length(unique(ue_df_plot$from))
  
  if (separate_self) {
    glist <- g2 + plot_spacer() + g + 
      plot_layout(heights = c(1, 0.3, nspecies),
                  ncol = 1)
    
    # Format axes
    glist[[1]] <- glist[[1]] + 
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x.bottom = element_blank())
    
    if (!is.na(title)) { 
      glist <- glist + title(title)
    }
    glist
  } else {
    g
  }
}

plot_background_rate <- function(ue_df,
                                 title = NA,
                                 textsize = 10,
                                 silhouettes = NA,
                                 write_label = FALSE,
                                 nudge_label = 0.3) {
  # Plot the background rate for species.
  # ### Inputs
  # ue_df: dataframe with results of UnitEvents inference. Must have columns:
  #   time
  #   excitefunc
  #   from
  #   to
  # title: plot title
  # textsize: text minimal size (for x and y axes)
  # silhouettes: optional labels with animal silhouettes to replace 
  #   default labels.
  #   if it exists, must be a named vector with each name 
  #   corresponding to a species name in ue_df. The elements are 
  #   markdown codes containing a <img src='pato_to_image'/> element.
  # write_label: write background rates values besides the points?
  # nudge_label: if the background rates are written, by how much should they 
  #   be nudged on the x-axis?
  # ### Output
  # a ggplot object, representing background rates for each species.
  
  spont_plot <- ue_df %>% group_by(to) %>%
    summarise(spont = unique(spont)) %>%
    rename("species" = "to")
  
  g <- ggplot(spont_plot, aes(x = species, y = spont)) + 
    geom_point() +
    theme_linedraw() +
    theme(axis.title.x = element_blank(),
          axis.text = element_text(size = textsize),
          axis.title.y = element_text(size = textsize*1.2),
          title = element_text(size = textsize*1.2)) +
    ylab("Background rate (day⁻¹)")
  
  if(write_label) {
    g <- g + 
      geom_text(aes(label = round(spont, 3)), nudge_x = nudge_label)
  }
  
  if(!is.na(title)) {
    g <- g + ggtitle(title)
  }
  
  if(!all(is.na(silhouettes))){
    g <- g + scale_x_discrete(labels = labs) +
      theme(axis.text.x = element_markdown(color = "black",
                                           size = textsize))
  }
  
  g
}

plot_graph <- function(g, layout = c(), 
                       repel = FALSE,
                       coledges = 'grey', 
                       colnodes = 'cadetblue3', 
                       coltext = "black",
                       textsize = 5, s=8,
                       arrsize = 6,
                       nudge_x = NA,
                       parse_labels = FALSE,
                       use_labels_column = FALSE){
  # Plots a graph g.
  # ### Inputs
  # layout: optional layout (defaults to layout_in_circle)
  # repel: repel node labels?
  # coledges: color of edges (a string color name)
  # colnodes: color of nodes (named vector named as species or a string
  #   color name). If it is NULL, species will be colored automatically.
  # coltext: color of text for the node labels. It can be a vector
  #   (will correspond to species alphabetical order),
  #   a named vector named as species or a string color name. 
  # textsize: node label size
  # s: size of nodes
  # arrsize: size of arrow
  # nudge_x: optional vector to move x labels for node labels.
  # parse_labels: parse text labels?
  # use_labels_column: use column 'label' for nodes labels?
  # ### Output
  # A ggplot object representing the graph of the Hawkes interaction model.
  
  if(length(layout) == 0){
    layout <- layout_in_circle(g)
  }
  if(is.na(nudge_x)){
    nudge_x <- rep(0, length(V(g)))
  }
  ar <- arrow(angle=30,length = unit(arrsize,"mm"),
              ends="last", type = "closed")
  alpha <- 1
  
  gr <- ggraph(g, layout = layout) + 
    geom_edge_fan(aes(width = weight),
                  arrow = ar,
                  end_cap = circle(s/2.1,unit="mm"),
                  alpha = alpha, 
                  colour = coledges) +
    geom_edge_loop(aes(width = weight,
                       span = 90,
                       direction = 45), 
                   arrow = ar,
                   end_cap = circle(s/2.2,unit="mm"),
                   alpha = alpha, colour = coledges)
  
  if(!all(is.null(colnodes))){ # not null
    if (length(colnodes) == 1) { # unique colour
      gr <- gr +
        geom_node_point(size = s, colour = colnodes)
    } else {
      gr <- gr +
        geom_node_point(size = s, aes(colour = names),
                        show.legend = FALSE) +
        scale_color_manual(values = colnodes)
    }
  } else{
    gr <- gr +
      geom_node_point(size = s, aes(colour = names), 
                      show.legend = FALSE) 
  }
  
  # Reorder coltext
  if(length(coltext) != 0 & !is.null(names(coltext))) {
    coltext <- coltext[sort(names(coltext))]
  }
  
  if(use_labels_column) {
    gr <- gr +
      geom_node_text(aes(label = labels), 
                     size = textsize,
                     nudge_x = nudge_x,
                     color = coltext,
                     repel = repel,
                     parse = parse_labels)
  } else {
    gr <- gr +
      geom_node_text(aes(label = names), 
                     size = textsize,
                     nudge_x = nudge_x,
                     color = coltext,
                     repel = repel,
                     parse = parse_labels,
                     fontface = "italic")
  }
  
  mar <- 0.2
  
  gr <- gr +
    scale_edge_width(range = c(0.2, s/2), guide = "none") +
    ylim(c(min(layout[,2])-mar, max(layout[,2]) + mar)) +
    xlim(c(min(layout[,1])-mar, max(layout[,1]) + mar)) +
    theme_void()
  
  gr
}

plot_observed_rate <- function(rates, 
                               data, 
                               timestep = 2, 
                               textsize = 20,
                               ptsize = 1,
                               lwd = 1,
                               t1, t2,
                               hlambda, hpoints,
                               minor_spacing,
                               major_spacing,
                               max_lambda,
                               ybreaks,
                               cols,
                               ylabel = TRUE,
                               xlab = "Time (days)"){
  # Plot the rate of an observed sequence of events.
  # ### Inputs
  # rates is a dataframe with pre-computed rates. Must have columns:
  #   time
  #   lambda
  #   species
  # data is the corresponding occurrence data.
  #   columns:  stamp and species.
  # timestep: timestep to plot
  # textsize: base text size (axes text, 
  #   axes labels are a 2 units bigger)
  # ptsize: point sizes for the events
  # t1 and t2 are time bounds for subsetting data. If missing, then all data will be plotted.
  # lwd: linewidth for intensities
  # hlambda and hpoints are plot relative heights resp. for lambda and points
  # minor_spacing and major_spacing: the space to add between plots. 
  #   minor_spacing is the space between the intensity function 
  #   and the correponding points.
  #   major_spacing is the space between points and the following intensity.
  # Default resp. to hpoints/2 and hpoints.
  # max_lambda: max value for y-axis for rates
  # ybreaks: breaks for y-axis (for lambda plots)
  # cols: named vector for colors: names are species names and
  #   contain colors 
  # ylabel: display labels ?
  # xlab: xlabel to display (optional)
  # ### Output 
  # A ggplot object generated with patchwork.
  #   Multiple plots in the same column, where all plots are paired,
  #   the top plot representing the intensity and the bottom plot the actual
  #   occurrences for one species.
  
  if(missing(t1)){
    t1 <- floor(min(data$stamp))
  }
  if(missing(t2)){
    t2 <- max(data$stamp)
  }
  
  # Filter rates and data 
  res_sub <- rates %>% filter(time <= t2 & time >= t1)
  obs_sub <- data %>% filter(stamp <= t2 & stamp >= t1)
  
  # Add min and max times
  if(min(res_sub$time) > t1){
    # Create baseline table
    tab <- res_sub %>% group_by(species) %>% 
      filter(time == min(time))
    tab$time <- t1
    res_sub <- res_sub %>% bind_rows(tab) %>%
      arrange(species, time)
  }
  
  if(max(res_sub$time) < t2){
    # Create baseline table
    tab <- res_sub %>% group_by(species) %>% 
      filter(time == max(time))
    tab$time <- t2
    res_sub <- res_sub %>% bind_rows(tab) %>%
      arrange(species, time)
  }
  
  # Rename data column stamp -> time
  obs_sub <- obs_sub %>% rename("time" = "stamp")
  
  # Get species list (only species seen during time interval)
  species_list <- unique(obs_sub$species)
  
  # Colour list
  if(missing(cols)){
    cols <- rep("black", length(species_list))
    names(cols) <- species_list
  }
  
  # Get maximum value for rate (common axis)
  if(missing(max_lambda)) {
    max_lambda <- max(res_sub$lambda)
  }
  
  
  # Get maximum time bounds
  tmin <- t1
  tmax <- t2
  
  # Initialise graphs list
  glist <- list()
  
  for(i in 1:length(species_list)){ # for each species
    sp <- species_list[i]
    
    # Select species data
    sub_spp <- res_sub %>% filter(species == sp)
    obs_spp <- obs_sub %>% filter(species == sp)
    
    l <- ggplot(sub_spp) + 
      geom_step(aes(x = time, y = lambda), col = cols[sp],
                linewidth = lwd) +
      theme_linedraw() +
      theme(text = element_text(size = textsize), 
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_text(angle = 0, vjust = 0.5,
                                        size = textsize*1.1))
    
    if(ylabel){
      l <- l + 
        ylab(as.expression(substitute("\u03BB"[i], list(i = as.name(sp)))))
    }
    else{
      l <- l + theme(axis.title.y = element_blank())
    }
    
    if(!missing(ybreaks)){
      l <- l + 
        scale_y_continuous(breaks = ybreaks,
                           limits = c(0, max_lambda),
                           expand = c(0.01, 0))
    }else{
      l <- l + ylim(c(0, max_lambda))
    }
    pt <- ggplot(obs_spp) + 
      geom_point(aes(x=time, y = 1), col = cols[sp], size = ptsize) +
      theme_linedraw() +
      ylab(sp)
    
    if(i != length(species_list)){ # it's not the last species
      pt <- pt + theme(text = element_text(size = textsize), 
                       axis.title.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid.major.y = element_blank(),
                       panel.grid.minor.y = element_blank())
    }else{
      pt <- pt + theme(text = element_text(size = textsize), 
                       axis.title.x = element_text(size = textsize*1.1),
                       axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),
                       panel.grid.major.y = element_blank(),
                       panel.grid.minor.y = element_blank()) +
        xlab(xlab)
    }
    
    if(length(glist) == 0){
      glist <- plot_spacer() + l + plot_spacer() + pt + plot_spacer()
    }else{
      glist <- glist + l + plot_spacer() + pt + plot_spacer()
    }
  }
  
  if (missing(minor_spacing)) {
    minor_spacing <- hpoints/2
  }
  
  if (missing(major_spacing)) {
    major_spacing <- hpoints
  }
  
  # Once all species-specific plots have been done, assemble them
  (glist & theme(plot.margin = margin(t = 0, b = 0, l = 3, r = 5)) 
    & scale_x_continuous(limits = c(tmin, tmax),
                         breaks = seq(tmin, tmax, by = timestep),
                         expand = c(0.001, 0))) +
    plot_layout(ncol = 1, 
                heights = c(0.5*hpoints, # Top spacer
                            rep(c(hlambda, minor_spacing, hpoints, major_spacing), length(species_list) -1),
                            c(hlambda, minor_spacing, hpoints, 0.5*hpoints)) # Last plot with bottom spacer
                )
  
}


# Simulation --------------------------------------------------------------
plot_interactions_simu <- function(df, 
                       title = NA, 
                       level = 0.05){
  # Plots the inferred interaction functions vs the real function
  # ### Inputs
  # df: a dataframe containing values for the inferred function(s)
  #   and the true function. It has columns:
  #   rep: repetition ID (s... for simul, true for true data)
  #   excitefunc: value of the excitation function
  #   from/to: ID of the species
  #   time: time in days
  # title: optional plot title
  # level: confidence level to use around simulations
  # ### Output
  # A ggplot of the interaction functions where the true function 
  #   is in red and the inferred function in blue with a confidence 
  #   interval (if there were several inferred models.)
  
  simul <- df %>% filter(rep != "true")
  true <- df %>% filter(rep == "true")
  
  simul <- simul %>% group_by(from, to, time)
  
  simul <- simul %>%
    summarise(inf = quantile(excitefunc, level/2)[[1]],
              sup = quantile(excitefunc, 1 - (level/2))[[1]],
              median = median(excitefunc),
              .groups = "drop")
  
  g <- ggplot(data = simul, aes(x = time))
  
  g <- g + geom_stepribbon(aes(ymin = inf, ymax = sup), 
                           fill = "cornflowerblue", alpha = 0.5) +
    geom_step(aes(y = median), col = "blue") +
    geom_step(data = true, 
              aes(y = excitefunc), col = "red")
  
  g <-  g +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed") +
    facet_grid(from ~ to, switch="y") +
    theme_linedraw() +
    theme(panel.grid.minor = element_blank()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_rect(fill = "white", colour = "black", 
                                          linewidth = 1),
          strip.text = element_text(color = "black")) +
    ylab("Excitation functions") + xlab("Time (days)")
  
  if(!is.na(title)){
    g <- g + ggtitle(title)
  }
  g
}

plot_background_rate_simu <- function(df, alpha = 0.05, title = NA){
  # Plots the inferred background rates vs the real rates
  # ### Inputs
  # df: the dataframe with true and estimated intensities. It has columns:
  #   rep: repetition ID (s... for simul, true for true data)
  #   excitefunc: value of the excitation function
  #   from/to: ID of the species
  #   time: time in days
  # alpha: confidence level for plotting 
  # title: optional plot title
  # ### Output
  # A ggplot of the background rates where the true rate 
  #   is in red and the inferred rate in blue with a confidence 
  #   interval (if there were several inferred models.)
  
  # Get estimates
  spont_est <- df %>% filter(rep != "true") 
  # Get true
  spont_true <- df %>% filter(rep == "true")
  
  # Get only unique values for true
  spont_true <- spont_true %>%
    select(to, spont) %>%
    unique()
  
  # Summarise intercepts estimates
  spont_est <- spont_est %>%
    select(to, spont, rep) %>%
    unique() %>% # get only one coeff per repetition (else repeated for each value of the filter function)
    group_by(to) %>%
    summarise(spont_median = median(spont),
              spont_upper = quantile(spont, 
                                     1-alpha/2),
              spont_lower = quantile(spont,
                                     alpha/2))
  
  g <- ggplot() +
    geom_errorbar(col = "cornflowerblue",
                  aes(x = to,
                      ymin = spont_lower, 
                      ymax = spont_upper), 
                  data = spont_est,
                  show.legend = FALSE) +
    geom_point(col = "blue", 
               aes(x = to, y = spont_median),
               data = spont_est) +
    geom_point(col = "red", 
               aes(x = to, y = spont),
               data = spont_true) +
    theme_linedraw()
  
  if(!is.na(title)){
    g <- g + ggtitle(title)
  }
  g
  
}

plot_perf <- function(d, 
                      xaxis = "trapping_days", 
                      thr, 
                      vline,
                      psize = 1){
  # Plots the performance of the inference.
  # d: the dataframe of observed sensitivity and specificity.
  #   It has columns:
  #   value (sensitivity or specificity value)
  #   type (sensitivity or specificity as "sensi" or "speci")
  #   "xaxis" (a measure of trapping days)
  #     possibly qinf, qsup (then quantiles are plotted arould sensitivity
  #     and specificity values.)
  # xaxis: name of the x-axis to choose for plotting 
  #   (must be present in d)
  # thr: optional threshold where to plot a horizontal line
  # vline: optional vline to plot to draw attention to a specific time
  # psize: point sizes
  # ### Output
  # a ggplot object with the sensitivity and specificity displayed along
  #   xaxis.
  
  cols <- c("#1B9E77", "#7570B3")
  names(cols) <- c("sensi", "speci")
  
  gbase <- ggplot(d, aes(x = !!sym(xaxis), 
                         col = type, group = type)) +
    geom_point(aes(y = value, shape = type), size = psize) +
    geom_line(aes(y = median, linetype = type)) +
    theme_linedraw() +
    theme(axis.title.y.left = element_blank(),
          strip.background = element_rect(fill = "white", colour = "black", 
                                          linewidth = 1),
          strip.text = element_text(color = "black"),
          legend.position = "bottom") +
    scale_color_manual(name = "Performance",
                       labels = c("true positive rate",
                                  "true negative rate"),
                       values = cols) +
    scale_linetype(name = "Performance",
                   labels = c("true positive rate",
                              "true negative rate")) +
    scale_shape_manual(name = "Performance",
                       labels = c("true positive rate",
                                  "true negative rate"),
                       values = c(19, 1))
  
  # Add quartiles
  if("qinf" %in% colnames(d) & "qsup" %in% colnames(d)){
    gbase <- gbase +
      geom_ribbon(aes(ymin=qinf, ymax=qsup, fill = type),
                  alpha = 0.2,
                  color = NA,
                  show.legend = FALSE)  +
      scale_fill_manual(values = cols)
  }
  
  if(!missing(thr)){
    gbase <- gbase + 
      geom_hline(yintercept = thr, linetype = "dashed", color = "darkgrey")
  }
  if(!missing(vline)){
    gbase <- gbase + 
      geom_vline(xintercept = vline, linetype = "dashed")
  }
  gbase
}


# Bias --------------------------------------------------------------------

plot_bias <- function(bias_df, fill = "valprop", textsize = 12) {
  # Plots the bias for each interaction from several inferences.
  # ### Inputs
  # bias_df: a dataframe containing the summarized results of the 
  #   inference for several repetitions. Must have columns:
  #   from and to (interacting species)
  #   the column specified with 'fill'
  # fill: name of the column to use for values to fill
  #   the plot.
  # textsize: size of the text
  ggplot(bias_df) +
    geom_tile(aes(y = from, x = to, fill = !!sym(fill))) +
    scale_fill_viridis(limits = c(0, 1), name = "Proportion") +
    coord_equal(expand = FALSE) +
    theme_linedraw() +
    theme(text = element_text(size = textsize),
          axis.title = element_text(size = textsize*1.05),
          legend.title = element_text(size = textsize*1.05),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill="white"))
}
