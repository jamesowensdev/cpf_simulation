# library(tidyverse)
# library(terra)
# library(R6)
# library(gganimate)
# library(gifski)
# library(ggforce)
# library(checkmate)

load_dependencies <- function() {
  required_packages <- c("tidyverse", "terra", "R6", "gganimate", "gifski", "ggforce", "checkmate")

  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

  if(length(missing_packages) > 0) {
    stop("The following required packages are missing: ", 
         paste(missing_packages, collapse = ", "), 
         ".\nPlease install them manually using install.packages() to run this simulation.")
  }
  
  message("Loading dependencies...")
  invisible(lapply(required_packages, library, character.only = TRUE))
  message("Dependencies loaded.")
}

load_dependencies()

ForagerAgent <- R6Class("ForagerAgent",
  public = list(
    x = NULL,
    y = NULL,
    nest_x = NULL,
    nest_y = NULL,
    inertia = NULL,
    velocity_x = 0,
    velocity_y = 0,
    state =  "Foraging",
    return_speed = 0.8,
  
    initialize = function(nest_x, nest_y, inertia) {
      self$nest_x <- nest_x
      self$nest_y <- nest_y
      self$x <- nest_x
      self$y <- nest_y
      self$inertia <- inertia

    },
    
    step = function(env_raster) {
      coords <- cbind(self$x, self$y)
      if(self$state == "Foraging") {
        self$move()
        self$check_forage(env_raster)
      } else
        self$move()
    },
    
    move = function() {
      if(self$state == "Returning") {
        self$move_return()
      } 
      else 
        self$move_randomly_inertia()
      
    },

    check_forage = function(env_raster) {
    # Extracts cell_id from the environment raster, checks cell_id for a 0 or 1 value
    # If a value of 1 is found at the cell_id, ofrage has been found, the forage is
    # removed from the cell and ther agent is changes state.
      cell_id <- cellFromXY(env_raster, cbind(self$x, self$y))
      
      if (!is.na(cell_id)) {
        food_val <- env_raster[cell_id][1,1] 
        
        if (!is.na(food_val) && food_val == 1) {
          self$state <- "Returning"
          env_raster[cell_id] <- 0 
        }
      }
    },

    move_randomly_inertia = function() {
      inertia <- self$inertia 
      magnitude <- 0

      brownian_x <- rnorm(1, mean = 0, sd = 1)
      brownian_y <- rnorm(1, mean = 0, sd = 1)

      # If the agent is not moving we contribute the movement vector solely to brownian motion
      if(self$velocity_x == 0 && self$velocity_y == 0) {
        new_dx <- brownian_x
        new_dy <- brownian_y
      } else {
        new_dx <- (self$velocity_x * inertia) + (brownian_x * (1 - inertia))
        new_dy <- (self$velocity_y * inertia) + (brownian_y * (1 - inertia))
      }
      
      magnitude <- sqrt(new_dx^2 + new_dy^2)
      
      # panic clause to prevent rate case where  
      if(magnitude < 0.0001) {
        random_angle <- runif(1, 0, 2 * pi)
      
        dx <- cos(random_angle)
        dy <- sin(random_angle)
      } else {
        dx <- (new_dx / magnitude)
        dy <- (new_dy / magnitude)
      }

      self$x <- self$x + dx
      self$y <- self$y + dy

      self$velocity_x <- dx
      self$velocity_y <- dy

      self$boundary_check()
  
    },

    move_return = function() {
      # Calculates the dx/dy from the agents current coords to the nest site.
      diff_x <- self$nest_x - self$x
      diff_y <- self$nest_y - self$y
      
      dist <- sqrt(diff_x^2 + diff_y^2)

      # Resets agents state when it comes within a step of the nest site.
      if (dist < 0.8) {
        self$x <- self$nest_x
        self$y <- self$nest_y

        self$state <- "Foraging" 

      } else {
        speed <- self$return_speed
        self$x <- self$x + (diff_x / dist) * speed
        self$y <- self$y + (diff_y / dist) * speed
      }
    },

    boundary_check = function() {
      # If x or y pos is out of bounds (OOB), sets OOB parameter to 0 or 100 dependent on OOB border 
      # Direction of travel is reversed with random noise addition.

      hit_wall <- FALSE
      
      if (self$x < 0) {
        self$x <- 0
        self$velocity_x <- abs(self$velocity_x) 
        hit_wall <- TRUE
        
      } else if (self$x > 100) {
        self$x <- 100
        self$velocity_x <- -abs(self$velocity_x) 
        hit_wall <- TRUE
      }

      if (self$y < 0) {
        self$y <- 0
        self$velocity_y <- abs(self$velocity_y) 
        hit_wall <- TRUE
      } else if (self$y > 100) {
        self$y <- 100
        self$velocity_y <- -abs(self$velocity_y)
        hit_wall <- TRUE
      }

      if (hit_wall) {
        self$velocity_x <- self$velocity_x + rnorm(1, 0, 0.5)
        self$velocity_y <- self$velocity_y + rnorm(1, 0, 0.5)
      }
    }
  )
)

generate_animation <- function(results, camera_coords, detections_df, radius, sim_name) {
  message("Generating plot for animation")
  plot_data <-
    results |>
    as_tibble() |>
    mutate(step = as.integer(step))

  if (!is.null(detections_df) && nrow(detections_df) > 0) {
    plot_data <- 
      plot_data |>
      mutate(
        is_detected = as.integer(step %in% detections_df$step),
        total_detections = cumsum(is_detected)
      )
  } else
    plot_data$total_detections <- 0

  animation_plot <- ggplot(plot_data, aes(x = x, y = y)) +
                      geom_path(aes(color = state, linewidth = state, group = 1), alpha = 0.7, ) +
                      geom_point(color = "red", size = 4) +
                      geom_circle(data = camera_coords, aes(x0 = x, y0 = y, r = radius, linetype = "Camera Radius"), 
                        inherit.aes = FALSE, color = "blue", alpha = 0.1, show.legend = FALSE ) +
                        geom_point(data = camera_coords, aes(x = x, y = y, shape = "Camera Trap"), 
                         size = 1, color = "blue", stroke = 1.5, inherit.aes = FALSE) +
                      annotate("point", x = 50, y = 50, size = 3, stroke = 2, color = "green") +
                      geom_text(aes(x = 2, y = 98, label = paste("Detections:", total_detections)), 
                        hjust = 0, 
                        size = 4, 
                        fontface = "bold") +
                      scale_color_manual(values = c("Foraging" = "grey", 
                               "Returning" = "red")) +
                      scale_linewidth_manual(values = c("Foraging" = 0.7, "Returning" = 0.5)) +
                      scale_linetype_manual(values = c("Camera Radius" = "dotted")) +
                      scale_shape_manual(values = c("Camera Trap" = 4)) +
                      coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
                      theme_minimal() +
                      transition_reveal(step) +
                      labs(
                        title = sim_name,                       
                        subtitle = "Simulation Step: {frame_along}", 
                        x = "X Coordinate", 
                        y = "Y Coordinate",
                        shape = "Camera"
                      ) +
                      theme(plot.title = element_text(size = 10, face = "bold"))
  
  folder_name <- "simulation_gifs"
  if (!dir.exists(folder_name)) {
    dir.create(folder_name)
  }
  filename <- paste(folder_name,"/",sim_name,".gif", sep = "")

  max_frames <- nrow(results)

  message("Animating plot")
  
  suppressMessages({
    animate(animation_plot, max_frames, fps = 30, renderer = gifski_renderer())
  })

  anim_save(filename, animation = last_animation())
  message(paste(filename, " was saved successfully to", getwd()))
}

set_environment <- function(seed, p_food){
  # Setting up the world, declaring a seed allows for reproducible random movements between runs
  env <- rast(nrows = 100, ncols = 100, xmin = 0, xmax = 100, ymin = 0, ymax = 100)
  set.seed(seed)
  values(env) <- sample(c(0, 1), size = ncell(env), replace = TRUE, prob = c(1 - p_food, p_food))
  return(env)
}

CameraTrap <- R6Class("CameraTrap",
          public = list(
            x = NULL,
            y = NULL,
            radius = NULL,
            p_detection_max = NULL,
            inf_prop = NULL,
            steepness = NULL,
            last_triggered = -9999,
            cooldown = NULL,

            initialize = function(x, y, radius, p_detection_max, inf_prop, steepness, cooldown) {
              self$x <- x
              self$y <- y
              self$radius <- radius
              self$p_detection_max <- p_detection_max
              self$inf_prop <- inf_prop
              self$steepness <- steepness
              self$cooldown <- cooldown
    },

    agent_detection = function(camera_id, agent_x,agent_y, current_step) {
      dx = self$x - agent_x
      dy = self$y - agent_y
      dist_sq = dx^2 + dy^2
      inf_dist <- self$radius * self$inf_prop
      timer <- current_step - self$last_triggered

      # If camera is on cooldown we skip all probability calcs
      if(timer < self$cooldown) return(NULL)

      # sqrt() is computationally expensive, we don't want to be doing it unless required.
      # If agent is outside radius no need to calcaute p_detection, skip.
      if(dist_sq > self$radius^2) return(NULL)
      
      dist <- sqrt(dist_sq)

      # Creates a reverse sigmoid function of detection probability against distance from camera
      # inf_dist equates to the distance from the camera where p_sigmoid = 0.5
      p_sigmoid <- 1 / (1 + exp(self$steepness * (dist - inf_dist)))
    
      p_detection <- self$p_detection_max * p_sigmoid

      # Random number between 0-1 generated to determine if detected or not.
      if(runif(1) < p_detection){
        self$last_triggered <- current_step
        return(list(camera_id = camera_id ,camera_x = self$x, camera_y = self$y))
      } else {
        return(NULL)
      }
    }
  )
)

  setup_camera_traps <- function(coords_df, radius, p_max, inf_prop, steepness, cooldown) {

    trap_list <- list()
    
    for(i in 1:nrow(coords_df)) {
      trap_list[[i]] <- CameraTrap$new(
        x = coords_df$x[i],
        y = coords_df$y[i],
        radius = radius,
        p_detection_max = p_max,
        inf_prop = inf_prop,
        steepness = steepness,
        cooldown = cooldown
      )
    }
    return(trap_list)
  }

  check_all_traps <- function(agent, trap_list, current_step) {
    step_detections <- list()
    detection_count <- 0

    for(i in seq_along(trap_list)) {
      trap <- trap_list[[i]]
      results <- trap$agent_detection(i, agent$x, agent$y, current_step)

      if(!is.null(results)) {
        detection_count <- detection_count + 1

        step_detections[[detection_count]] <- list(
          camera_id = i,
          step = current_step,
          camera_x = results$camera_x,
          camera_y = results$camera_y
        )
      }
    }
      if(detection_count > 0) return(step_detections) else return(NULL)
  
  }

run_simulation <- function(
  sim_name = "cpf_simulation",
  steps = 1000,
  seed = 619,
  nest_coords = c(50,50),
  inertia = 0.7,
  p_food = 0.02,
  camera_coords = NULL,
  radius = NULL,
  p_detection_max = 0.9,
  inf_prop = 0.5,
  steepness = 0.8,
  cooldown = NULL,
  create_gif = FALSE){
  
  start_time <- Sys.time()
  message("Setting up R Enviroment...")
  load_dependencies()
  message("")

  coll <- makeAssertCollection()
  assert_data_frame(camera_coords, min.rows = 1, add = coll)
  assert_number(radius, lower = 0, finite = TRUE, add = coll) 
  assert_number(cooldown, lower = 0, add = coll)
  reportAssertions(coll)


  

  message("**************************************************")
  message(paste("--- Starting Simulation:", sim_name, "---"))
  message("**************************************************")

  ##############################################################
  #                      WORLD SETUP                           #
  ##############################################################
  message(paste("Setting up environment..."))
  env <- set_environment(seed, p_food)
  
  ##############################################################
  #                     FORAGER SETUP                          #
  ##############################################################
  message(paste("Initalising Foraging Agent..."))
  forager <- ForagerAgent$new(nest_x = nest_coords[1], nest_y = nest_coords[2], 
                              inertia = inertia)
  
  # Preallocating results vectors to be combined into tibble after simulation
  pos_x <- numeric(steps + 1)
  pos_y <- numeric(steps + 1)
  agent_state <- character(steps + 1)

  pos_x[1] <- forager$x
  pos_y[1] <- forager$y
  agent_state[1] <- "Foraging"

  ##############################################################
  #                   CAMERA TRAP SETUP                        #
  ##############################################################
  message("Initalising Camera Traps...")
  traps <- setup_camera_traps(camera_coords, radius, p_detection_max, inf_prop, steepness, cooldown)

  detection_log <- vector("list", steps)

  ##############################################################
  #                   Execute Simulation                       #
  ##############################################################
  message("Executing Simulation...")
  # Execute event loop
  for(i in 1:steps) {
    forager$step(env)


    if (!exists(".Random.seed", .GlobalEnv)) {
      stop("CRITICAL ERROR: Random Number Generator state (.Random.seed) not found. 
          Please run set.seed() before starting the simulation loop.")
    }

    old_seed <- .Random.seed
    set.seed(seed + i)
    detections <- check_all_traps(forager, traps, i)
    
    detection_log[[i]] <- detections
    
    pos_x[i+1] <- forager$x
    pos_y[i+1] <- forager$y
    agent_state[i+1] <- forager$state
    # Reset seed to original each loop to avoid variation in paths between simulations based on detection rate
    .Random.seed <<- old_seed
  }
  message("Simulation Complete!")

  end_time <- Sys.time()
  run_duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  ##############################################################
  #                     Post-Processing                        #
  ##############################################################
  message("Processing data...")

  results <- tibble::tibble(
    step = 0:steps,
    x = pos_x,
    y = pos_y,
    state = agent_state
  )

  detection_log <- detection_log[!sapply(detection_log, is.null)]

  if(length(detection_log) > 0) {
    detections_df <- dplyr::bind_rows(unlist(detection_log, recursive = FALSE))
  } else {
    detections_df <- NULL
  }
  
  ##############################################################
  #                   Plotting & Animation                     #
  ##############################################################
  message("Generating plots and/or Animations...")

  sim_plot <- ggplot(results, aes(x = x, y = y)) +
                geom_path(alpha = 0.7) +
                geom_point(data = tail(results, 1), color = "red",  size = 3) +
                geom_circle(data = camera_coords, aes(x0 = x, y0 = y, r = radius), 
                        inherit.aes = FALSE, color = "blue", alpha = 0.1, linetype = "dotted") +
                geom_point(data = camera_coords, aes(x = x, y = y), 
                        shape = 4, size = 3, color = "blue", stroke = 1.5, inherit.aes = FALSE) +
                annotate("point", x = nest_coords[1], y = nest_coords[2] , size = 3, color = "green") +
                coord_fixed(xlim = c(0, 100), ylim = c(0, 100)) +
                theme_minimal()
  
  anim_time <- NULL

  if(create_gif) {
    anim_start <- Sys.time()
    generate_animation(results, camera_coords, detections_df, radius, sim_name)
    anim_end <- Sys.time()
    anim_time <- as.numeric(difftime(anim_end, anim_start, units = "secs"))
  }

  ##############################################################
  #                       Exporting                            #
  ##############################################################

  export <- list(
    metadata = list(
      sim_name = sim_name,
      seed = seed,
      steps = steps,
      date_run = start_time,
      sim_duration = run_duration,
      animation_time = anim_time
    ),
    movement = results, 
    detections = detections_df,
    plot = sim_plot)
  message("**************************************************")
  message("Export Successful!")
  message("**************************************************")
  message("")

  return(export)
}



