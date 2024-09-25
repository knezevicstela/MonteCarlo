library(shiny)
library(shinythemes)

# Linearni kongruencijalni generator (LCG)
LCG <- function(seed, n, a = 1664525, c = 1013904223, m = 2^32) {
  random_numbers <- numeric(n)
  random_numbers[1] <- seed
  
  for (i in 2:n) {
    random_numbers[i] <- (a * random_numbers[i - 1] + c) %% m
  }
  
  random_numbers <- random_numbers / m
  return(random_numbers)
}

# Multiplikativni kongruencijalni generator (MCG)
MCG <- function(seed, n, a = 1664525, m = 2^32) {
  random_numbers <- numeric(n)
  random_numbers[1] <- seed
  
  for (i in 2:n) {
    random_numbers[i] <- (a * random_numbers[i - 1]) %% m
  }
  
  random_numbers <- random_numbers / m
  return(random_numbers)
}

# Mersenne Twister generator
MersenneTwister <- function(seed, n) {
  set.seed(seed, kind = "Mersenne-Twister")
  return(runif(n))
}

# Haltonov sekvencijalni generator
Halton <- function(n, base = 2) {
  seq <- numeric(n)
  for (i in 1:n) {
    seq[i] <- 0
    f <- 1 / base
    j <- i
    while (j > 0) {
      seq[i] <- seq[i] + f * (j %% base)
      j <- floor(j / base)
      f <- f / base
    }
  }
  return(seq)
}

# Uniformni generator slučajnih brojeva
UniformGenerator <- function(seed, n) {
  set.seed(seed)
  return(runif(n))
}

# Eksponencijalni generator slučajnih brojeva
ExponentialGenerator <- function(seed, n) {
  set.seed(seed)
  return(rexp(n))
}

# Normalni generator slučajnih brojeva
NormalGenerator <- function(seed, n) {
  set.seed(seed)
  return(rnorm(n))
}

# Klasična Monte Carlo integracija
MonteCarloIntegracija <- function(f, a, b, n, generator, seed) {
  if (is.infinite(a) || is.infinite(b)) {
    # Za beskonačne intervale koristimo normalnu distribuciju
    u <- generator(seed, n)
    x <- qnorm(u)
    integral_approx <- mean(f(x))
  } else {
    u <- generator(seed, n)
    x <- a + (b - a) * u
    integral_approx <- mean(f(x)) * (b - a)
  }
  
  return(integral_approx)
}


# Stratificirana Monte Carlo integracija
StratificiranaMonteCarlo <- function(f, a, b, n, generator, seed) {
  if (is.infinite(a) || is.infinite(b)) {
    k <- max(1, floor(sqrt(n)))
    n_per_stratum <- floor(n / k)
    u <- generator(seed, n)
    stratumi <- rep(1:k, each = n_per_stratum)
    x <- qnorm(u)
    integral_approx <- mean(f(x))
  } else {
    k <- max(1, floor(sqrt(n)))
    n_per_stratum <- floor(n / k)
    u <- generator(seed, n)
    stratumi <- rep(1:k, each = n_per_stratum)
    x <- a + (stratumi - 1 + u) * (b - a) / k
    integral_approx <- (b - a) * mean(f(x))
  }
  
  return(integral_approx)
}


# Antitetička Monte Carlo integracija
AntitetickaMonteCarlo <- function(f, a, b, n, generator, seed) {
  if (is.infinite(a) || is.infinite(b)) {
    u <- generator(seed, n/2)
    x1 <- qnorm(u)
    x2 <- qnorm(1 - u)
    integral_approx <- mean(f(c(x1, x2)))
  } else {
    u <- generator(seed, n/2)
    x1 <- a + (b - a) * u
    x2 <- a + (b - a) * (1 - u)
    integral_approx <- (b - a) * mean(f(c(x1, x2)))
  }
  
  return(integral_approx)
}


# Hit-or-Miss Monte Carlo integracija
HitOrMissMonteCarlo <- function(f, a, b, n, generator, seed) {
  if (is.infinite(a) || is.infinite(b)) {
    u_x <- generator(seed, n)
    u_y <- generator(seed + 1, n)
    x <- qnorm(u_x)
    f_max <- max(f(x))  # Procjena maksimuma funkcije
    y <- f_max * u_y
    hits <- sum(y <= f(x))
    area <- f_max
    integral_approx <- area * (hits / n)
  } else {
    u_x <- generator(seed, n)
    u_y <- generator(seed + 1, n)
    x <- a + (b - a) * u_x
    f_max <- max(f(x))  # Pretpostavka da imamo procjenu maksimuma funkcije
    y <- f_max * u_y
    hits <- sum(y <= f(x))
    area <- (b - a) * f_max
    integral_approx <- area * (hits / n)
  }
  
  return(integral_approx)
}


# Importance Sampling Monte Carlo integracija
ImportanceSamplingMonteCarlo <- function(f, a, b, n, generator, seed) {
  if (is.infinite(a) || is.infinite(b)) {
    u <- generator(seed, n)
    x <- qnorm(u)
    g <- dnorm(x)
    w <- 1 / g
    integral_approx <- mean(f(x) * w)
  } else {
    u <- generator(seed, n)
    x <- a + (b - a) * u
    g <- dnorm(x, mean = (a + b) / 2, sd = (b - a) / 6)
    w <- 1 / g
    integral_approx <- mean(f(x) * w)
  }
  
  return(integral_approx)
}


# Random Walk Monte Carlo integracija
RandomWalkMonteCarlo <- function(f, a, b, n, generator, seed) {
  if (is.infinite(a) || is.infinite(b)) {
    u <- generator(seed, n)
    x <- numeric(n)
    x[1] <- qnorm(u[1])
    for (i in 2:n) {
      x[i] <- x[i - 1] + rnorm(1)
    }
    integral_approx <- mean(f(x))
  } else {
    u <- generator(seed, n)
    x <- numeric(n)
    x[1] <- a + (b - a) * u[1]
    for (i in 2:n) {
      x[i] <- x[i - 1] + (b - a) * (u[i] - 0.5)
      if (x[i] < a) x[i] <- a
      if (x[i] > b) x[i] <- b
    }
    integral_approx <- (b - a) * mean(f(x))
  }
  
  return(integral_approx)
}


# Numerička integracija
NumerickaIntegracija <- function(f, a, b) {
  if (is.infinite(a) || is.infinite(b)) {
    rezultat <- integrate(f, a, b, rel.tol = .Machine$double.eps^0.5)
  } else {
    rezultat <- integrate(f, a, b)
  }
  return(rezultat$value)
}


# Funkcija za izračunavanje aproksimacije za različite vrijednosti n i prikazivanje rezultata
compute_examples <- function(f, a, b, seed, generator) {
  ns <- c(10, 100, 1000, 10000, 100000)
  numeric_result <- NumerickaIntegracija(f, a, b)
  
  results <- data.frame(
    N = ns,
    Monte_Carlo = numeric(length(ns)),
    Numericka_Integracija = numeric(length(ns)),  # Dodano
    Apsolutna_Pogreska = numeric(length(ns)),
    Relativna_Pogreska_Postotak = numeric(length(ns))
  )
  
  for (i in 1:length(ns)) {
    mc_result <- NA
    if (is.infinite(a) || is.infinite(b)) {
      mc_result <- ImportanceSamplingMonteCarlo(f, a, b, ns[i], generator, seed + i)
    } else {
      mc_result <- MonteCarloIntegracija(f, a, b, ns[i], generator, seed + i)
    }
    
    abs_error <- abs(mc_result - numeric_result)
    rel_error <- abs_error / abs(numeric_result) * 100
    
    results$Monte_Carlo[i] <- mc_result
    results$Numericka_Integracija[i] <- numeric_result  # Dodano
    results$Apsolutna_Pogreska[i] <- abs_error
    results$Relativna_Pogreska_Postotak[i] <- rel_error
  }
  
  return(results)
}

# UI dio Shiny aplikacije
ui <- fluidPage(
  theme = shinytheme("flatly"),
  titlePanel("Monte Carlo Integracija"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel(
          "Unos",
          textInput("funkcija", "Unesite funkciju za integraciju:", value = "sin(x)"),
          
          # Promjena granica s mogućnošću odabira beskonačnosti
          selectizeInput("a", "Donja granica integracije (a):",
                         choices = c("-∞" = "-Inf", 
                                     "0" = "0", 
                                     "1" = "1", 
                                     "2" = "2"),
                         selected = "0",
                         options = list(create = TRUE)), # Omogućava unos vlastitih vrijednosti
          
          selectizeInput("b", "Gornja granica integracije (b):",
                         choices = c("+∞" = "Inf", 
                                     "pi" = "3.14159", 
                                     "1" = "1", 
                                     "2" = "2"),
                         selected = "3.14159",
                         options = list(create = TRUE)), # Omogućava unos vlastitih vrijednosti
          
          numericInput("n", "Broj slučajnih točaka (n):", value = 1000, min = 1),
          numericInput("seed", "Seed za generator slučajnih brojeva:", value = 1234),
          selectInput("generator", "Odaberite generator slučajnih brojeva:",
                      choices = list("Linearni kongruencijalni generator (LCG)" = "LCG",
                                     "Multiplikativni kongruencijalni generator (MCG)" = "MCG",
                                     "Mersenne Twister generator" = "MersenneTwister",
                                     "Haltonov sekvencijalni generator" = "Halton",
                                     "Uniformni generator" = "Uniform",
                                     "Eksponencijalni generator" = "Exponential",
                                     "Normalni generator" = "Normal")),
          selectInput("metoda", "Odaberite Monte Carlo metodu:",
                      choices = list("Klasična Monte Carlo integracija" = "MC",
                                     "Stratificirana Monte Carlo integracija" = "StratMC",
                                     "Antitetička Monte Carlo integracija" = "AntitetickaMC",
                                     "Hit-or-Miss Monte Carlo integracija" = "HitOrMissMC",
                                     "Importance Sampling Monte Carlo integracija" = "ImportanceMC",
                                     "Random Walk Monte Carlo integracija" = "RandomWalkMC")),
          actionButton("calculate", "Izračunaj"),
          style = "padding: 20px;"  # Stil za povećanje razmaka oko elemenata
        ),
        tabPanel(
          "Primjeri Integrala",
          selectInput("example", "Odaberite primjer integrala:",
                      choices = list(
                        "sin(x) od 0 do π" = "sin",
                        "exp(-x^2) od -∞ do +∞" = "exp",
                        "1/(1+x^2) od -∞ do +∞ (integral arctangenta)" = "arctan",
                        "cos(x) od 0 do π/2" = "cos"
                      )),
          uiOutput("example_results"),
          style = "padding: 20px;"  # Stil za povećanje razmaka oko elemenata
        )
      )
    ),
    mainPanel(
      uiOutput("rezultat_latex"),
      verbatimTextOutput("rezultat"),
      style = "padding: 20px;"  # Stil za povećanje razmaka oko elemenata
    )
  )
)

# Server dio Shiny aplikacije
server <- function(input, output, session) {
  results <- reactiveValues(
    mc_result = NULL,
    num_result = NULL,
    abs_error = NULL,
    rel_error = NULL
  )
  
  observeEvent(input$calculate, {
    f <- function(x) eval(parse(text = input$funkcija))
    a <- as.numeric(input$a)  # Pretvorba unosa u numeričke vrijednosti, uz obradu beskonačnosti
    b <- as.numeric(input$b)  
    n <- input$n
    seed <- input$seed
    generator <- switch(input$generator, 
                        LCG = LCG, 
                        MCG = MCG, 
                        MersenneTwister = MersenneTwister, 
                        Halton = function(seed, n) Halton(n),
                        Uniform = function(seed, n) UniformGenerator(seed, n),
                        Exponential = function(seed, n) ExponentialGenerator(seed, n),
                        Normal = function(seed, n) NormalGenerator(seed, n))
    
    rezultat_mc <- switch(input$metoda,
                          MC = MonteCarloIntegracija(f, a, b, n, generator, seed),
                          StratMC = StratificiranaMonteCarlo(f, a, b, n, generator, seed),
                          AntitetickaMC = AntitetickaMonteCarlo(f, a, b, n, generator, seed),
                          HitOrMissMC = HitOrMissMonteCarlo(f, a, b, n, generator, seed),
                          ImportanceMC = ImportanceSamplingMonteCarlo(f, a, b, n, generator, seed),
                          RandomWalkMC = RandomWalkMonteCarlo(f, a, b, n, generator, seed))
    
    rezultat_num <- NumerickaIntegracija(f, a, b)
    
    results$mc_result <- round(rezultat_mc, 5)
    results$num_result <- round(rezultat_num, 5)
    results$abs_error <- round(abs(rezultat_mc - rezultat_num), 5)
    results$rel_error <- round(results$abs_error / abs(rezultat_num) * 100, 5)
    
    output$rezultat_latex <- renderUI({
      withMathJax(
        paste0(
          "$$ \\int_{", a, "}^{", b, "} ", gsub("x", "x", input$funkcija), " \\, dx $$ ")
      )
    })
    
    output$rezultat <- renderText({
      paste("Monte Carlo aproksimacija integrala je:", results$mc_result, "\n",
            "Numerička integracija daje rezultat:", results$num_result, "\n",
            "Apsolutna pogreška je:", results$abs_error, "\n",
            "Relativna pogreška (%):", results$rel_error)
    })
  })
  
  observeEvent(input$tabs, {
    if (input$tabs == "Primjeri Integrala") {
      
      results$mc_result <- NULL
      results$num_result <- NULL
      results$abs_error <- NULL
      results$rel_error <- NULL
      output$rezultat_latex <- renderUI(NULL)
      output$rezultat <- renderText(NULL)
    }
  })
  
  output$example_results <- renderUI({
    example_choice <- input$example
    seed <- 1234  
    generator <- UniformGenerator  # za primjere koristimo uniformni generator
    
    if (example_choice == "sin") {
      f <- function(x) sin(x)
      a <- 0
      b <- pi
      integral_latex <- "\\int_{0}^{\\pi} \\sin(x) \\, dx"
    } else if (example_choice == "exp") {
      f <- function(x) exp(-x^2)
      a <- -Inf
      b <- Inf
      integral_latex <- "\\int_{-\\infty}^{\\infty} e^{-x^2} \\, dx"
    } else if (example_choice == "arctan") {
      f <- function(x) 1 / (1 + x^2)
      a <- -Inf
      b <- Inf
      integral_latex <- "\\int_{-\\infty}^{\\infty} \\frac{1}{1 + x^2} \\, dx"
    } else if (example_choice == "cos") {
      f <- function(x) cos(x)
      a <- 0
      b <- pi / 2
      integral_latex <- "\\int_{0}^{\\frac{\\pi}{2}} \\cos(x) \\, dx"
    }
    
    results <- compute_examples(f, a, b, seed, generator)
    
    
    results_latex <- paste0(
      "\\begin{array}{|c|c|c|c|c|} \\hline\n",
      " N & \\text{Monte Carlo} & \\text{Numerička Integracija} & \\text{Apsolutna Pogreška} & \\text{Relativna Pogreška} \\\\ \\hline",
      paste0(apply(results, 1, function(row) {
        paste(row[1], "&", round(row[2], 5), "&", round(row[3], 5), "&", round(row[4], 5), "&", round(row[5], 5), "\\\\ \\hline")
      }), collapse = "\n"),
      "\\end{array}"
    )
    
    
    withMathJax(
      paste0(
        "$$ ", integral_latex, " $$",
        results_latex
      )
    )
  })
}

# Pokretanje aplikacije
shinyApp(ui = ui, server = server)
