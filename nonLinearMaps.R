#' @importFrom rpatrec noise unbox
#'
#' nonLinearMaps
#'
#' A collection of nonlinear systems exhibiting a variety of interesting dynamical
#' behaviors - chaos, nonlinear periodicity, and nonlinear stochasticity. 
#' Includes options for adding noise to the generated time series.
#' 
#' All functions have a common structure: system(n, noise.level, a, b ...)
#' 
#' where:
#' n - number of time-points to generate
#' noise.level - level of white noise to be added to the final series
#' a, b ... - system-specific parameters
#' 
#' Cubic map
#' @description
#' Generates a time series using a quasiperiodically forced cubic map, exhibiting
#' periodic (one-frequency torus) dynamics by default.
#' @details
#' The cubic map described by Venkatesan & Lakshmanan (2001) is defined as
#' follows:
#' \deqn{x_{i+1}=Q+f \cos \left(2 \pi \theta_{i}\right)-A x_{i}+x_{i}^{3}}
#' \deqn{\theta_{i+1}=\theta_{i}+\omega(\bmod 1)}#' 
#' Depending on the choice of parameters, the cubic map can exhibit a wide range
#' of periodic, chaotic, and strange non-chaotic dynamics. The default parameter
#' settings, \emph{f} = 0, \emph{Q} = 0, and \emph{A} = 1, produce periodic 
#' (one-frequency torus) dynamics. Chaotic dynamics can be produced by setting 
#' \emph{f} = -0.8, \emph{Q} = 0, and \emph{A} = 1.5. 
#' #' The implementation of the function is after Toker et al. (2020).
#' @param n Length of the generated time series. Default: 5000 samples.
#' @param noise.level Defines the variance of white observational noise that will
#' be added to the generated time series. Default: 0.
#' @param f The \emph{f} parameter. Default: 0.
#' @param Q The \emph{Q} parameter. Default: 0.
#' @param A The \emph{A} parameter. Default: 1.
#' @param x0 Optional initial value of \emph{x}. If left unspecified, 
#' it is generated randomly.
#' @param theta0 Optional initial value of \theta. If left unspecified, 
#' it is generated randomly.
#' @return Three vectors containing the values of \emph{x}, \theta, and a linear 
#' combination of the two, \emph{series}.
#' @note Setting the parameters to \emph{f} = 0.7, \emph{Q} = 0, 
#' and \emph{A} = 1.88697 is known to produce a strange non-chaotic regime via the
#' Heagy-Hammel route. Parameters \emph{f} = 0.35, \emph{Q} = 0, and \emph{A} = 
#' 0.35 produce a strange non-chaotic regime via Type-3 intermittency. Setting 
#' \emph{f} = -0.18, \emph{Q} = 0, and \emph{A} = 1.1 produces period-doubled 
#' dynamics.
#' @references Toker, D., F. T. Sommer  and P. R. Ruffino (2020). A simple 
#' method for detecting chaos in nature. Communications Biology, 3, 11.
#' Venkatesan, A. & M. Lakshmanan (2001) Interruption of torus 
#' doubling bifurcation and genesis of strange nonchaotic attractors in a 
#' quasiperiodically forced map: Mechanisms and their characterizations. 
#' Physical Review E, 63, 026219.
#' @author Erik Tihelka
#' @seealso \code{\link{freitasMap}, \link{henonMap}, \link{logisticMap},
#' \link{poincareOscillator}, \link{sineMap}
#' @examples
#' \dontrun{
#' cubicMap(n=1000, noise.level=10, f=-0.8, Q=0, A=1.5)
#' }
#' @export cubicMap
cubicMap = function(n=5000, noise.level=0, f=0, Q=0, A=1, 
                    x0=abs(runif(n=1, min=0, max=1)), 
                    theta0=abs(runif(n=1, min=0, max=1))) {
  n = n+100
  w = GoldenRatio=(1+sqrt(5))/2
  x = theta = series = vector(mode = "numeric", length = n)
  x[[1]] = x0
  theta[[1]] = theta0

  # compute x and theta
  sampleVector = seq(2, n)
  for (i in sampleVector) {
    x[[i]] = Q + f*cos(2*pi*theta[[i - 1]]) - A*x[[i - 1]] + x[[i - 1]]^3
    theta[[i]] = (theta[[i - 1]]+w) %% 1
  }
  x = x[101:n]
  theta = theta[101:n]
  
  # take linear combination of x and theta
  series = x/6+theta/10

  # add noise
  if (noise.level>0) {
    x <- rpatrec::noise(x, type = "white", final_level = noise.level)
    theta <- rpatrec::noise(theta, type = "white", final_level = noise.level)
    series <- rpatrec::noise(series, type = "white", final_level = noise.level)
  }
  else {x <- x
  theta <- theta
  series <- series
  }
  list(series)
}
#' 
#' Freitas map
#' @description
#' Generates a time series using the Freitas map, a nonlinear stochastic system.
#' @details
#' The nonlinear correlated noise system described by Freitas et al. (2009) is as
#' follows:
#' \deqn{v_{i+1}=3 v_{i}+4 v_{i-1}\left(1-v_{i}\right)}
#' where \emph{v_{n}} is a uniform random variable distributed between 0 and 1.
#' The authors have shown that this form of noise is missclassified as chaos by
#' some methods.
#' The implementation of the function is after Toker et al. (2020).
#' @param n Length of the generated time series. Default: 5000 samples.
#' @param noise.level Defines the variance of white observational noise that will
#' be added to the generated time series. Default: 0.
#' @param y0 Optional 2-dimensional vector indicating the starting values of 
#' \emph{y}. If left unspecified, it is generated randomly.
#' @param v0 Optional vector of length \emph{n} indicating the starting values
#' of  \emph{v}. If left unspecified, it is generated randomly.
#' @return A vector containing the values of the time series that has been 
#' generated.
#' @references Freitas, U. S., C. Letellier and L. A. Aguirre, L. A. (2009).
#' Failure in distinguishing colored noise from chaos using the "noise titration"
#' technique. Physical Review E, 79(3), 035201.
#' Toker, D., F. T. Sommer  and P. R. Ruffino (2020). A simple method for 
#' detecting chaos in nature. Communications Biology, 3, 11.
#' @author Erik Tihelka
#' @seealso \code{\link{cubicMap}, \link{henonMap}, \link{logisticMap}, 
#' \link{poincareOscillator}, \link{sineMap}
#' @examples
#' \dontrun{
#' freitasMap(n=1000, noise.level=10)
#' }
#' @export freitasMap
freitasMap = function(n=5000, noise.level=0, y0=runif(n=2, min=0, max=1), 
                  v0=runif(n=n, min=0, max=1)) {
  n = n+2
  y = series = vector(mode="numeric", length=n)
  y[1:2] = y0
  v = v0
  
  sampleVector = seq(3, n)
  for (i in sampleVector) {
    y[[i]] = 3*v[[i-1]]+4*v[[i-2]]*(1-v[[i-1]])
  }
  y <- y[3:n]
  
  # add noise
  if (noise.level>0) {
    series <- rpatrec::noise(y, type = "white", final_level = noise.level)
  }
  else series <- y
  list(series)
}
#' 
#' Hénon map
#' @description
#' A classical two-dimensional dynamical system exhibiting. Exhibits chaotic 
#' behavior by default.
#' @details
#' The Hénon map , sometimes also known as the Hénon-Pomeau map, is defined as 
#' follows:
#' \deqn{ x_n = 1 - a \cdot x_{n - 1}^2 + y_{n - 1}}{x[n] = 1 - a*x[n - 1]^2 + y[n - 1]}
#' \deqn{ y_n = b \cdot x_{n - 1}}{y[n] = b*x[n - 1].}
#' The default selection for parameters, (\emph{a}=1.4 and \emph{b}=0.3), is
#' known to produce chaotic behaviour.
#' The implementation of the function is after Toker et al. (2020) and 
#' Garcia (2022).
#' #' @param n Length of the generated time series. Default: 5000 samples.
#' @param noise.level Defines the variance of white observational noise that will
#' be added to the generated time series. Default: 0.
#' @param a The \emph{a} parameter. Default: 1.4
#' @param b The \emph{b} parameter. Default: 0.3
#' @param x0 Optional initial value of \emph{x}. If left unspecified, 
#' it is generated randomly.
#' @param y0 Optional initial value of \emph{y}. If left unspecified, 
#' it is generated randomly.
#' @param n.transient Number of transient samples to be discarded. 
#' Default: 20 samples.
#' @return Three vectors containing the values of \emph{a}, \emph{b}, and a linear 
#' combination of the two, \emph{series}
#' @note Setting (\emph{a}=1.25 and \emph{b}=0.3) is known to produce 
#' a periodic system. a chaotic system. Some initial values will tend to 
#' @references Garcia, C. A. (2022). nonlinearTseries: Nonlinear Time Series
#' Analysis. R package version 0.2.12.
#' https://CRAN.R-project.org/package=nonlinearTseries
#' Hénon, M. (1976). A two-dimensonal mapping with a strange attractor. 
#' Communications in Mathematical Physics, 50, 376-392.
#' Toker, D., F. T. Sommer  and P. R. Ruffino (2020). A simple method for 
#' detecting chaos in nature. Communications Biology, 3, 11.
#' @author Erik Tihelka
#' @seealso \code{\link{cubicMap}, \link{freitasMap}, \link{logisticMap},
#' \link{poincareOscillator}, \link{sineMap}
#' @examples
#' \dontrun{
#' henonMap(n=1000, noise.level=10, n.transient=100)
#' }
#' @export henonMap
henonMap = function(n=5000, noise.level=0, a=1.4, b=0.3, 
                       x0=runif(n=1, min=-0.5, max=0.5), 
                       y0=runif(n=1, min=-0.5, max=0.5), n.transient=20) {
  n = as.numeric(n + n.transient)
  x = y = series = vector(mode = "numeric", length = n)
  x[[1]] = x0
  y[[1]] = y0
  
  sampleVector = seq(2, n)
  for (i in sampleVector) {
    x[[i]] = y[[i - 1]] + 1 - a * x[[i - 1]] ^ 2
    y[[i]] = b * x[[i - 1]]
  }
  
  # eliminate transients
  transientVector = seq(n.transient)
  x = x[-transientVector]
  y = y[-transientVector]
  
  # take linear combination of x and y
  series = x+y
  
  # add noise
  if (noise.level>0) {
    x <- rpatrec::noise(x, type = "white", final_level = noise.level)
    y <- rpatrec::noise(y, type = "white", final_level = noise.level)
    series <- rpatrec::noise(series, type = "white", final_level = noise.level)
  }
  else {x <- x
  y <- y
  series <- series
  }
  list(series)
}
#' 
#' Logistic map
#' @description
#' Generates a time series using the logistic map, exhibiting chaotic dynamics 
#' by default.
#' @details
#' The logistic map is defined as follows:
#' \deqn{x_n = r  \cdot  x_{n-1}   \cdot  (1 - x_{n-1})}{x[n] = r * x[n-1]  * (1 - x[n-1])}
#' The default selection for the \emph{r} parameter is known to
#' produce a chaotic time series.
#' The implementation of the function is after Toker et al. (2020) and 
#' Garcia (2022).
#' @param n Length of the generated time series. Default: 5000 samples.
#' @param noise.level Defines the variance of white observational noise that will
#' be added to the generated time series. Default: 0.
#' @param r The \emph{r} parameter. Default: 3.9
#' @param x0 Optional initial value of \emph{r}. If left unspecified, 
#' it is generated randomly.
#' @param n.transient Number of transient samples to be discarded. 
#' Default: 20 samples.
#' @return A vector containing the values of the time series that has been 
#' generated.
#' @note \emph{r} = 3.9 is known to produce a chaotic system.
#' Some initial values may lead to an unstable system that will tend to 
#' infinity.
#' @references Garcia, C. A. (2022). nonlinearTseries: Nonlinear Time Series
#' Analysis. R package version 0.2.12.
#' https://CRAN.R-project.org/package=nonlinearTseries
#' May, R. M. (1976). Simple mathematical models with very complicated dynamics.
#' Nature, 261, 459-467.
#' Toker, D., F. T. Sommer  and P. R. Ruffino (2020). A simple method for 
#' detecting chaos in nature. Communications Biology, 3, 11.
#' @author Erik Tihelka
#' @seealso \code{\link{cubicMap}, \link{freitasMap}, \link{henonMap},
#' \link{poincareOscillator}, \link{sineMap}
#' @examples
#' \dontrun{
#' logisticMap(n=1000, noise.level=10, n.transient=10)
#' }
#' @export logisticMap
logisticMap = function(n=5000, noise.level=0, r=3.9, x0=runif(n = 1, min = 0,
                                                                  max = 1), 
                          n.transient=20) {
  n = n + n.transient
  x = vector(mode = "numeric", length = n)
  x[[1]] = x0
  sampleVector = seq(2, n)
  for (i in sampleVector) {
    x[[i]] = r * x[[i - 1]]  * (1 - x[[i - 1]])
  }

  # eliminate transients
  transientVector = seq(n.transient)
  x = x[-transientVector]
  
  # add noise
  if (noise.level>0) {
    series <- rpatrec::noise(x, type = "white", final_level = noise.level)
  }
  else series <- x
  list(series)
}

#' Poincaré oscillator
#' @description Generates a time series using the Poincaré oscillator, 
#' exhibiting nonlinear periodic behaviour by default.
#' @details
#' The Poincaré oscillator is a classical model system, often used to study 
#' biological oscillators subjected to periodic stimulation. It can exhibit a
#' wide range of dynamical behaviours, including nonlonear periodicity,
#' quasi-periodicity, and chaos. The oscillator is described as follows:
#' \deqn{\frac{d r}{d t}=k r(1-r)}
#' \deqn{\frac{d \Phi}{d t}=2 \Phi}
#' where \emph{k} t controls the oscillator's relaxation rate. The phase of the
#' system is described by \Phi, its angular coordinate in a unit cycle. 
#' Perturbation of magnitude \emph{b} lead the system away from the unit cycle, 
#' resetting the phase of the oscillator as follows:
#' \deqn{g(\phi)=\frac{1}{2 \pi} \arccos \frac{\cos 2 \pi \phi+b}{\sqrt{1+b^{2}+2 b \cos 2 \pi \phi}}(\bmod 1)}
#' The period of stimulation is determined by parameter \tau.
#' The implementation of the function is after Wanzhen et al. (1992), Glass 
#' (2003), and Toker et al. (2020).
#' Copyright Leon Glass 2003, Centre for Nonlinear Dynamics in Physiology and Medicine.
#' @param n Length of the generated time series. Default: 5000 samples.
#' @param noise.level Defines the variance of white observational noise that will
#' be added to the generated time series. Default: 0.
#' @param b The \emph{b} parameter, representing stimulation strength. Default: 
#' 1.113
#' @param tau The \tau parameter, representing the period of the stimulation. 
#' Default: 0.69
#' @return 'Series', a vector containing the values of the time series that has been 
#' generated. A second vector, 'beats', lists the number of beats between 
#' successive stimuli.
#' @note Parameter values \emph{b} = 1.13 and \tau = 0.65 are known to produce a
#' chaotic system. Quasiperiodic dynamics are produced by \emph{b} = 0.95 and 
#' \tau = 0.75.
#' @references Glass, L. (2003) Resetting and entraining biological rhythms. 
#' in Nonlinear Dynamics in Physiology and Medicine (pp. 123-148). Springer, 
#' New York.
#' Toker, D., F. T. Sommer  and P. R. Ruffino (2020). A simple method for 
#' detecting chaos in nature. Communications Biology, 3, 11.
#' Wanzhen, Z., L. Glass and A. Shrier (1992) The topology of phase response 
#' curves induced by single and paired stimuli in spontaneously oscillating 
#' chick heart cell aggregates. Journal of Biological Rhythms, 7, 89-104.
#' @author Erik Tihelka
#' @seealso \code{\link{cubicMap}, \link{freitasMap}, \link{henonMap},
#' \link{logisticMap}, \link{sineMap}
#' @examples
#' \dontrun{
#' poincareOscillator(n=1000, noise.level=10)
#' }
#' #' @export poincareOscillator
poincareOscillator = function(n=5000, noise.level=0, b=1.113, tau=0.69) {
  n = n+100
  phi = beats = vector(mode="numeric", length=n)
  phi <- rep(0, n)
  phi[[1]] <- runif(1)
  
  sampleVector = seq(2, n)
  for (i in sampleVector) {
    angle = 2*pi*phi[[i-1]]
    rprime=sqrt(1+b^2+2*b*cos(angle))
    argument=(cos(angle)+b)/rprime
    phi[[i]]=acos(argument)/(2*pi)
    if (phi[[i-1]] > 0.5) {
      phi[[i]]=1-phi[[i]]
    }
    phi[[i]]=phi[[i]]+tau
    beats[[i]]=phi[[i]]-phi[[i]] %% 1
    phi[[i]]=phi[[i]] %% 1
  }
  phi = phi[101:n]
  beats = beats[101:n]
  
  # add noise
  if (noise.level>0) {
    series <- rpatrec::noise(phi, type = "white", final_level = noise.level)
    beats <- rpatrec::noise(beats, type = "white", final_level = noise.level)
  }
  else {series <- phi
  beats <- beats
  }
  list(series)
}
#' 
#' Sine map
#' @description
#' Generates a time series using a noise-driven sine map described
#' by Freitas et al. (2009).
#' @details
#' The nonlinear stochastic system described by Freitas et al. (2009) is as
#' follows:
#' \deqn{x_{i+1}=\mu \sin \left(x_{i}\right)+Y_{i} \eta_{i}}
#' where \mu = 2.4, \emph{Y_{i}} is a random Bernoulli process that equals 1 
#' with probability 0.01 and 0 with probability 0.99, and \eta_{i} is a random
#' variable with a uniform distribution between ???2 and 2. Freitas et al. (2009)
#' have shown that this map is sometimes missclassfied as chaotic by some 
#' methods.
#' The implementation of the function is after Toker et al. (2020).
#' @param n Length of the generated time series. Default: 5000 samples.
#' @param noise.level Defines the variance of white observational noise that will
#' be added to the generated time series. Default: 0.
#' @param mu The \mu parameter. Default: 2.4.
#' @param x0 Optional initial value of \emph{x}. If left unspecified, 
#' it is generated randomly.
#' @return A vector containing the values of the time series that has been 
#' generated.
#' @references Freitas, U. S., C. Letellier and L. A. Aguirre, L. A. (2009).
#' Failure in distinguishing colored noise from chaos using the "noise titration"
#' technique. Physical Review E, 79(3), 035201.
#' Toker, D., F. T. Sommer  and P. R. Ruffino (2020). A simple method for 
#' detecting chaos in nature. Communications Biology, 3, 11.
#' @author Erik Tihelka
#' @seealso \code{\link{cubicMap}, \link{freitasMap}, \link{henonMap},
#' \link{logisticMap}, \link{poincareOscillator} 
#' @examples
#' \dontrun{
#' sineMap(n=1000, noise.level=10)
#' }
sineMap = function(n=5000, noise.level=0, mu=2.4, x0=runif(1)) {
  x = vector(mode="numeric", length=n)
  x[[1]] = x0
  
  sampleVector = seq(2, n)
  for (i in sampleVector) {
    # random Bernoulli process
    q=runif(1)
    if (q<0.01) {
      Y=1
    }
    else Y=0
    x[[i]]=mu*sin(x[[i-1]])+Y*(4*runif(1)-2)
  }
  
  # add noise
  if (noise.level>0) {
    series <- rpatrec::noise(x, type = "white", final_level = noise.level)
  }
  else series <- x
  list(series)
}