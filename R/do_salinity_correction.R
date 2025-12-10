#' Salinity correction factor for dissolved oxygen concentration measurements
#'
#' @details Dissolved oxygen concentration measured by HOBO sensors needs to be
#'   corrected for salinity (see manual:
#'   \url{https://www.onsetcomp.com/files/manual_pdfs/15603-B-MAN-U26x.pdf}).
#'
#'   The salinity correction factor can be calculated using the Benson and
#'   Krause 1984 equations or the Garcia and Gordon 1992 equations. These
#'   equations should only be used when 0 < Temperature < 40 degrees Celcius, 0
#'   < Salinity < 40 PSU and 0.5 < Pressure < 1.1 atm.
#'
#'   \bold{Benson & Krause, 1984}:
#'
#'   Used by the USGS DOTables (USGS, 2011).
#'
#'   Results from this equation should match those from the USGS DOTABLES (part
#'   C) at \url{https://water.usgs.gov/water-resources/software/DOTABLES/}.
#'
#'   Final term in Equation 32:
#'
#'   \deqn{F_{s} = exp(-Salinity * (0.017674 - 10.754 / T_{Kelvin} + 2140.7 /
#'   T_{Kelvin}^2))}
#'
#'   \bold{Garcia &  Gordon, 1992:}
#'
#'   Garcia & Gordon re-fit the Benson & Krause data with a higher order
#'   polynomial and defined a scaled temperature (\eqn{T_{s}}).
#'
#'   This correction factor is used in the SCOR WG 142 (Part C).
#'
#'   Results from this equation are very similar to the Benson & Krause
#'   correction factor (to ~4 decimal places).
#'
#'   \deqn{B0 = -6.24523E-3} \deqn{B1 = -7.37614E-3} \deqn{B2 = -1.03410E-2}
#'   \deqn{B3 = -8.17083E-3} \deqn{C0 = -4.88682E-7}
#'
#'   \deqn{F_{s} = exp(Salinity * (B0 + B1 * T_{s} + B2 * T_{s}^2 + B3 *
#'   T_{s}^3) + C0 * Salinity^2)}
#'
#'   \deqn{T_{s} = log((298.15 - Temperature) / (273.15 + Temperature))}
#'
#'   Note: the HOBO Dissolved Oxygen Assistant Software salinity correction
#'   factor uses the same form as the Garcia and Gordon equation, but different
#'   coefficients (although the results are similar).
#'
#'   \bold{References}
#'
#'   Benson, Bruce B., Krause, Daniel, (1984), \emph{The concentration and
#'   isotopic fractionation of oxygen dissolved in freshwater and seawater in
#'   equilibrium with the atmosphere, Limnology and Oceanography}, 3, doi:
#'   10.4319/lo.1984.29.3.0620.
#'
#'   Bittig, H. &.-J. (2016). \emph{SCOR WG 142: Quality Control Procedures for
#'   Oxygen and Other Biogeochemical Sensors on Floats and Gliders.
#'   Recommendations on the conversion between oxygen quantities for Bio-Argo
#'   floats and other autonomous sensor platforms.}
#'   \url{https://repository.oceanbestpractices.org/handle/11329/417}.
#'
#'   Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater:
#'   Better fitting equations}, Limnol. Oceanogr., 37(6).
#'
#'   USGS. \emph{Change to Solubility Equations for Oxygen in Water.} Technical
#'   Memorandum 2011.03. USGS Office of Water Quality, 2011.
#'
#' @param dat_wide Data frame with at least one column:
#'   \code{temperature_degree_c}. Corresponding salinity (psu) data may be
#'   included in column \code{salinity_psu}. Additional columns will be ignored
#'   and returned.
#'
#' @param sal A single value of salinity (psu). This value must be specified if
#'   there is no \code{salinity_psu} column in \code{dat_wide}. Default is \code{Sal
#'   = NULL}. Note: if \code{Sal} is specified when there is a \code{salinity_psu}
#'   column in \code{dat_wide}, the function will stop with an error.
#'
#' @param method Equation to use to calculate salinity correction factor.
#'   Options are \code{method = "garcia-gordon"} (the default) and \code{method
#'   = "benson-krause"}.
#'
#' @return Returns \code{dat_wide} with additional column(s), \code{F_s} and
#'   \code{salinity_psu}.
#'
#' @importFrom dplyr %>% mutate
#'
#' @export


do_salinity_correction <- function(
  dat_wide,
  sal = NULL,
  method = "garcia-gordon"
){

# Error Messages ----------------------------------------------------------
  cols <- colnames(dat_wide)

  if(!(tolower(method) %in% c("garcia-gordon", "benson-krause"))){
    stop("Method argument not recognized.
         \nHINT: method must be 'garcia-gordon' or 'benson-krause'")
  }

  if("salinity_psu" %in% cols & !is.null(sal)){
    stop("Conflicting salinity values.
         \nHINT: Remove column Salinity or set Sal argument to NULL.")
  }

# calculate F_s -----------------------------------------------------------

  if(!is.null(sal)) dat_wide <- mutate(dat_wide, salinity_psu = sal)

# Benson-Krause Method ----------------------------------------------------

  if(tolower(method) == "benson-krause"){

    B0_BK <- 0.017674
    B1_BK <- -10.754
    B2_BK <-  2140.7

    dat_out <- dat_wide %>%
      mutate(
        temperature_degree_c = as.numeric(temperature_degree_c),

        # temperature in Kelvin
        T_Kelvin = temperature_degree_c  + 273.15,

        # correction factor
        F_s = exp(-salinity_psu * (B0_BK + B1_BK / T_Kelvin + B2_BK / T_Kelvin^2)),
        F_s = round(F_s, digits = 4)
      ) %>%
      select(-T_Kelvin)
  }

# Garcia-Gordon Method ----------------------------------------------------

  if(tolower(method) == "garcia-gordon"){

    B0_GG = -6.24523E-3
    B1_GG = -7.37614E-3
    B2_GG = -1.03410E-2
    B3_GG = -8.17083E-3
    C0_GG = -4.88682E-7

    dat_out <- dat_wide %>%
      mutate(
        temperature_degree_c = as.numeric(temperature_degree_c),

        # scaled temperature
        T_s = log((298.15 - temperature_degree_c) / (273.15 + temperature_degree_c)),

        # correction factor
        F_s = exp(salinity_psu * (B0_GG + B1_GG * T_s + B2_GG * T_s^2 + B3_GG * T_s^3)
                    + C0_GG * salinity_psu^2),
        F_s = round(F_s, digits = 4)
      ) %>%
      select(-T_s)

  }

  dat_out

}

