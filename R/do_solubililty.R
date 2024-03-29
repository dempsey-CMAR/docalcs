#' Calculate the solubility of dissolved oxygen
#'
#' @description Calculate the solubility of dissolved oxygen based on
#'   temperature, salinity and barometric pressure, following the equations in
#'   Benson and Krause 1984 or Garcia and Gordon 1992.
#'
#' @details Solubility of dissolved oxygen (mg/L) is calculated using the
#'   equations in Benson and Krause 1984, or the or Garcia and Gordon 1992
#'   refit. These equations should only be used when 0 < Temperature < 40
#'   degrees Celcius, 0 < Salinity < 40 PSU and 0.5 < Pressure < 1.1 atm.
#'
#'   Results from this function should closely match those from the USGS
#'   DOTABLES Table B (output from GG method may differ in the second decimal
#'   place) at \url{https://water.usgs.gov/water-resources/software/DOTABLES/}.
#'
#'   For more information see Equations 24 and 32, and Table 2 from Benson and
#'   Krause 1984 or Equation 8 and Table 1 from Garcia and Gordon 1992.
#'
#'   For the Garcia and Gordon equation, coefficients from the first column of
#'   Table 1 were used. Conversion factor of 1.42905 was applied to convert from
#'   cm^3/dm^3 (mL/L) to mg/L (USGS 2011).
#'
#'   \bold{References}
#'
#'   Benson, Bruce B., Krause, Daniel, (1984), The concentration and isotopic
#'   fractionation of oxygen dissolved in freshwater and seawater in equilibrium
#'   with the atmosphere, Limnology and Oceanography, 3, doi:
#'   10.4319/lo.1984.29.3.0620.
#'
#'   Garcia, H., and L. Gordon (1992), \emph{Oxygen solubility in seawater:
#'   Better fitting equations}, Limnol. Oceanogr., 37(6).
#'
#'   USGS. \emph{Change to Solubility Equations for Oxygen in Water.} Technical
#'   Memorandum 2011.03. USGS Office of Water Quality, 2011.
#'
#' @inheritParams do_salinity_correction
#' @inheritParams do_pressure_correction
#'
#' @param dat Data frame with at least one column: \code{temperature_degree_c}.
#'   Corresponding salinity (psu) and pressure (atm) data may be included in
#'   columns \code{salinity_psu} and \code{pressure_atm}. Additional columns will be
#'   ignored and returned.
#'
#' @param method Equation to use to calculate dissolved oxygen solubility.
#'   Options are \code{method = "garcia-gordon"} (the default) and \code{method
#'   = "benson-krause"}.
#'
#' @param return_factors Logical parameter. If \code{TRUE} the function returns
#'   \code{dat.wide} with additional columns for the parameters used to
#'   calculate \code{F_s} (salinity correction factor), \code{F_p} (pressure
#'   correction factor), and \code{C_p} (DO solubility).
#'
#'   If \code{FALSE}, the function returns \code{dat.wide} with additional
#'   column(s) \code{salinity_psu}, \code{pressure_atm}, and \code{C_p}.
#'
#' @importFrom dplyr %>% mutate select
#'
#' @export



do_solubility <- function(dat,
                         sal = NULL, p_atm = NULL,
                         method = "garcia-gordon",
                         return_factors = FALSE){

# Benson Krause ------------------------------------------------------------

  if(tolower(method) == "benson-krause"){

    dat_out <- dat %>%
      mutate(
        temperature_degree_c = as.numeric(temperature_degree_c),

        # temperature in Kelvin
        T_Kelvin = temperature_degree_c + 273.15,

        # Equation 32 (modified to return units of mg / L)
        C_star = exp(
          -139.34411 +
            1.575701e5 / T_Kelvin -
            6.642308e7 / T_Kelvin^2 +
            1.243800e10 / T_Kelvin^3 -
            8.621949e11 / T_Kelvin^4
          # this term is accounted for in DO_salinity_correction()
          # -Salinity * (0.017674 - 10.754 / T_Kelvin + 2140.7 / T_Kelvin^2)
        )
      ) %>%
      do_salinity_correction(sal = sal, method = "benson-krause") %>%
      do_pressure_correction(p_atm = p_atm) %>%
      mutate(
        C_p = C_star * F_s * F_p,
        C_p = round(C_p, digits = 2)
      )

  }


# Garcia Gordon -----------------------------------------------------------

  if(tolower(method) == "garcia-gordon"){

    # coefficients
    A0_GG <- 2.00907
    A1_GG <- 3.22014
    A2_GG <- 4.05010
    A3_GG <- 4.94457
    A4_GG <- -2.56847E-1
    A5_GG <- 3.88767

    mg_L <- 1.42905   # to convert from mL/L to mg/L

    dat_out <- dat %>%
      mutate(
        temperature_degree_c = as.numeric(temperature_degree_c),

        # scaled temperature
        T_s = log((298.15 - temperature_degree_c) / (273.15 + temperature_degree_c)),

        # oxygen solubility (equation 8, coefficients from first col Table 1)
        C_star = mg_L * exp(
          A0_GG
          + A1_GG * T_s
          + A2_GG * T_s^2
          + A3_GG * T_s^3
          + A4_GG * T_s^4
          + A5_GG * T_s^5
          # these terms are accounted for in DO_salinity_correction()
         # + Salinity * (B0_GG + B1_GG * T_s + B2_GG * T_s^2 + B3_GG * T_s^3)
         # + C0_GG * Salinity^2
        )
      ) %>%
      do_salinity_correction(sal = sal, method = "garcia-gordon") %>%
      do_pressure_correction(p_atm = p_atm) %>%
      mutate(
        C_p = C_star * F_s * F_p,
        C_p = round(C_p, digits = 2)
      )

  }

  if(isFALSE(return_factors)){

    dat_out <- dat_out %>%
      select(-C_star, -F_s, -T_Kelvin, -theta, -P_wv, - alt_correction, -F_p)

  }

  dat_out

}
