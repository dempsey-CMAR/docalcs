#' Pressure correction factor for dissolved oxygen measurements
#'
#' @details Dissolved oxygen concentration (and partial pressure?) measurements
#'   should be corrected for atmospheric pressure if the pressure deviates from
#'   substantially from 1 atm (Benson & Krause, 1984).
#'
#'   This function calculates pressure correction factor following Benson and
#'   Krause (1984), as suggested by USGS (2011). This is similar to what is
#'   described in the HOBO manual HOBO U26 Percent Saturation Calculation.pdf
#'   (except this function accounts for salinity in the water vapour
#'   calculation).
#'
#'   Equation 24 in Benson & Krause 1984:
#'
#'   \deqn{C_{p} = C_{star} * Pressure * (((1 - P_{wv} / Pressure) * (1 - theta
#'   * Pressure)) / ((1 - P_{wv}) * (1 - theta)))}
#'
#'   \deqn{F_{p} = Pressure * (((1 - P_{wv} / Pressure) * (1 - theta *
#'   Pressure)) / ((1 - P_{wv}) * (1 - theta)))}
#'
#'   \eqn{P_{wv}} is water vapour pressure, which depends on temperature and
#'   salinity and \eqn{theta} depends on the second virial coefficient of
#'   oxygen, and is calculated using temperature (see Table 2 of Benson and
#'   Krause 1984).
#'
#'   \bold{References}
#'
#'   Benson, Bruce B., Krause, Daniel, (1984), \emph{The concentration and
#'   isotopic fractionation of oxygen dissolved in freshwater and seawater in
#'   equilibrium with the atmosphere, Limnology and Oceanography}, 3, doi:
#'   10.4319/lo.1984.29.3.0620.
#'
#'   USGS. \emph{Change to Solubility Equations for Oxygen in Water.} Technical
#'   Memorandum 2011.03. USGS Office of Water Quality, 2011.
#'
#' @param dat_wide Data frame with at least one column:
#'   \code{temperature_degree_c}. Corresponding pressure (atm) data may be
#'   included in column \code{pressure_atm}. Additional columns will be ignored
#'   and returned.
#'
#' @param p_atm A single value of barometric pressure (atm). This value should
#'   be specified if there is no \code{pressure_atm} column in \code{dat_wide}.
#'   Default is \code{p_atm = NULL}. Note: if \code{p_atm} is specified when
#'   there is a \code{pressure_atm} column in \code{dat_wide}, function will
#'   stop with an error.
#'
#' @param sal A single value of salinity (psu). This value must be specified if
#'   there is no \code{salinity_psu} column in \code{dat_wide}. Default is
#'   \code{sal = NULL}. Note: if \code{sal} is specified when there is a
#'   \code{salinity_psu} column in \code{dat_wide}, the function will stop with
#'   an error.
#'
#' @return Returns \code{dat_wide} with additional column(s), \code{F_p} and
#'   \code{Pressure}.
#'
#' @importFrom dplyr %>% mutate
#'
#' @export


do_pressure_correction <- function(dat_wide, sal = NULL, p_atm = NULL){

  # Error Messages ----------------------------------------------------------

  cols <- colnames(dat_wide)

  if("pressure_atm" %in% cols & !is.null(p_atm)){

    stop("Conflicting pressure values.
         \nHINT: Remove column Pressure or set p_atm argument to NULL.")

  }

  if("salinity_psu" %in% cols & !is.null(sal)){

    stop("Conflicting salinity values.
         \nHINT: Remove column Salinity or set sal argument to NULL.")

  }

# Calculate F_p -----------------------------------------------------------

  if(!is.null(p_atm)) dat_wide <- mutate(dat_wide, pressure_atm = p_atm)
  if(!is.null(sal)) dat_wide <- mutate(dat_wide, salinity_psu = sal)

  dat_wide %>%
    mutate(
      temperature = as.numeric(temperature_degree_c),

      # temperature in Kelvin
      T_Kelvin = temperature_degree_c + 273.15,

      # Coefficient that depends on the Second Virial Coefficient of oxygen
      theta = 0.000975 - (1.426e-5) * temperature_degree_c  +
        (6.436e-8) * temperature_degree_c^2,

      # partial pressure of water vapour in atm
      P_wv = (1 - 5.370e-4 * salinity_psu) *
        exp(
          18.1973 * (1 - 373.16 / T_Kelvin) +
            3.1813e-7 * (1 - exp(26.1205 * (1 - T_Kelvin / 373.16))) -
            1.8726e-2 * (1 - exp(8.03945 * (1 - 373.16 / T_Kelvin))) +
            5.02802 * log(373.16 / T_Kelvin)
        ),

      # square brackets term in Equation 24
      alt_correction = (((1 - P_wv /pressure_atm) * (1 - theta * pressure_atm)) /
        ((1 - P_wv) * (1 - theta))),

      F_p = pressure_atm * alt_correction,
      F_p = round(F_p, digits = 4)

    )

}















