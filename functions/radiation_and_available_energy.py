import os
from pyTSEB import TSEB
from pyTSEB import energy_combination_ET as pet
from pyTSEB import meteo_utils as met
from pyTSEB import net_radiation as rad
from pyTSEB import clumping_index as ci
from pypro4sail import machine_learning_regression as inv
from pypro4sail import prospect as pro
from pypro4sail import four_sail as sail
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import ipywidgets as w
from glob import glob

print("Thanks! libraries imported")
print("You can now continue with the following tasks")

slide_kwargs = {"continuous_update": False}
FIGSIZE = (12.0, 8.0)

np.seterr(all="ignore")
INPUT_FOLDER = os.path.join(os.path.dirname(os.path.dirname(__file__)), "input")
OUTPUT_FOLDER = os.path.join(os.path.dirname(os.path.dirname(__file__)), "output")
SOIL_FOLDER = os.path.join(INPUT_FOLDER, "soil_spectral_library")
ALBEDO_WEIGHT_FILE = os.path.join(INPUT_FOLDER, "albedo_weights.csv")
N_SIM = 100
pet.ITERATIONS = 5
# Generate the list with VZAs (from 0 to 89)
FLUX_LIMS = 0, 500
ET_LIMS = 0, 10
LAI = 1.5
MIN_LAI = 0
MAX_LAI = 4
H_C = 0.12  # Canopy height
LEAF_WIDTH = 0.1  # Leaf width (m)
# Create list of heights
EMIS_C = 0.98
EMIS_S = 0.05
# Calculations are carried out for the following meteorological conditions: R, =
# 300 W m-2; D = 0,10,20 mb; T, = 25 “C; and u = 2ms-’. Such meteorological conditions
# might be considered typical for midday in the middle of a growing season at a subtropical
# site. However, the objective is not to make detailed predictions for particular meteoro-
# logical conditions, it is rather to illustrate the general features of the theoretical treatment
# described.
SDN = 300.  # We assume sn=300 to try to approach rn=400 when computing ln in the module
TAIR = 25 + 273.15  # K
PRESS = 1013.  # mb
TIME_START = 0
TIME_END = 24
LATITUDE = 40
HOURS = np.linspace(TIME_START, TIME_END, N_SIM)
ALBEDO_WEIGHTS = pd.read_csv(ALBEDO_WEIGHT_FILE)


w_kappa = w.FloatSlider(min=0, max=1, value=0.5, step=0.01,
                        description='$\kappa$ (--)',
                        description_tooltip="Extinction coefficient",
                        **slide_kwargs)

w_cab = w.FloatSlider(value=inv.MEAN_CAB,
                      min=inv.MIN_CAB,
                      max=inv.MAX_CAB,
                      step=1, description='Cab ($\mu$g/cm²)',
                      description_tooltip="Chlorophyll content",
                      **slide_kwargs)

w_cw = w.FloatSlider(value=inv.MEAN_CW,
                     min=inv.MIN_CW,
                     max=inv.MAX_CW,
                     step=0.001, description='Cw (g/cm²)',
                     description_tooltip="Leaf water content",
                     readout_format=".3f",
                     **slide_kwargs)

w_lai = w.FloatSlider(value=inv.MEAN_LAI,
                      min=MIN_LAI,
                      max=MAX_LAI,
                      step=0.1, description='LAI (m²/m²)',
                      description_tooltip="Landscape Leaf Area Index",
                      **slide_kwargs)


w_leaf_angle = w.FloatSlider(value=57, min=1, max=90, step=1,
                             description='LIDF (deg.)',
                             description_tooltip="Dominmant leaf zenith angle",
                             **slide_kwargs)

w_sza = w.FloatSlider(value=37, min=0, max=89, step=1,
                      description='SZA (deg.)',
                      description_tooltip="Solar zenith angle",
                      **slide_kwargs)

w_saa = w.FloatSlider(value=180, min=0, max=359, step=1,
                      description='SAA (deg.)',
                      description_tooltip="Solar azimuth angle",
                      **slide_kwargs)

SOIL_FILES = sorted(glob(os.path.join(SOIL_FOLDER, "*.txt")))
SOIL_TYPES = [os.path.splitext(os.path.basename(i))[0] for i in SOIL_FILES]

w_soil = w.Dropdown(options=SOIL_TYPES,
                    value=SOIL_TYPES[0],
                    description='Soil type',
                    description_tooltip='Select type of soil',)

w_hc = w.FloatSlider(min=0.01, max=8, value=4, step=0.01, description='$h_c$ (m)',
                     description_tooltip="Canopy height",
                     **slide_kwargs)

w_hb_ratio = w.FloatSlider(min=0.0, max=0.9, value=0.5, step=0.01,
                           description="$h_{bottom}$ (--)",
                           description_tooltip='Bottom of the canopy relative height',
                           **slide_kwargs)

w_wc = w.FloatSlider(min=0.1, max=3, value=1, step=0.01, description='Forma',
                     description_tooltip="Canopy shape as the ratio between canopy width and canopy height",
                     **slide_kwargs)

w_fc = w.FloatSlider(min=0.05, max=1, value=0.2, step=0.01,
                     description='$f_c$ (--)',
                     description_tooltip='Canopy fraction',
                     **slide_kwargs)

w_interrow = w.FloatSlider(min=0.1, max=6, value=1, step=0.1,
                           description='L (m)',
                           description_tooltip='Distance between rows',
                           **slide_kwargs)

w_psi = w.FloatSlider(min=-90., max=90, value=0, description='Row azimuth (deg.)',
                      description_tooltip='Row azimuth angle',
                      **slide_kwargs)

w_skyl = w.FloatSlider(min=0, max=1, value=0.1, step=0.01,
                       description='Skyl (--)',
                       description_tooltip='Ratio of diffuse radiation',
                       **slide_kwargs)

w_leaf_abs = w.FloatSlider(min=0.8, max=1, value=0.9, step=0.01,
                           description='Abs. hoja',
                           description_tooltip='Leaf absorptance',
                           **slide_kwargs)

w_sdn = w.FloatSlider(min=0., max=400, value=300,
                      description='S$\downarrow$ (W/m²)',
                      description_tooltip='Daily shortwave irradiance',
                      **slide_kwargs)


w_tair = w.FloatSlider(min=273, max=313, value=293, step=0.1,
                       description='T$_{air}$ (K)',
                       description_tooltip="Air temperature",
                       **slide_kwargs)

w_hr = w.FloatSlider(min=0, max=100, value=50, step=0.1,
                     description='HR (%)',
                     description_tooltip="Relative humidity",
                     **slide_kwargs)

w_emiss = w.FloatSlider(min=0.9, max=0.99, value=0.98, step=0.01,
                        description='$\epsilon$',
                        description_tooltip="Surface emissivity",
                        **slide_kwargs)

w_lat = w.FloatSlider(min=0.0, max=90, value=40, step=1,
                      description='Lat. (deg.)',
                      description_tooltip="Site latitude",
                      **slide_kwargs)


def beer_lambert_law(kappa, length):
     tau = np.exp(-kappa * length)
     plot_attenuation(length, tau)

def plot_attenuation(length, tau):
    fig, axs = plt.subplots(figsize=FIGSIZE)
    # Create an empty list with the returning f_c for each canopy
    axs.plot(length, tau, 'k-')
    axs.set_xlabel('Path length (arbitrary units)')
    axs.set_ylabel('Trasmittance')
    axs.set_ylim((0, 1))
    axs.set_xlim((np.min(length), np.max(length)))
    plt.tight_layout()
    plt.show()


def plot_fipar(lai, leaf_angle):
    fig, axs = plt.subplots(figsize=FIGSIZE)
    # Create an empty list with the returning f_c for each canopy
    chi = rad.leafangle_2_chi(leaf_angle)
    szas = np.arange(0, 90)
    fipar = TSEB.calc_F_theta_campbell(szas, lai, x_LAD=chi)
    # Plot the results
    axs.plot(szas, fipar, 'r')
    axs.set_xlabel('Solar Zenith Angle (degrees)')
    axs.set_ylabel('Fraction of Interecptad Radiation')
    axs.set_ylim((0, 1.05))
    axs.set_xlim((0, 90))
    plt.tight_layout()
    plt.show()


def solar_irradiance(sdn, szas):
    rdirvis, rdifvis, rdirnir, rdifnir = rad.calc_potential_irradiance_weiss(szas,
                                                                             press=np.full(N_SIM, 1013.15))
    sdn_pot = rdirvis + rdifvis + rdirnir + rdifnir
    sdn = sdn * np.size(sdn_pot) * sdn_pot / np.sum(sdn_pot)
    return sdn


def relative_azimuth(angle1, angle2, degrees=True):
    if degrees:
        ref = 180.
        divisor = 360.
    else:
        ref = np.pi
        divisor = 2 * np.pi
    psi = ref - np.abs(ref - np.mod(np.abs(angle1 - angle2), divisor))
    return psi


def plot_net_radiation_rows(lai, leaf_angle, h_c, f_c, row_distance, row_direction,
                            sdn_day=1, skyl=0.1, fvis=0.55, lat=LATITUDE,
                            cab=inv.MEAN_CAB, cw=inv.MEAN_CW, soil_type=SOIL_TYPES[0]):

    rho_leaf, tau_leaf, rho_soil = get_surface_spectra(cab, cw, soil_type)

    chi = rad.leafangle_2_chi(leaf_angle)
    local_lai = np.full(N_SIM, lai) / f_c
    szas, saas = met.calc_sun_angles(np.full(N_SIM, lat),
                                     np.zeros(N_SIM),
                                     np.zeros(N_SIM),
                                     np.full(N_SIM, 180),
                                     HOURS)
    sundown = szas > 89
    sdn_dir = (1. - skyl)
    sdn_dif = skyl
    w_c = row_distance * f_c / np.full(N_SIM, h_c)
    psi = relative_azimuth(row_direction, saas)
    omega = ci.calc_omega_rows(np.full(N_SIM, lai),
                               np.full(N_SIM, f_c),
                               theta=szas,
                               psi=psi,
                               w_c=w_c,
                               x_lad=np.full(N_SIM, chi))

    lai_eff = local_lai * omega
    # Compute radiation for row crop canopy
    sn_v, sn_s = rad.calc_Sn_Campbell(lai,
                                      szas,
                                      sdn_dir,
                                      sdn_dif,
                                      fvis,
                                      1. - fvis,
                                      np.full(N_SIM, rho_leaf[0]),
                                      np.full(N_SIM, tau_leaf[0]),
                                      np.full(N_SIM, rho_leaf[1]),
                                      np.full(N_SIM, tau_leaf[1]),
                                      np.full(N_SIM, rho_soil[0]),
                                      np.full(N_SIM, rho_soil[1]),
                                      x_LAD=chi,
                                      LAI_eff=lai_eff)

    sn_v[sundown] = np.nan
    sn_s[sundown] = np.nan
    # Compute radiation for homogeneous canopy
    sn_v_0, sn_s_0 = rad.calc_Sn_Campbell(lai,
                                          szas,
                                          sdn_dir,
                                          sdn_dif,
                                          fvis,
                                          1. - fvis,
                                          np.full(N_SIM, rho_leaf[0]),
                                          np.full(N_SIM, tau_leaf[0]),
                                          np.full(N_SIM, rho_leaf[1]),
                                          np.full(N_SIM, tau_leaf[1]),
                                          np.full(N_SIM, rho_soil[0]),
                                          np.full(N_SIM, rho_soil[1]),
                                          x_LAD=chi,
                                          LAI_eff=None)

    sn_v_0[sundown] = np.nan
    sn_s_0[sundown] = np.nan
    sn = sn_v + sn_s
    sn_0 = sn_v_0 + sn_s_0
    fig, axs = plt.subplots(figsize=FIGSIZE)
    # Plot the results
    axs.plot(HOURS, sn, 'black', label="Row crop, total")
    axs.plot(HOURS, sn_0, 'black', ls="--", label="Homogenous crop, total")
    axs.plot(HOURS, sn_v, 'green', label="Row crop, canopy")
    axs.plot(HOURS, sn_v_0, 'green', ls="--",label="Homogenous crop, canopy")
    axs.set_xlabel('Time (hours)')
    if sdn_day != 1:
        axs.set_ylabel('Absorbed radiation (--)')
    else:
        axs.set_ylabel('Fraction of absorbed radiation (--)')
    axs.set_xlim((TIME_START, TIME_END))
    axs.set_ylim((0, 1))
    axs.legend()
    plt.show()


def plot_net_solar_radiation(lai, leaf_angle, h_c, f_c, row_distance, row_direction,
                             sdn_day, skyl=0.1, fvis=0.55, lat=LATITUDE,
                             cab=inv.MEAN_CAB, cw=inv.MEAN_CW, soil_type=SOIL_TYPES[0]):

    rho_leaf, tau_leaf, rho_soil = get_surface_spectra(cab, cw, soil_type=soil_type)
    chi = rad.leafangle_2_chi(leaf_angle)
    local_lai = np.full(N_SIM, lai) / f_c
    szas, saas = met.calc_sun_angles(np.full(N_SIM, lat),
                                     np.zeros(N_SIM),
                                     np.zeros(N_SIM),
                                     np.full(N_SIM, 180),
                                     HOURS)

    sdn = solar_irradiance(sdn_day, szas)
    sdn_dir = (1. - skyl) * sdn
    sdn_dif = skyl * sdn
    if f_c == 1:
        omega = 1.0
    else:
        w_c = row_distance * f_c / np.full(N_SIM, h_c)
        psi = relative_azimuth(row_direction, saas)
        omega = ci.calc_omega_rows(np.full(N_SIM, lai),
                                   np.full(N_SIM, f_c),
                                   theta=szas,
                                   psi=psi,
                                   w_c=w_c,
                                   x_lad=np.full(N_SIM, chi))

    lai_eff = local_lai * omega
    sn_v, sn_s = rad.calc_Sn_Campbell(lai,
                                      szas,
                                      sdn_dir,
                                      sdn_dif,
                                      fvis,
                                      1. - fvis,
                                      np.full(N_SIM, rho_leaf[0]),
                                      np.full(N_SIM, tau_leaf[0]),
                                      np.full(N_SIM, rho_leaf[1]),
                                      np.full(N_SIM, tau_leaf[1]),
                                      np.full(N_SIM, rho_soil[0]),
                                      np.full(N_SIM, rho_soil[1]),
                                      x_LAD=chi,
                                      LAI_eff=lai_eff)

    sundown = szas > 89
    sn_v[sundown] = np.nan
    sn_s[sundown] = np.nan
    sn = sn_v + sn_s
    albedo = 1. - sn / sdn
    fig, axs = plt.subplots(figsize=FIGSIZE)
    # Plot the results
    axs.plot(HOURS, sn, 'black', label="Net Radiation")
    axs.plot(HOURS, sn_v, 'green', label="Net Canopy Radiation")
    axs.plot(HOURS, sn_s, 'orange', label="Net Soil Radiation")
    axs.set_xlabel('Time')
    axs.set_ylabel('Solar Radiation (W/m²)')
    axs.set_xlim((TIME_START, TIME_END))
    axs.set_ylim((0, 1000))
    axs.legend()
    secax = axs.twinx()
    secax.plot(HOURS, albedo, "b")
    secax.set_ylabel('Albedo', color="blue")
    secax.tick_params(axis='y', colors='blue')
    secax.set_ylim((0, 1))

    if f_c < 1:
        # Compute radiation for homogeneous canopy
        sn_v_0, sn_s_0 = rad.calc_Sn_Campbell(lai,
                                              szas,
                                              1 - skyl,
                                              skyl,
                                              fvis,
                                              1. - fvis,
                                              np.full(N_SIM, rho_leaf[0]),
                                              np.full(N_SIM, tau_leaf[0]),
                                              np.full(N_SIM, rho_leaf[1]),
                                              np.full(N_SIM, tau_leaf[1]),
                                              np.full(N_SIM, rho_soil[0]),
                                              np.full(N_SIM, rho_soil[1]),
                                              x_LAD=chi,
                                              LAI_eff=None)

        sn_v_0[sundown] = np.nan
        sn_s_0[sundown] = np.nan
        sn_0 = sn_v_0 + sn_s_0
        fig, axs = plt.subplots(figsize=FIGSIZE)
        # Plot the results
        axs.plot(HOURS, sn / sdn, 'black', label="Row crop, total")
        axs.plot(HOURS, sn_0, 'black', ls="--", label="Homogenous crop, total")
        axs.plot(HOURS, sn_v / sdn, 'green', label="Row crop, canopy")
        axs.plot(HOURS, sn_v_0, 'green', ls="--", label="Homogenous crop, canopy")
        axs.set_xlabel('Time (hours)')
        axs.set_ylabel('Fraction of absorbed radiation (--)')
        axs.set_xlim((TIME_START, TIME_END))
        axs.set_ylim((0, 1))
        axs.legend()
        plt.show()

    plt.show()


def plot_apar(local_lai, leaf_angle, h_c, f_c, row_distance,
              row_direction, sdn_day, skyl=0.1, fvis=0.55,
              lat=LATITUDE, cab=inv.MEAN_CAB, cw=inv.MEAN_CW, soil_type=SOIL_TYPES[0]):

    rho_leaf, tau_leaf, rho_soil = get_surface_spectra(cab, cw, soil_type=soil_type)

    szas, saas = met.calc_sun_angles(np.full(N_SIM, lat),
                                     np.zeros(N_SIM),
                                     np.zeros(N_SIM),
                                     np.full(N_SIM, 180),
                                     HOURS)

    chi = rad.leafangle_2_chi(leaf_angle)
    lai = np.full(N_SIM, local_lai) * f_c
    sdn = solar_irradiance(sdn_day, szas)
    sdn_dir = (1. - skyl) * sdn
    sdn_dif = skyl * sdn
    w_c = row_distance * f_c / np.full(N_SIM, h_c)
    omega = ci.calc_omega_rows(np.full(N_SIM, lai),
                               np.full(N_SIM, f_c),
                               theta=szas,
                               psi=row_direction - saas,
                               w_c=w_c,
                               x_lad=np.full(N_SIM, chi))
    lai_eff = local_lai * omega
    albb, albd, taubt, taudt = rad.calc_spectra_Cambpell(lai,
                                                         szas,
                                                         rho_leaf[0],
                                                         tau_leaf[0],
                                                         np.full(N_SIM, rho_soil[0]),
                                                         x_lad=chi,
                                                         lai_eff=lai_eff)
    apar = (1.0 - taubt) * (1.0 - albb) * sdn_dir * fvis + \
           (1.0 - taudt) * (1.0 - albd) * sdn_dif * fvis

    akb = rad.calc_K_be_Campbell(szas, chi, radians=False)
    taub = np.exp(-akb * lai_eff)
    taud = rad._calc_taud(chi, lai)
    ipar = (1.0 - taub) * sdn_dir * fvis + (1.0 - taud) * sdn_dif * fvis
    par = sdn_dir * fvis + sdn_dif * fvis
    # Plot the results
    fig, axs = plt.subplots(figsize=FIGSIZE)
    axs.plot(HOURS, par, 'black', label="Entrante")
    axs.plot(HOURS, apar, 'green', label="Absorbida")
    axs.plot(HOURS, ipar, 'red', label="Interceptada")
    axs.set_xlabel('Time')
    axs.set_ylabel('PAR (W/m²)')
    axs.set_xlim((TIME_START, TIME_END))
    axs.set_ylim((0, 1000))
    axs.legend()
    plt.show()


def relative_azimuth(angle1, angle2, degrees=True):
    if degrees:
        ref = 180.
        divisor = 360.
    else:
        ref = np.pi
        divisor = 2 * np.pi
    psi = ref - np.abs(ref - np.mod(np.abs(angle1 - angle2), divisor))
    return psi


def get_surface_spectra(cab, cw=inv.MEAN_CW, soil_type=SOIL_TYPES[0]):
    rho_leaf = []
    tau_leaf = []
    rho_soil = []
    wls, r, t = pro.prospectd(inv.MEAN_N_LEAF, cab, inv.MEAN_CAR,
                              inv.MEAN_CBROWN, cw, inv.MEAN_CM,
                              inv.MEAN_ANT)

    rsoil = np.genfromtxt(os.path.join(SOIL_FOLDER, f"{soil_type}.txt"))
    rho_leaf.append(np.sum(r * ALBEDO_WEIGHTS["w_PAR"]))
    tau_leaf.append(np.sum(t * ALBEDO_WEIGHTS["w_PAR"]))
    rho_soil.append(np.sum(rsoil[:, 1] * ALBEDO_WEIGHTS["w_PAR"]))
    rho_leaf.append(np.sum(r * ALBEDO_WEIGHTS["w_NIR"]))
    tau_leaf.append(np.sum(t * ALBEDO_WEIGHTS["w_NIR"]))
    rho_soil.append(np.sum(rsoil[:, 1] * ALBEDO_WEIGHTS["w_NIR"]))
    rho_leaf.append(np.sum(r * ALBEDO_WEIGHTS["w_SW"]))
    tau_leaf.append(np.sum(t * ALBEDO_WEIGHTS["w_SW"]))
    rho_soil.append(np.sum(rsoil[:, 1] * ALBEDO_WEIGHTS["w_SW"]))

    return rho_leaf, tau_leaf, rho_soil


def bidirectional_reflectance(cab, cw, lai, leaf_angle, sza=35., saa=180., skyl=0.1,
                              soil_type=SOIL_TYPES[0]):
    wls_sim = ["PAR", "SW"]
    vzas = np.linspace(0, 90, num=50)
    vaas = np.linspace(0, 360, num=60)
    step_vza = np.radians(vzas[1] - vzas[0])
    step_psi = np.radians(vaas[1] - vaas[0])
    vzas, vaas = np.meshgrid(vzas, vaas, indexing="ij")
    n_sims = np.size(vzas)
    cosvza = np.cos(np.radians(vzas.reshape(-1)))
    sinvza = np.sin(np.radians(vzas.reshape(-1)))
    lidf = sail.calc_lidf_campbell_vec(np.full(n_sims, leaf_angle))
    rho_leaf, tau_leaf, rho_soil = get_surface_spectra(cab, cw, soil_type=soil_type)
    rho_leaf2 = np.tile([rho_leaf[0], rho_leaf[2]], (n_sims, 1))
    tau_leaf2 = np.tile([tau_leaf[0], tau_leaf[2]], (n_sims, 1))
    rho_soil2 = np.tile([rho_soil[0], rho_soil[2]], (n_sims, 1))
    psi = relative_azimuth(vaas.reshape(-1), saa)

    [_,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     _,
     rdot,
     _,
     _,
     rsot,
     _,
     _,
     _] = sail.foursail_vec(np.full(n_sims, lai),
                            np.full(n_sims, inv.MEAN_HOTSPOT),
                            lidf,
                            np.full(n_sims, sza),
                            vzas.reshape(-1),
                            psi,
                            rho_leaf2.T,
                            tau_leaf2.T,
                            rho_soil2.T)


    r2 = rdot * skyl + rsot * (1 - skyl)
    if lai == 0:
        r2 = rho_soil.T
    albedo = np.sum(r2 * cosvza * sinvza * step_vza * step_psi / np.pi, axis=1)
    polar_plot(vzas, vaas, r2, sza, saa, wls_sim, albedo)

    return wls_sim, r2


def polar_plot(vzas, vaas, rhos, sza, saa, wls_sim, albedo):
    fig, ax = plt.subplots(ncols=len(wls_sim), figsize=FIGSIZE,
                           subplot_kw={'projection': 'polar'})
    for i, wl in enumerate(wls_sim):
        ax[i].set_theta_direction(-1)
        ax[i].set_theta_offset(np.pi / 2.0)
        im = ax[i].pcolormesh(np.radians(vaas), vzas, rhos[i].reshape(vzas.shape),
                               vmin=0, vmax=0.25, cmap="inferno", shading='auto')
        ax[i].plot(np.radians(saa), sza, marker="*", markersize=12,
                   markerfacecolor="none", markeredgecolor="black")
        ax[i].grid(True)
        ax[i].set_title(f"{wls_sim[i]}\nalbedo={albedo[i]:.2f}")

    plt.tight_layout()
    cbar_ax = fig.add_axes([0.1, 0.1, 0.8, 0.03])
    fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
    plt.show()


def plot_longwave_radiation(t_air, hr, delta_t, emiss=0.98):
    ea = hr * met.calc_vapor_pressure(t_air)
    emiss_atm = rad.calc_emiss_atm(ea, t_air)
    l_sky = emiss_atm * met.calc_stephan_boltzmann(t_air)
    lst = t_air + delta_t
    ln = emiss * (l_sky - met.calc_stephan_boltzmann(lst))
    plt.figure(figsize=FIGSIZE)
    plt.plot(lst, ln, "k-")
    plt.xlabel("LST (K)")
    plt.xlim((np.min(lst), np.max(lst)))
    plt.ylabel("Net longwave radiation (W/m²)")
    plt.ylim((-200, 400))
    plt.tight_layout()
    plt.show()
    return ln




