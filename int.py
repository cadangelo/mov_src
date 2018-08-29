import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad, romberg, romb, trapz
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from pylab import imshow, show, get_cmap
from math import atan
from decimal import Decimal

d_0 = 100.0 #original distance from source to detector
s_0 = 100 #original source strength
#tot_d = 500.0 #total distance from orig to final src position
tot_d = 10000.0 #total distance from orig to final src position 
E = 1.6e-13 #source energy is 1 MeV = 1.6e-13 J
mu_rho = 31.08 #mass absorption coeff for tissue
velocity = [1, 5]
velocity.extend(range(10, 210, 10)) #range of velocities of moving source


def flux(t, d, l, v, s):
    return  s/4/np.pi*np.exp(-l*t)/(d**2 + v**2 * t**2)

def get_num_samples(order_num_samples, int_type):
    if int_type == 'fixed_romberg':
        return 2**order_num_samples + 1
    elif int_type == 'trapz':
        return order_num_samples*10
 
def hand_calc_fluence(s, v, d, t):
    return s/4/np.pi/v/d*atan(v*t/d)

def calc_fluence(tot_t, d_0, l, v, s_0, order_num_samples, int_type):
    samples = []
    # romberg integration of flux over time
    if int_type == 'romberg':
        return romberg(flux, 0, tot_t, args = (d_0, l, v, s_0))
    # hand calculated fluence
    elif int_type == 'hand_calc':
        return hand_calc_fluence(s_0, v, d_0, tot_t)
    elif int_type == 'trapz':
        num_samples = get_num_samples(order_num_samples, int_type)
        t_pts = np.linspace(0, tot_t, num_samples, endpoint=True)
        dt = t_pts.item(1) - t_pts.item(0)
        for t in t_pts:
            sample = flux(t, d_0, l, v, s_0)
            samples.append(sample)
        return trapz(samples, x=t_pts, dx=dt)
    # fixed sample romberg integration of flux over time
    elif int_type == 'fixed_romberg':
        num_samples = get_num_samples(order_num_samples, int_type)
        t_pts = np.linspace(0, tot_t, num_samples, endpoint=True)
        dt = t_pts.item(1) - t_pts.item(0)
        for t in t_pts:
            sample = flux(t, d_0, l, v, s_0)
            samples.append(sample)
        return romb(samples, dt)

def calc_dose(fluence, E, mu_rho):
    return float(fluence)*E*mu_rho


def plt_err_discrete_from_cont():

    # Varying number of time points for fixed-sample numerical integration 
    # No source decay
    order_num_samples = np.linspace(1, 8, 8, endpoint=True)
    error = []
    num_samples_list = []
    time = []
    error_matrix = np.empty([len(order_num_samples), len(velocity)])
    cont_int_type = 'hand_calc'# or 'romberg'
    fixed_int_type = 'fixed_romberg' #or trapz
    for i, order in enumerate(order_num_samples):
        del time[:]
        num_samples_list.append(get_num_samples(order, fixed_int_type))
        for j, v in enumerate(velocity):
            tot_t = tot_d/v #total time for src to move from orig to final pos
            time.append('%.2f'%tot_t)
            cont_flu = calc_fluence(tot_t, d_0, 0, v, s_0, 1, cont_int_type)
            cont_d = calc_dose(cont_flu, E, mu_rho)
            fixed_flu = calc_fluence(tot_t, d_0, 0, v, s_0, order, fixed_int_type)
            fixed_d = calc_dose(fixed_flu, E, mu_rho)
            error_matrix[i, j] = '%.1E'%Decimal(100*(fixed_d-cont_d)/cont_d)

    # Plot error map
    fig, ax = plt.subplots()
    im = ax.imshow(error_matrix, cmap=get_cmap("coolwarm"), origin='lower',interpolation='nearest', vmin=-10, vmax=10)

    ax.set_xticks(np.arange(len(time)))
    ax.set_yticks(np.arange(len(num_samples_list)))
    ax.set_xticklabels(time)
    ax.set_yticklabels(num_samples_list)
    for i in range(len(num_samples_list)):
        for j in range(len(time)):
            text = ax.text(j, i, error_matrix[i,j], ha="center", va="center", color="black")

    ax.set_title("Error according to number of points chosen for integration")
    ax.set_xlabel("Total time of geometry movement [s]")
    ax.set_ylabel("Number of time points")
    fig.tight_layout()
    cbar = fig.colorbar(im, ticks=[-10,0,10])
    cbar.ax.set_yticklabels(['-10.0', '0', '10.0'])
    cbar.set_label('% Error from analytical calculation')
    plt.show()

def plt_err_ignoring_decay():
    filename = open("test.txt", 'w+')
    dose_v_t = [] #list of total dose as a function of total time
    time = [] #list of total times from orig to final pos
    # Varying half-lives: 1 min, 10 min, 30 min, 1 hr, 8 hr, 24 hr, 48 hr, 72 hr, 1 wk, 2 wk, 3 wk, 4 wk
#    half_life = [60, 600, 1800, 3600, 28800, 86400, 172800, 259200, 604800, 1210000, 1814000, 241900, 999]
    labels = ['1 min', '5 min', '1 hr', '8 hr', '1 wk', '1 yr', 'no decay'] 
    half_life = [60, 300, 3600, 28800, 604800, 31540000, 999]
    error_matrix = np.empty([len(half_life), len(velocity)])
    for i, (hl, label) in enumerate(zip(half_life, labels)):
        del dose_v_t[:]
        del time[:]
        if hl == 999:
            lam = 0
        else:
            lam = 0.693/hl #decay constant
        for j, v in enumerate(velocity):
            tot_t = tot_d/v #total time for src to move from orig to final pos
            time.append('%.2f'%tot_t)
            flu = calc_fluence(tot_t, d_0, lam, v, s_0, 1, 'romberg')
            d = calc_dose(flu, E, mu_rho)
            no_dec_flu = calc_fluence(tot_t, d_0, 0, v, s_0, 1, 'romberg')
            no_dec_d = calc_dose(no_dec_flu, E, mu_rho)
            dose_v_t.append(d)
            error_matrix[i, j] = '%.1E'%Decimal(100*(no_dec_d-d)/d)

        plt.plot(time, dose_v_t, label=label)
        print >> filename, str(label)
        for t, dose in zip(time, dose_v_t):
            print >> filename, "{0:<15s} {1}".format(str(t), str(dose))
   
    # Plot dose w/ source decay at various HL
    plt.title("Total Dose with Moving Source Decaying at Varying Half-Lives")
    plt.xlabel("Total time of Source Movement [s]") 
    plt.ylabel("Total Dose [Gy]")
    plt.yscale('log')
    plt.legend()
    #plt.show()

    # Plot error map
    fig, ax = plt.subplots()
    im = ax.imshow(error_matrix, cmap=get_cmap("coolwarm"), origin='lower',interpolation='nearest', vmin=-10, vmax=10)

    ax.set_xticks(np.arange(len(time)))
    ax.set_yticks(np.arange(len(half_life)))
    ax.set_xticklabels(time)
    ax.set_yticklabels(labels)
    for i in range(len(half_life)):
        for j in range(len(time)):
            text = ax.text(j, i, error_matrix[i,j], ha="center", va="center", color="black")
    ax.set_title("Error by ignoring decay")
    ax.set_xlabel("Total time of geometry movement [s]")
    ax.set_ylabel("Half-life of source [s]")
    fig.tight_layout()
    cbar = fig.colorbar(im, ticks=[-10,0,10])
    cbar.ax.set_yticklabels(['-10', '0', '10'])
    cbar.set_label('% Error from ignoring decay')

    plt.show()

def main():

    plt_err_ignoring_decay()
   
    plt_err_discrete_from_cont()

        


if __name__ == '__main__':
    main()
