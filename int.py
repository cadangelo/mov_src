import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad, romberg, romb
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt
from pylab import imshow, show, get_cmap
from math import atan

d_0 = 100.0 #original distance from source to detector
s_0 = 1.0 #original source strength
#tot_d = 500.0 #total distance from orig to final src position
tot_d = 1000.0 #total distance from orig to final src position 
E = 1.6e-13 #source energy is 1 MeV = 1.6e-13 J
mu_rho = 31.08 #mass absorption coeff for tissue
velocity = [1, 5]
velocity.extend(range(10, 210, 10)) #range of velocities of moving source


def flux(t, d, l, v, s):
    return  s/4/np.pi*np.exp(-l*t)/(d**2 + v**2 * t**2)

def get_num_samples(order_num_samples, int_type):
    if int_type == 'romberg':
        return 2**order_num_samples + 1
    #elif int_type == 'trapz':
    #    return order_num_samples*10
 
def hand_calc_fluence(s, v, d, t):
    return s/4/np.pi/v/d*atan(v*t/d)

def calc_fluence(tot_t, d_0, l, v, s_0, order_num_samples):
    samples = []
    # romberg integration of flux over time
    if order_num_samples == -1:
        return romberg(flux, 0, tot_t, args = (d_0, l, v, s_0))
    # hand calculated fluence
    elif order_num_samples == 0:
        return hand_calc_fluence(s_0, v, d_0, tot_t)
   # elif order_num_samples == -2:
   #     num_samples = get_num_samples(order_num_samples, 'trapz')
   #     t_pts = np.linspace(0, tot_t, num_samples, endpoint=True)
   #     dt = t_pts.item(1) - t_pts.item(0)
   #     return trapz(samples, x=t_pts, dx=dt)
    # fixed sample romberg integration of flux over time
    else:
    #    del samples[:]
        num_samples = get_num_samples(order_num_samples, 'romberg')
        #print "tot_t, order, num_samples", tot_t, order_num_samples, num_samples
        t_pts = np.linspace(0, tot_t, num_samples, endpoint=True)
        dt = t_pts.item(1) - t_pts.item(0)
        for t in t_pts:
            sample = flux(t, d_0, l, v, s_0)
            samples.append(sample)
        return romb(samples, dt)

def calc_dose(fluence, E, mu_rho):
    return float(fluence)*E*mu_rho


def plt_err_discrete_from_cont():

    # Varying number of fixed time points for romberg integration 
    order_num_samples = np.linspace(1, 8, 8, endpoint=True)
    #error = np.empty([len(velocity)])
    error = []
    num_samples_list = []
    time = []
    error_matrix = np.empty([len(order_num_samples), len(velocity)])
    print "size error mat", error_matrix.size
    i = -1
    for order in order_num_samples:
        del time[:]
        del error[:]
        num_samples_list.append(2**order+1)
        i = i+1
        j = -1
        for v in velocity:
            tot_t = tot_d/v #total time for src to move from orig to final pos
            time.append('%.2f'%tot_t)
            flu = calc_fluence(tot_t, d_0, 0, v, s_0, order)
            d = calc_dose(flu, E, mu_rho)
            romb_flu = calc_fluence(tot_t, d_0, 0, v, s_0, -1)
            romb_d = calc_dose(romb_flu, E, mu_rho)
            hand_flu = calc_fluence(tot_t, d_0, 0, v, s_0, 0)
            hand_d = calc_dose(hand_flu, E, mu_rho)
            trapz_flu = calc_fluence(tot_t, d_0, 0, v, s_0, -2)
            trapz_d = calc_dose(hand_flu, E, mu_rho)
            #error.append((d-cont_d)/cont_d)
            j = j+1
            error_matrix[i, j] = '%.2f'%(100*(d-hand_d)/hand_d)
            print "n, t, error" , 2**order+1, tot_t, error_matrix[i,j]
            #np.append(error, (d-cont_d)/cont_d)
            #print "num pts, time, error", 2**order+1, tot_t, error
#        print "error", error
        #print "order, len error", order, np.prod(error.shape)
        #np.append(error_matrix, error, axis=0)

    print error_matrix

    fig, ax = plt.subplots()
    im = ax.imshow(error_matrix, cmap=get_cmap("coolwarm"), interpolation='nearest')

    ax.set_xticks(np.arange(len(time)))
    ax.set_yticks(np.arange(len(num_samples_list)))
    ax.set_xticklabels(time)
    ax.set_yticklabels(num_samples_list)
    for i in range(len(num_samples_list)):
        for j in range(len(time)):
            text = ax.text(j, i, error_matrix[i,j], ha="center", va="center", color="w")

    ax.set_title("Error according to number of points chosen for integration")
    fig.tight_layout()
    cbar = fig.colorbar(im, ticks=[-1,0,1])
    cbar.ax.set_yticklabels(['<-1', '0', '>1'])
    plt.show()

def plt_err_ignoring_decay():
    filename = open("test.txt", 'w+')
    dose_v_t = [] #list of total dose as a function of total time
    time = [] #list of total times from orig to final pos
    # Varying half-lives: 1 min, 10 min, 30 min, 1 hr, 8 hr, 24 hr, 48 hr, 72 hr, 1 wk, 2 wk, 3 wk, 4 wk
#    half_life = [60, 600, 1800, 3600, 28800, 86400, 172800, 259200, 604800, 1210000, 1814000, 241900, 999]
    labels = ['1 min', '5 min', '1 hr', '8 hr', '1 wk', '1 yr', 'no decay'] 
    half_life = [60, 300, 3600, 28800, 86400, 604800, 31540000, 999]
    for hl, label in zip(half_life, labels):
        del dose_v_t[:]
        del time[:]
        if hl == 999:
            lam = 0
        else:
            lam = 0.693/hl #decay constant
        for v in velocity:
            tot_t = tot_d/v #total time for src to move from orig to final pos
            time.append(tot_t)
            order_num_samples = -1
            flu = calc_fluence(tot_t, d_0, lam, v, s_0, order_num_samples)
        #    if order_num_samples < 0:
            print "flu", flu
            d = calc_dose(flu, E, mu_rho)
        #    else:
        #        d = calc_dose(flu, E, mu_rho)
            dose_v_t.append(d)

        plt.plot(time, dose_v_t, label=label)
        print >> filename, str(label)
        for t, dose in zip(time, dose_v_t):
            print >> filename, "{0:<15s} {1}".format(str(t), str(dose))
   
    plt.title("Total Dose with Moving Source Decaying at Varying Half-Lives")
    plt.xlabel("Total time of Source Movement [s]") 
    plt.ylabel("Total Dose [Gy]")
    plt.yscale('log')
    plt.legend()
    #plt.show()


def main():

#    plt_err_ignoring_decay()
   
    plt_err_discrete_from_cont()

        


if __name__ == '__main__':
    main()
