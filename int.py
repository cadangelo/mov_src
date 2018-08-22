import scipy.integrate as integrate
import scipy.special as special
from scipy.integrate import quad, romberg
import numpy as np
from numpy import pi
import matplotlib.pyplot as plt

def flux(t, d, l, v, s):
    return  s/4/np.pi/np.exp(-l*t)/(d**2 + v**2 * t**2)

def calc_fluence(tot_t, d, l, v, s):
    return quad(flux, 0, tot_t, args = (d, l, v, s))

def calc_dose(fluence, E, mu_rho):
    return float(fluence)*E*mu_rho


def main():
    filename = open("test.txt", 'w+')
    d_0 = 100.0 #original distance from source to detector
    s_0 = 1.0 #original source strength
    #tot_d = 500.0 #total distance from orig to final src position
    tot_d = 1000.0 #total distance from orig to final src position 
    E = 1.6e-13 #source energy is 1 MeV = 1.6e-13 J
    mu_rho = 31.08 #mass absorption coeff for tissue
    dose_v_t = [] #list of total dose as a function of total time
    time = [] #list of total times from orig to final pos
    velocity = [1, 5]
    velocity.extend(range(10, 210, 10)) #range of velocities of moving source
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
            flu = calc_fluence(tot_t, d_0, lam, v, s_0)
            d = calc_dose(flu[0], E, mu_rho)
            #d = calc_dose(flu, E, mu_rho)
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
    plt.show()


if __name__ == '__main__':
    main()
