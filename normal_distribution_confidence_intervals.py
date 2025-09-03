import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import t
from scipy.stats import chi2
from scipy.stats import norm


k_values = [5, 20, 50]  # Number of samples
mean     = 20           # Mean value of the normal distribution
sigma    = 2            # Sigma value of the normal distribution

CL = 0.90               # The Confidence Level

random_sets = int(input("Insert the random sets (iterations): "))

b_values = [[] for _ in k_values]  # For merged histograms in Question (b)

for k in k_values:

    # Τhe requested quantities from Exercise 4.
    # For every k, we save the values in these arrays.
    a_values               = []  # For the question (a)
    b_values_seperately    = []  # For the question (b)

    #For the third question of the Exercise 4:
    accepted_mean  = [] # Mean  Values which are into the CB
    accepted_sigma = [] # Sigma Values which are into the CB

    # Percentile of normal distribution for a given CL and a given k:
    t_value = t.ppf(1 - (1 - CL) / 2, k - 1)

    for _ in range(random_sets):

        # Generation of the sample:
        x = np.random.normal(mean, sigma, k)

        # Τhe requested quantities from Exercise 4:
        a = np.sqrt(k) * (np.mean(x) - mean) / np.sqrt(np.var(x))
        b = (k - 1) * np.var(x) / sigma ** 2

        a_values.append(a)                    # For showing the histograms of Question (a)
        b_values[k_values.index(k)].append(b) # For showing the histograms merged.
        b_values_seperately.append(b)         # For showing the histogram seperately.

        # Confidence band calculation for mean:
        low_mean = np.mean(x) - t_value * np.sqrt(np.var(x)) / np.sqrt(k)
        up_mean  = np.mean(x) + t_value * np.sqrt(np.var(x)) / np.sqrt(k)

        # Confidence band calculation for sigma:
        low_sigma = np.sqrt((k - 1) * np.var(x) / chi2.ppf(1 - (1 - CL) / 2, k - 1))
        up_sigma  = np.sqrt((k - 1) * np.var(x) / chi2.ppf(    (1 - CL) / 2, k - 1))

        if random_sets == 1:
            print("===========================================================")
            print(f"For k = {k} samples: CB (Mean) = [{low_mean }, {up_mean }]")
            print(f"For k = {k} samples: CB(Sigma) = [{low_sigma}, {up_sigma}]")
        
        # Min mean and max mean should be into the CB:
        if low_mean <= mean <= up_mean:
            accepted_mean.append(1)  #'1' if it is in the CB
        
        # Min mean and max sigma should be into the CB:
        if low_sigma <= sigma <= up_sigma:
            accepted_sigma.append(1) #'1' if it is in the CB

    #RESULTS FOR QUESTION (a):
    plt.figure(figsize = (8, 6))
    plt.hist(a_values, bins='auto', color='skyblue', edgecolor='black', density=True, label = f"Histogram for k={k}")
    plt.title(r"Parameter $\frac{{\sqrt{{k}} \cdot (\bar{{x}} - \mu)}}{{s}}$ for k = {} and {} random sets".format(k,random_sets), 
              fontsize=16)
    plt.plot(np.linspace(min(a_values), max(a_values), 1000), 
             t.pdf(np.linspace(min(a_values), max(a_values), 1000), k), label = f"Student distribution fit for k = {k}")
    
    plt.xlabel(r'$\frac{\sqrt{k} \cdot (\bar{x} - \mu)}{s}$', fontsize=16)
    plt.ylabel("Density", fontsize=16)
    plt.xlim(-7, 7)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    #RESULTS FOR QUESTION (b):
    plt.figure(figsize = (8, 6))
    plt.hist(b_values_seperately, bins = 'auto', color='red', edgecolor='black', density=True, label = f"Histogram for k={k}")
    plt.plot(np.linspace(min(b_values_seperately), max(b_values_seperately), 1000), 
             chi2.pdf(np.linspace(min(b_values_seperately), max(b_values_seperately), 1000), k), label = f"Chi-Square fit for k = {k}")
    plt.title(r"Parameter $\frac{{(k - 1) \cdot s^2}}{{\sigma^2}}$ for k = {} and {} random sets".format(k, random_sets), 
              fontsize=16)
    plt.xlabel(r"$\frac{{(k - 1) \cdot s^2}}{{\sigma^2}}$", fontsize=16)
    plt.ylabel("Density", fontsize=16)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    #RESULTS FOR QUESTION (c):
    print("===========================================================================================================")
    print(f"MEAN  PERCENTAGE FOR k = {k} AND {random_sets} RANDOM SETS: ", len(accepted_mean ) / random_sets * 100, "%")
    print(f"SIGMA PERCENTAGE FOR k = {k} AND {random_sets} RANDOM SETS: ", len(accepted_sigma) / random_sets * 100, "%")

# EXTRA HISTOGRAM:
# Question b: Merged histograms:
plt.figure(figsize = (8, 6))

plt.hist(b_values[0], bins = 40,  color='salmon', edgecolor='black', label = f"k = {k_values[0]}", density=True)
plt.plot(np.linspace(min(b_values[0]), max(b_values[0]), 1000), 
            chi2.pdf(np.linspace(0, max(b_values[0]), 1000), 3), label = f"Chi-Square fit for k = 5", color = "black")

plt.hist(b_values[1], bins = 40,  color='blue',   edgecolor='black', label = f"k = {k_values[1]}", density=True)
plt.plot(np.linspace(min(b_values[1]), max(b_values[1]), 1000), 
            chi2.pdf(np.linspace(min(b_values[1]), max(b_values[1]), 1000), 20), label = f"Chi-Square fit for k = 20", color = "pink")

plt.hist(b_values[2], bins = 40,  color='green',  edgecolor='black', label = f"k = {k_values[2]}", density=True)
plt.plot(np.linspace(min(b_values[2]), max(b_values[2]), 1000), 
            chi2.pdf(np.linspace(min(b_values[2]), max(b_values[2]), 1000), 50), label = f"Chi-Square fit for k = 50", color = "cyan")

plt.title(r"Parameter $\frac{{(k - 1) \cdot s^2}}{{\sigma^2}}$ for all k values and {} random sets".format(random_sets), 
          fontsize=16)
plt.xlabel(r"$\frac{{(k - 1) \cdot s^2}}{{\sigma^2}}$", fontsize=16)
plt.ylabel("Density", fontsize=16)
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
