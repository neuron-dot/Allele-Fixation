# Allele Fixation 
import matplotlib.pyplot as plt

def simulate_selection(p_init, s, dominant, max_gen=30000):
    p = p_init
    generations = []
    frequencies = []
    threshold_gen = None

    for gen in range(max_gen):
        generations.append(gen)
        frequencies.append(p)
        if p < 0.05 and threshold_gen is None:
            threshold_gen = gen
            print(f'Threshold reached at generation {threshold_gen}, frequency {p:.5f}')
        if p < 0.001:
            break
        if dominant:
            w_bar = 1 - s * (p**2 + 2*p*(1-p))
            p = (p*(1-s)) / w_bar
        else:
            w_bar = 1 - s * (p**2)
            p = ((p**2)*(1-s) + p*(1-p)) / w_bar
    return generations, frequencies, threshold_gen

# Initial conditions
p_init = 0.99
strong_s = 0.25
weak_s = 0.008

# Simulations
scenarios = [
    (simulate_selection(p_init, strong_s, False), 'Recessive - Strong Selection'),
    (simulate_selection(p_init, weak_s, False), 'Recessive - Weak Selection'),
    (simulate_selection(p_init, strong_s, True), 'Dominant - Strong Selection'),
    (simulate_selection(p_init, weak_s, True), 'Dominant - Weak Selection')
]

# Separate Plots
for (gens, freqs, th_gen), title in scenarios:
    plt.figure(figsize=(12, 6))
    plt.plot(gens, freqs, label=title)
    plt.axhline(y=0.05, color='r', linestyle='--', label='Threshold (0.05)')
    if th_gen is not None:
        plt.scatter(th_gen, 0.05, color='red', label=f'p<0.05 at generation {th_gen}')
    plt.xlabel('Generations')
    plt.ylabel('Allele Frequency (p)')
    plt.title(title)
    plt.legend()
    plt.show()
