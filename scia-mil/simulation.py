import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist

def simulate_grid(grid_size, num_dots, spacing=15):
    locations = np.random.uniform(0, grid_size, size=(1,2))
    iter = 0
    while locations.shape[0]<num_dots:
        iter += 1
        if iter >= 100000:
            print(
                'Max iter reached, final dataset size is {}'.format(
                    locations.shape[0]
                    )
                )
            break
        tmp = np.random.uniform(0, grid_size, size=(1,2))
        dist = cdist(tmp, locations)[0]
        if (dist<=spacing).any():
            continue
        else:
            locations = np.append(locations, tmp, axis=0)
    return locations

def estimate_radius(v, center, r_min, r_max, shape = 'blob'):
    shape = {'blob':'blob', 'ring':'ring'}[shape]
    if shape == 'blob':
        r = np.sqrt(sum((v-center)**2))
        return r < r_max
    elif shape == 'ring':
        r = np.sqrt(sum((v-center)**2))
        return (r > r_min) & (r < r_max)

def simulate_pip(receivers, senders, locations, dist_cutoff = 30):
    pip = []
    senders = np.array(senders)
    for r in receivers:
        dist_to_senders = cdist(locations[[r]], locations[senders])[0]
        crit = dist_to_senders <= dist_cutoff
        pip.append(senders[crit])
    return np.array(pip,dtype=object)

def plot_pip(pip, receivers, locations):
    pip_all = np.unique([x for p in pip for x in p])
    receivers_all = [receivers[i] for i in range(len(receivers)) if pip[i].shape[0]>0]
    _ = plt.figure(figsize=(8,8))
    sns.scatterplot(
        x = locations[:,0], y = locations[:,1], color='gray', s=10,
        label = 'Others')
    sns.scatterplot(
        x = locations[receivers,0], y = locations[receivers,1], 
        color = 'g', s = 10, label = 'Out-of-range Receivers')
    sns.scatterplot(
        x = locations[receivers_all,0], y = locations[receivers_all,1], 
        color = 'g', s = 25, label = 'Receivers')
    sns.scatterplot(
        x = locations[senders,0], y = locations[senders,1], 
        color = 'r', s = 10, label = 'Non-PIP')
    sns.scatterplot(
        x = locations[pip_all,0], y = locations[pip_all,1], color = 'r', s=50,
        label = 'PIP')
    plt.legend(bbox_to_anchor = (1,0.5), loc = 'center left')
    for i in range(len(receivers)):
        if len(pip[i]) == 0:
            continue
        else:
            r = receivers[i]
            p = pip[i]
            for indiv_p in p:
                start = locations[indiv_p]
                end = locations[r]
                dx, dy = end-start
                plt.arrow(
                    *start, dx, dy,
                    head_width = 10,
                    length_includes_head=True)

def simulate_senders_receivers():
    np.random.seed(0)
    centers = np.array([[292, 799],[894, 149],[855, 774],[611, 391],[297, 366]])
    radius = np.random.randint(75, 150, size = 5)
    senders = []
    receivers = []
    for c, R in zip(centers,radius):
        _sender = [
            i for i in range(num_dots) if estimate_radius(
                locations[i], c, None, R)
                ]
        _receiver = [
            i for i in range(num_dots) if estimate_radius(
                locations[i], c, R,R+50,shape = 'ring')
                ]
        senders += _sender
        receivers += _receiver
    senders = list(set(senders))
    receivers = list(set(receivers))
    return senders, receivers

output_path = '/project/shared/xiao_wang/projects/cell2cell_inter/simulation'
grid_size = 1000
num_dots = 2500
dist_cutoff = 50
locations = simulate_grid(grid_size, num_dots)
senders, receivers = simulate_senders_receivers()
pip = simulate_pip(receivers, senders, locations, dist_cutoff)
total_ip = simulate_pip(receivers, senders, locations, dist_cutoff * 2)

_ = plt.figure(figsize=(8,8))
sns.scatterplot(x = locations[:,0], y = locations[:,1], color='gray')
sns.scatterplot(x = locations[senders,0], y = locations[senders,1], color = 'r')
sns.scatterplot(x = locations[receivers,0], y = locations[receivers,1], color = 'g')

pip_n = [len(x) for x in pip]
plt.hist(pip_n)
plot_pip(total_ip, receivers, locations)

# simulate sender expression with normaal distribution
# This also means in the real application, the pathways/gene expression need to be
# z-scored
exp_sender = pd.DataFrame(
    np.random.normal(size=[len(senders), 50]), index = senders)

# simulate betas
primary_beta = np.random.normal(0, 0.1, size=10)
trivial_beta = np.random.uniform(0 , 0.0001, 40)
betas = np.append(primary_beta, trivial_beta)

exp_receiver = np.zeros(shape=(len(receivers), 1))
for i in range(len(receivers)):
    i_senders = pip[i]
    if len(i_senders) == 0:
        i_senders = total_ip[i]
    exp_receiver[i] = np.matmul(exp_sender.loc[i_senders], betas).sum()

exp_receiver = (exp_receiver > 0) + 0

pd.DataFrame(exp_receiver).to_csv(
    output_path + '/exp_receiver', index=None,header=None)
# exp_sender needs to have index as to matrix send to spacia.r
exp_sender.to_csv(
    output_path + '/exp_sender', index=None,header=None)
pd.DataFrame(exp_receiver).to_csv(
    output_path + '/exp_receiver', index=None,header=None)